# core_functions.R -- Analysis engine for Epitope-Evaluator 2
#
# Everything in this file is a pure function of an ee_dataset plus parameters:
# no Shiny, no reactivity, no plotting. That is what makes tests/test_core.R
# possible without installing Shiny, and it is why the six tools share one
# implementation of the set logic instead of six near-copies.
#
# Complexity summary (n = peptide rows, k = alleles, c = cutoffs, p = proteins):
#
#   operation                v1                              v2
#   ------------------------------------------------------------------------
#   combine alleles (AND/OR) apply() over a data.frame       matrixStats rowMins/rowMaxs
#   distribution histogram   plotly binning in the browser   findInterval + tabulate, O(n)
#   allele x cutoff grid     O(k*c) full-column scans        one sort per allele, O(n log n + k*c)
#   UpSet intersections      O(2^k) full-table filters       one grouping pass, O(n*k)
#   protein x allele counts  melt() to n*k long rows         tabulate per column, O(n*k), no copy
#   Viewer lane packing      O(n^2 * L) pairwise overlap     first-fit sweep, O(n * lanes)

# ---------------------------------------------------------------------------
# Score combination
# ---------------------------------------------------------------------------

#' Combine several alleles into one score per peptide
#'
#' "Intersection" (AND) asks for peptides that bind *every* selected allele, so
#' the governing value is the worst (max) rank. "Union" (OR) asks for peptides
#' binding *at least one*, so it is the best (min) rank. NAs -- which MHCflurry
#' and IEDB produce for allele/peptide pairs that were never predicted -- are
#' excluded rather than poisoning the whole row, which is what v1 did.
#'
#' @param ds An ee_dataset.
#' @param alleles Character vector of allele names.
#' @param mode "Intersection" or "Union".
#' @param unique_only If TRUE, return one value per distinct peptide.
#' @return Numeric vector; NA where the peptide has no value for any selected allele.
ee_combine <- function(ds, alleles, mode = c("Intersection", "Union"), unique_only = FALSE) {
  mode <- match.arg(mode)
  alleles <- ee_check_alleles(ds, alleles)
  rows <- if (unique_only) ds$uniq_idx else seq_len(ds$n_rows)
  m <- ds$scores[rows, alleles, drop = FALSE]

  if (length(alleles) == 1L) return(as.numeric(m[, 1L]))
  out <- if (mode == "Intersection") matrixStats::rowMaxs(m, na.rm = TRUE)
         else                        matrixStats::rowMins(m, na.rm = TRUE)
  # rowMins/rowMaxs return +/-Inf for all-NA rows.
  out[!is.finite(out)] <- NA_real_
  out
}

#' Build a peptide x allele "is a binder" matrix, one column at a time
#'
#' `ds$scores[rows, alleles] <= cutoff` is the obvious spelling, but it
#' materialises a full numeric copy of the slice AND a full logical matrix
#' before either can be freed. At 10^6 peptides x 20 alleles that is ~240 MB of
#' simultaneous temporaries, which is the difference between fitting and not
#' fitting in a 1 GB container. Filling the result column by column keeps the
#' transient cost at one column.
#'
#' NA (a pair the predictor never scored) counts as "not a binder", never as 0.
ee_binder_matrix <- function(ds, rows, alleles, cutoff) {
  full <- ee_is_all_rows(ds, rows)
  out <- matrix(FALSE, length(rows), length(alleles), dimnames = list(NULL, alleles))
  for (j in seq_along(alleles)) {
    v <- if (full) ds$scores[, alleles[j]] else ds$scores[rows, alleles[j]]
    out[, j] <- !is.na(v) & v <= cutoff
  }
  out
}

#' Number of selected alleles each row binds, without building a matrix at all
ee_binder_count <- function(ds, rows, alleles, cutoff) {
  full <- ee_is_all_rows(ds, rows)
  n <- integer(length(rows))
  for (a in alleles) {
    v <- if (full) ds$scores[, a] else ds$scores[rows, a]
    n <- n + (!is.na(v) & v <= cutoff)
  }
  as.integer(n)
}

#' Does this row index select every row, in order?
#'
#' When every peptide is unique -- the common case for a whole-proteome run --
#' uniq_idx is just seq_len(n), and skipping the fancy indexing avoids copying
#' a column per allele.
ee_is_all_rows <- function(ds, rows) {
  length(rows) == nrow(ds$scores) && rows[1L] == 1L &&
    rows[length(rows)] == nrow(ds$scores) && !is.unsorted(rows, strictly = TRUE)
}

#' Validate an allele selection and fail with a useful message
ee_check_alleles <- function(ds, alleles) {
  if (!length(alleles)) {
    stop("Select at least one MHC allele.", call. = FALSE)
  }
  unknown <- setdiff(alleles, ds$alleles)
  if (length(unknown)) {
    stop("Allele(s) not present in the loaded data: ", ee_truncate_list(unknown), ".", call. = FALSE)
  }
  alleles
}

ee_check_proteins <- function(ds, proteins) {
  if (!length(proteins)) stop("Select at least one protein.", call. = FALSE)
  unknown <- setdiff(proteins, ds$proteins$ID)
  if (length(unknown)) {
    stop("Protein(s) not present in the loaded data: ", ee_truncate_list(unknown), ".", call. = FALSE)
  }
  proteins
}

#' Sanity-check a cutoff range
ee_check_range <- function(xmin, xmax, step = NULL) {
  if (!is.finite(xmin) || !is.finite(xmax)) stop("The %rank limits must be numbers.", call. = FALSE)
  if (xmax <= xmin) stop(sprintf("Max %%rank (%g) must be greater than min %%rank (%g).", xmax, xmin),
                         call. = FALSE)
  if (!is.null(step)) {
    if (!is.finite(step) || step <= 0) stop("Bin width must be a positive number.", call. = FALSE)
    if ((xmax - xmin) / step > 5000) {
      stop(sprintf("Bin width %g would create %.0f bins over [%g, %g]. Increase the bin width.",
                   step, (xmax - xmin) / step, xmin, xmax), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Tool 1 -- Distribution
# ---------------------------------------------------------------------------

#' Bin scores into a histogram server-side
#'
#' v1 handed every raw value to plotly and let the browser bin it, which means
#' shipping one number per peptide per redraw. Binning here sends a few hundred
#' numbers instead, and is what makes the tool usable at 10^6 peptides.
#'
#' @return data.table(bin_start, bin_end, bin_mid, count, cumulative, density)
ee_histogram <- function(values, xmin, xmax, step, cumulative = FALSE) {
  ee_check_range(xmin, xmax, step)
  values <- values[is.finite(values)]

  breaks <- seq(xmin, xmax, by = step)
  if (breaks[length(breaks)] < xmax - 1e-9) breaks <- c(breaks, breaks[length(breaks)] + step)
  nb <- length(breaks) - 1L
  if (nb < 1L) stop("The chosen range and bin width produce no bins.", call. = FALSE)

  # Half-open bins [lo, hi), with the last bin closed so xmax itself is counted.
  keep <- values >= xmin & values <= xmax
  idx  <- findInterval(values[keep], breaks, rightmost.closed = TRUE, all.inside = TRUE)
  cnt  <- tabulate(idx, nbins = nb)

  total <- sum(cnt)
  data.table(
    bin_start  = breaks[-length(breaks)],
    bin_end    = breaks[-1L],
    bin_mid    = (breaks[-length(breaks)] + breaks[-1L]) / 2,
    count      = cnt,
    cumulative = cumsum(cnt),
    density    = if (total > 0) cnt / (total * step) else rep(0, nb)
  )
}

#' Epitopes passing a filter, as a downloadable table
#'
#' @return data.table(Peptide, Pos, End, ID, PepLength, Score) sorted by score.
ee_distribution_table <- function(ds, alleles, mode, xmin, xmax) {
  ee_check_range(xmin, xmax)
  sc   <- ee_combine(ds, alleles, mode)
  keep <- which(is.finite(sc) & sc >= xmin & sc <= xmax)
  out  <- ds$peptides[keep, .(Peptide, Pos, End, PepLength, ID)]
  out[, Score := round(sc[keep], 4)]
  data.table::setorder(out, Score, ID, Pos)
  out[]
}

#' Number of epitopes per allele at a ladder of cutoffs
#'
#' v1 looped over cutoffs and rescanned every allele column each time
#' (100 cutoffs x 12 alleles = 1200 full-column comparisons). Sorting each
#' column once and using findInterval turns that into one sort per allele.
#'
#' @return data.table(allele, cutoff, count) in long form, ready for a heatmap.
ee_cutoff_grid <- function(ds, alleles, cutoffs, unique_only = TRUE) {
  alleles <- ee_check_alleles(ds, alleles)
  cutoffs <- sort(unique(cutoffs[is.finite(cutoffs)]))
  if (!length(cutoffs)) stop("No valid cutoff values.", call. = FALSE)

  rows <- if (unique_only) ds$uniq_idx else seq_len(ds$n_rows)
  res <- vector("list", length(alleles))
  for (j in seq_along(alleles)) {
    v <- ds$scores[rows, alleles[j]]
    v <- sort(v[is.finite(v)])
    res[[j]] <- data.table(
      allele = alleles[j],
      cutoff = cutoffs,
      count  = findInterval(cutoffs, v)   # number of values <= each cutoff
    )
  }
  out <- data.table::rbindlist(res)
  out[, allele := factor(allele, levels = alleles)]
  out[]
}

# ---------------------------------------------------------------------------
# Tools 2 & 6 -- Intersection (over alleles) and Conservation (over proteins)
#
# Both answer "which elements belong to which combination of sets"; only the
# definition of a set differs. v1 implemented them twice, and both times by
# enumerating all 2^k combinations and running a full-table filter per
# combination -- 1,023 table scans for 10 sets, using dplyr::filter_(), which
# has been defunct since dplyr 1.1 and so crashes outright on current R.
# ---------------------------------------------------------------------------

#' Build a peptide x set membership matrix
#'
#' @param by "allele": set j = peptides binding allele j at <= cutoff.
#'           "protein": set j = peptides occurring in protein j, restricted to
#'           those that pass the allele filter.
#' @return list(m = logical matrix, peptides = character vector of row names)
ee_membership <- function(ds, sets, cutoff, by = c("allele", "protein"),
                          alleles = NULL, mode = "Union") {
  by <- match.arg(by)
  if (!is.finite(cutoff)) stop("The cutoff must be a number.", call. = FALSE)
  if (length(sets) < 2L) {
    stop("Select at least two ", by, "s to compare.", call. = FALSE)
  }

  if (by == "allele") {
    sets <- ee_check_alleles(ds, sets)
    m <- ee_binder_matrix(ds, ds$uniq_idx, sets, cutoff)
    peps <- ds$uniq_pep
  } else {
    sets <- ee_check_proteins(ds, sets)
    sc   <- ee_combine(ds, alleles, mode)               # per row, not per peptide
    ok   <- which(is.finite(sc) & sc <= cutoff & ds$peptides$ID %chin% sets)
    if (!length(ok)) {
      return(list(m = matrix(FALSE, 0L, length(sets), dimnames = list(NULL, sets)),
                  peptides = character(0)))
    }
    pid  <- ds$pep_index[ok]
    keep <- sort(unique(pid))
    pos  <- match(pid, keep)
    m <- matrix(FALSE, length(keep), length(sets), dimnames = list(NULL, sets))
    m[cbind(pos, match(ds$peptides$ID[ok], sets))] <- TRUE
    peps <- ds$uniq_pep[keep]
  }

  keep <- which(rowSums(m) > 0L)
  list(m = m[keep, , drop = FALSE], peptides = peps[keep])
}

#' Collapse a membership matrix into observed set combinations
#'
#' One grouping pass over the rows. Only combinations that actually occur are
#' returned, which is also what makes the result usable for 20+ sets where
#' enumerating 2^k is impossible.
#'
#' @param max_peptides Cap on how many peptide sequences are pasted into the
#'   `peptides` column. The on-screen table shortens this further via
#'   ee_peptide_preview(); the download always writes this full string.
#' @return data.table(n_sets, sets, count, peptides) sorted by count descending.
ee_combinations <- function(mem, collect_peptides = TRUE, max_peptides = 2000L) {
  m <- mem$m
  if (!nrow(m)) {
    return(data.table(n_sets = integer(0), sets = character(0),
                      count = integer(0), peptides = character(0)))
  }
  cols <- colnames(m)
  dt <- as.data.table(m)
  dt[, .pep := mem$peptides]

  agg <- if (collect_peptides) {
    dt[, .(count = .N,
           peptides = paste(utils::head(.pep, max_peptides), collapse = ", ")),
       by = cols]
  } else {
    dt[, .(count = .N), by = cols]
  }

  flags <- as.matrix(agg[, cols, with = FALSE])
  agg[, n_sets := as.integer(rowSums(flags))]
  agg[, sets := apply(flags, 1L, function(r) paste(cols[r], collapse = " & "))]
  keep <- c("n_sets", "sets", "count", if (collect_peptides) "peptides")
  out <- agg[n_sets > 0L, ..keep]
  data.table::setorder(out, -count, -n_sets)
  out[]
}

#' Set sizes (the horizontal bars of an UpSet plot)
ee_set_sizes <- function(mem) {
  data.table(set = colnames(mem$m), size = as.integer(colSums(mem$m)))
}

# ---------------------------------------------------------------------------
# Tool 3 -- Density
# ---------------------------------------------------------------------------

#' Epitopes per protein and epitope density
#'
#' Behaviour change vs v1: ProtLength is the true amino-acid length from the
#' FASTA. v1 used nchar(sequence) - 8 for every predictor, so class I densities
#' were divided by L-8 and class II densities by L-8 as well (the correct window
#' count for 15-mers is L-14). Densities here are epitopes per amino acid.
#'
#' @return data.table(ID, ProtLength, Epitopes, Density) for every protein,
#'   including proteins with zero epitopes.
ee_density_table <- function(ds, alleles, mode, cutoff) {
  if (!is.finite(cutoff)) stop("The cutoff must be a number.", call. = FALSE)
  sc   <- ee_combine(ds, alleles, mode)
  hit  <- is.finite(sc) & sc <= cutoff

  cnt <- ds$peptides[hit, .(Epitopes = .N), by = ID]
  out <- merge(ds$proteins, cnt, by = "ID", all.x = TRUE)
  out[is.na(Epitopes), Epitopes := 0L]
  out[, Density := ifelse(ProtLength > 0, Epitopes / ProtLength, NA_real_)]
  data.table::setorder(out, -Epitopes)
  out[]
}

#' Epitope counts for every protein x allele pair
#'
#' v1 melted the whole table to n*k long rows through reshape::melt before
#' grouping. Here each allele column is counted with tabulate() directly against
#' the protein index, so nothing is reshaped and nothing is copied.
#'
#' @return data.table(ID, allele, Epitopes, ProtLength, Density) -- a complete
#'   grid, so zero-epitope combinations render as empty cells rather than holes.
ee_protein_allele_counts <- function(ds, alleles, cutoff) {
  alleles <- ee_check_alleles(ds, alleles)
  if (!is.finite(cutoff)) stop("The cutoff must be a number.", call. = FALSE)

  prot <- ds$proteins
  pidx <- match(ds$peptides$ID, prot$ID)
  np   <- nrow(prot)

  counts <- vapply(alleles, function(a) {
    v <- ds$scores[, a]
    as.numeric(tabulate(pidx[is.finite(v) & v <= cutoff], nbins = np))
  }, numeric(np))
  dim(counts) <- c(np, length(alleles))

  out <- data.table(
    ID       = rep(prot$ID, times = length(alleles)),
    allele   = rep(alleles, each = np),
    Epitopes = as.vector(counts)
  )
  out[, ProtLength := rep(prot$ProtLength, times = length(alleles))]
  out[, Density := ifelse(ProtLength > 0, Epitopes / ProtLength, NA_real_)]
  out[, allele := factor(allele, levels = alleles)]
  out[]
}

# ---------------------------------------------------------------------------
# Tool 4 -- Viewer
# ---------------------------------------------------------------------------

#' Assign non-overlapping display lanes to intervals
#'
#' First-fit sweep in start order: an interval goes into the first lane whose
#' last occupied position ends before it starts. This is the classic greedy
#' interval-graph colouring, so it uses the minimum possible number of lanes.
#'
#' v1 did this with a nested while/for over data.frame rows, materialising
#' seq(start, end) for both members of every pair and intersecting them -- which
#' is why the Viewer was documented as taking "a few minutes". Lane depth is
#' bounded by the peptide length (windows tile the protein one residue apart),
#' so the inner scan is over ~9-15 lanes regardless of protein size.
#'
#' @param start,end Integer vectors of equal length (inclusive coordinates).
#' @param gap Minimum empty residues required between two intervals in a lane.
#' @return Integer vector of 1-based lane numbers, aligned to the input order.
ee_pack_lanes <- function(start, end, gap = 1L) {
  n <- length(start)
  if (n == 0L) return(integer(0))
  ord <- order(start, end)
  lane_end <- numeric(0)
  lane <- integer(n)

  for (i in ord) {
    s <- start[i]
    placed <- FALSE
    if (length(lane_end)) {
      free <- which(lane_end + gap < s)
      if (length(free)) {
        L <- free[1L]
        lane[i] <- L
        lane_end[L] <- end[i]
        placed <- TRUE
      }
    }
    if (!placed) {
      lane_end <- c(lane_end, end[i])
      lane[i] <- length(lane_end)
    }
  }
  lane
}

#' Epitope layout for one protein
#'
#' @return list(protein = list(id, length), epitopes = data.table(...), n_lanes)
ee_viewer_layout <- function(ds, protein, alleles, mode, cutoff,
                             max_rects = EE$limits$viewer_rects) {
  protein <- ee_check_proteins(ds, protein)[1]
  alleles <- ee_check_alleles(ds, alleles)
  if (!is.finite(cutoff)) stop("The cutoff must be a number.", call. = FALSE)

  in_prot <- which(ds$peptides$ID == protein)
  plen <- ds$proteins$ProtLength[match(protein, ds$proteins$ID)]
  empty <- data.table(Peptide = character(0), Pos = integer(0), End = integer(0),
                      NAlleles = integer(0), Alleles = character(0), lane = integer(0))

  if (!length(in_prot)) return(list(protein = list(id = protein, length = plen),
                                    epitopes = empty, n_lanes = 0L, truncated = FALSE))

  bind <- ee_binder_matrix(ds, in_prot, alleles, cutoff)

  # An epitope must satisfy the AND/OR condition over the selected alleles.
  pass <- if (mode == "Intersection") rowSums(bind) == length(alleles) else rowSums(bind) > 0L
  sel <- which(pass)
  if (!length(sel)) return(list(protein = list(id = protein, length = plen),
                                epitopes = empty, n_lanes = 0L, truncated = FALSE))

  truncated <- FALSE
  if (length(sel) > max_rects) {
    # Keep the most promiscuous epitopes; drawing more than this locks the browser.
    sel <- sel[order(-rowSums(bind[sel, , drop = FALSE]))][seq_len(max_rects)]
    sel <- sort(sel)
    truncated <- TRUE
  }

  rows <- in_prot[sel]
  bsel <- bind[sel, , drop = FALSE]
  ep <- data.table(
    Peptide  = ds$peptides$Peptide[rows],
    Pos      = ds$peptides$Pos[rows],
    End      = ds$peptides$End[rows],
    NAlleles = as.integer(rowSums(bsel)),
    Alleles  = apply(bsel, 1L, function(r) paste(alleles[r], collapse = ", "))
  )
  ep[, lane := ee_pack_lanes(Pos, End)]
  data.table::setorder(ep, Pos)

  list(protein = list(id = protein, length = plen),
       epitopes = ep[], n_lanes = max(ep$lane), truncated = truncated)
}

# ---------------------------------------------------------------------------
# Tool 5 -- Promiscuity
# ---------------------------------------------------------------------------

#' Epitopes binding at least a minimum number of alleles
#'
#' Behaviour change vs v1: promiscuity is counted with <= against the weak
#' cutoff, matching both the documentation and the SB/WB colouring in the
#' heatmap. v1's table used a strict < while its plot used <=, so the two
#' disagreed for every epitope sitting exactly on the cutoff.
#'
#' @return list(table = data.table(Peptide, Pos, End, ID, NAlleles),
#'              matrix = character matrix of "SB"/"WB"/"" per peptide x allele,
#'              values = the underlying numeric scores)
ee_promiscuity <- function(ds, alleles, strong, weak, min_alleles,
                           max_rows = EE$limits$heatmap_rows) {
  alleles <- ee_check_alleles(ds, alleles)
  if (!is.finite(strong) || !is.finite(weak)) stop("Both cutoffs must be numbers.", call. = FALSE)
  if (strong > weak) {
    stop(sprintf("The strong-binding cutoff (%g) must not exceed the weak-binding cutoff (%g).",
                 strong, weak), call. = FALSE)
  }
  if (!is.finite(min_alleles) || min_alleles < 1) min_alleles <- 1

  u <- ds$uniq_idx
  # Counted without materialising the peptide x allele matrix: only the rows
  # that survive the filter are sliced out below, and there are usually few.
  prom <- ee_binder_count(ds, u, alleles, weak)

  sel <- which(prom >= min_alleles)
  if (!length(sel)) {
    return(list(table = data.table(Peptide = character(0), Pos = integer(0), End = integer(0),
                                   ID = character(0), NAlleles = integer(0)),
                matrix = matrix(character(0), 0L, length(alleles), dimnames = list(NULL, alleles)),
                values = matrix(numeric(0), 0L, length(alleles), dimnames = list(NULL, alleles)),
                truncated = FALSE, n_total = 0L))
  }
  sel <- sel[order(-prom[sel])]

  tab <- ds$peptides[u[sel], .(Peptide, Pos, End, ID)]
  tab[, NAlleles := prom[sel]]

  truncated <- length(sel) > max_rows
  show <- if (truncated) sel[seq_len(max_rows)] else sel

  # Only the rows that will actually be drawn are sliced out of the score
  # matrix -- at most EE$limits$heatmap_rows of them, however large the input.
  vals <- ds$scores[u[show], alleles, drop = FALSE]
  dimnames(vals) <- list(NULL, alleles)
  cls  <- matrix("", nrow(vals), ncol(vals), dimnames = list(ds$uniq_pep[show], alleles))
  cls[!is.na(vals) & vals <= weak & vals > strong] <- "WB"
  cls[!is.na(vals) & vals <= strong] <- "SB"

  list(table = tab[], matrix = cls, values = vals,
       truncated = truncated, n_total = length(sel))
}

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

#' Suggested defaults for the current dataset
ee_suggest <- function(ds) {
  d <- EE$defaults[[if (ds$mhc_class == "II") "II" else "I"]]
  if (ds$score_type == "Score") {
    # Binding-affinity scores are not percentile ranks; 500 nM is the classic
    # class I threshold and 1000 nM the usual class II one.
    d <- list(cutoff = if (ds$mhc_class == "II") 1000 else 500,
              strong = if (ds$mhc_class == "II") 500 else 50,
              weak   = if (ds$mhc_class == "II") 1000 else 500,
              step   = if (ds$mhc_class == "II") 50 else 25)
  }
  d$xmin <- 0
  d$xmax <- d$cutoff
  d$score_max <- suppressWarnings(max(ds$scores, na.rm = TRUE))
  if (!is.finite(d$score_max)) d$score_max <- d$cutoff
  d
}

#' Choose a "minimum number of alleles" that actually returns something
#'
#' v1 defaulted to (number of alleles - 1), inherited here at first. On real
#' data almost nothing binds 11 of 12 alleles at 2 %rank, so the tool opened
#' with an empty heatmap and the user had to guess a workable number.
#'
#' Instead, count how many alleles each peptide binds and pick the strictest
#' threshold that still leaves at least `target` epitopes -- so the tool opens
#' on the most promiscuous set worth looking at, whatever the dataset.
#'
#' @return An integer between 2 and the number of alleles.
ee_suggest_min_alleles <- function(ds, alleles, weak, target = 25L) {
  alleles <- ee_check_alleles(ds, alleles)
  k <- length(alleles)
  if (k < 2L) return(1L)
  prom <- ee_binder_count(ds, ds$uniq_idx, alleles, weak)
  # counts[i] = how many peptides bind at least i alleles
  at_least <- rev(cumsum(rev(tabulate(prom, nbins = k))))
  ok <- which(at_least >= target)
  if (!length(ok)) return(2L)
  max(2L, min(max(ok), k))
}

#' Human-readable label for the score axis
ee_score_label <- function(ds) {
  if (identical(ds$score_type, "Rank")) "% rank" else
    if (identical(ds$predictor, "NetMHC") || identical(ds$predictor, "IEDB Consensus"))
      "Binding affinity (nM)" else "Binding score"
}

#' Shorten a comma-separated peptide list for on-screen display
#'
#' A single intersection can hold thousands of peptides. Pasting all of them
#' into one cell makes the row thousands of pixels tall and the table unusable,
#' so the screen version shows the first few and says how many are hidden. The
#' download button still writes the complete list.
ee_peptide_preview <- function(x, n = 6L) {
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return("")
    v <- strsplit(s, ", ", fixed = TRUE)[[1]]
    if (length(v) <= n) return(paste(v, collapse = ", "))
    sprintf("%s  (+%s more)", paste(v[seq_len(n)], collapse = ", "),
            format(length(v) - n, big.mark = ","))
  }, character(1), USE.NAMES = FALSE)
}
