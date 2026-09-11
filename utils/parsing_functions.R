# parsing_functions.R -- Input readers for Epitope-Evaluator 2
#
# Every predictor output is normalised into one canonical object (see
# ee_dataset() at the bottom of this file):
#
#   peptides : data.table(Peptide, Pos, End, ID, ProtLength)   one row per (peptide, protein)
#   scores   : numeric matrix, nrow(peptides) x n_alleles      one column per MHC allele
#
# Design notes vs v1:
#   * data.table::fread everywhere (v1 used read.table / read.delim), 10-40x faster
#     and it does not coerce strings to factors or mangle allele names.
#   * The NetMHC family is parsed from its two header rows rather than from
#     hard-coded column strides. v1 assumed a stride of 3, 4 or 6 depending on a
#     fragile ncol arithmetic check and picked the wrong columns for
#     NetMHCpan/NetMHCIIpan "Score" when the -BA block was absent.
#   * Protein IDs are matched to the FASTA by name (with truncation-aware prefix
#     matching) and cross-checked against the length implied by the prediction
#     file. v1 matched by position in a unique() vector, which silently shifted
#     every ID whenever a protein was too short to produce a peptide.
#   * Scores are stored as one numeric matrix, not as data.frame columns, so the
#     row-wise min/max that every tool needs is a matrixStats call over
#     contiguous memory instead of apply() over a list of columns.

# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------

#' Read a multi-FASTA file
#'
#' Replaces seqinr::read.fasta + phylotools::get.fasta.name (two dependencies,
#' three passes over the file) with a single readLines pass.
#'
#' @param path Path to a FASTA file.
#' @return data.table(Name, Header, Length, Sequence); errors with a readable
#'   message if the file is not FASTA.
ee_read_fasta <- function(path) {
  if (!file.exists(path)) stop("FASTA file not found: ", path, call. = FALSE)

  lines <- tryCatch(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    error = function(e) stop("Could not read the FASTA file: ", conditionMessage(e), call. = FALSE)
  )
  lines <- lines[nzchar(stri_trim_both(lines))]
  if (!length(lines)) stop("The FASTA file is empty.", call. = FALSE)

  hdr_idx <- which(stri_startswith_fixed(lines, ">"))
  if (!length(hdr_idx)) {
    stop("No '>' header lines found: this does not look like a FASTA file.", call. = FALSE)
  }
  if (hdr_idx[1] != 1L) {
    stop("The FASTA file must start with a '>' header line.", call. = FALSE)
  }

  headers <- stri_sub(lines[hdr_idx], 2L)
  # Sequence block for record i runs from hdr_idx[i]+1 to hdr_idx[i+1]-1.
  block  <- rep.int(seq_along(hdr_idx), diff(c(hdr_idx, length(lines) + 1L)))
  is_seq <- !seq_along(lines) %in% hdr_idx
  seqs   <- vapply(
    split(lines[is_seq], block[is_seq]),
    function(x) stri_flatten(x),
    character(1)
  )
  # split() drops records with no sequence lines; re-align by name.
  full <- character(length(hdr_idx))
  full[as.integer(names(seqs))] <- unname(seqs)
  full <- stri_replace_all_charclass(stri_trans_toupper(full), "[\\s*-]", "")

  if (!any(nzchar(full))) {
    stop("Every FASTA record is empty: no sequences to analyse.", call. = FALSE)
  }

  out <- data.table(
    Name     = stri_split_regex(headers, "\\s+", n = 2L, simplify = TRUE)[, 1L],
    Header   = headers,
    Length   = stri_length(full),
    Sequence = full
  )
  if (anyDuplicated(out$Name)) {
    dup <- unique(out$Name[duplicated(out$Name)])
    warning(sprintf("Duplicated FASTA identifiers (%s). Keeping the first occurrence of each.",
                    ee_truncate_list(dup)), call. = FALSE)
    out <- out[!duplicated(Name)]
  }
  out[]
}

# ---------------------------------------------------------------------------
# Protein-ID reconciliation
# ---------------------------------------------------------------------------

#' Normalise an identifier for comparison
#'
#' Predictors rewrite FASTA names: NetMHC turns "sp|P0DTC2|SPIKE_SARS2" into
#' "sp_P0DTC2_SPIKE" (pipes to underscores, truncated to a fixed width).
ee_norm_id <- function(x) {
  x <- stri_replace_all_charclass(as.character(x), "[|:;,/\\\\ ]", "_")
  stri_trans_toupper(stri_replace_all_regex(x, "_+", "_"))
}

#' Candidate lookup keys for one FASTA record
#'
#' Predictors derive their protein column from the FASTA header in several
#' incompatible ways, then truncate the result to a fixed width. For
#' "tr|A0A663DJA2|A0A663DJA2_SARS2" the observed variants are the whole name,
#' the name without the database prefix, and the bare entry name -- so every one
#' of them is registered as a key and prefix matching is tried against all.
ee_id_keys <- function(names_vec) {
  n  <- ee_norm_id(names_vec)
  bits <- stri_split_fixed(n, "_")
  lapply(seq_along(n), function(i) {
    b <- bits[[i]]
    keys <- n[i]
    if (length(b) > 1L) keys <- c(keys, stri_flatten(b[-1], "_"))   # drop db prefix
    if (length(b) > 2L) keys <- c(keys, stri_flatten(b[-(1:2)], "_")) # bare entry name
    keys <- c(keys, b[nzchar(b)])                                    # individual fields
    unique(keys[nzchar(keys)])
  })
}

#' Map prediction-file protein IDs onto FASTA records
#'
#' Three passes, each only filling in what the previous one left unresolved:
#'   1. exact match against any lookup key of any record
#'   2. the prediction ID is a prefix of a key (handles fixed-width truncation);
#'      ties are broken by the protein length implied by the prediction file
#'   3. length-only match against the still-unused records
#' A positional fallback is used last, and only when the counts agree exactly.
#'
#' The implied length is max(start + peptide length - 1) over the protein:
#' predictors emit every sliding window, so this equals the true protein length
#' and makes the mapping self-checking. v1 had no such check and silently
#' misassigned every ID when a protein was too short to yield a peptide.
#'
#' @param pred_ids Unique protein IDs, in order of first appearance.
#' @param fasta data.table from ee_read_fasta().
#' @param implied_len Named numeric: implied protein length per pred_id, or NULL.
#' @return list(name, length, notes), each aligned to pred_ids.
ee_match_proteins <- function(pred_ids, fasta, implied_len = NULL) {
  pred_ids <- as.character(pred_ids)
  np <- ee_norm_id(pred_ids)
  keys <- ee_id_keys(fasta$Name)

  # Flatten to a key -> record lookup table.
  owner <- rep.int(seq_len(nrow(fasta)), lengths(keys))
  flat  <- unlist(keys, use.names = FALSE)

  il <- if (!is.null(implied_len)) as.numeric(implied_len[pred_ids]) else rep(NA_real_, length(np))
  idx  <- rep(NA_integer_, length(np))
  how  <- rep(NA_character_, length(np))
  notes <- character(0)

  take <- function(cands, i) {
    cands <- unique(cands[!cands %in% idx[-i]])          # a record maps to one ID only
    if (length(cands) == 1L) return(cands)
    if (length(cands) > 1L && !is.na(il[i])) {           # break ties on protein length
      byl <- cands[fasta$Length[cands] == il[i]]
      if (length(byl) == 1L) return(byl)
    }
    NA_integer_
  }

  # Pass 1 -- exact.
  for (i in seq_along(np)) {
    hit <- take(owner[flat == np[i]], i)
    if (!is.na(hit)) { idx[i] <- hit; how[i] <- "exact name match" }
  }
  # Pass 2 -- truncation-aware prefix.
  for (i in which(is.na(idx))) {
    hit <- take(owner[stri_startswith_fixed(flat, np[i])], i)
    if (!is.na(hit)) { idx[i] <- hit; how[i] <- "truncated-name prefix match" }
  }
  # Pass 3 -- unique length among the records nothing has claimed.
  free <- setdiff(seq_len(nrow(fasta)), idx[!is.na(idx)])
  for (i in which(is.na(idx))) {
    if (is.na(il[i])) next
    hit <- free[fasta$Length[free] == il[i]]
    if (length(hit) == 1L) {
      idx[i] <- hit; how[i] <- "protein-length match"
      free <- setdiff(free, hit)
    }
  }
  # Pass 4 -- positional, only if nothing at all resolved and the counts agree.
  if (all(is.na(idx)) && length(np) == nrow(fasta)) {
    idx <- seq_along(np); how[] <- "positional fallback"
    notes <- c(notes, paste(
      "Protein identifiers could not be matched to the FASTA by name, so they were matched",
      "by order of appearance. Verify that this FASTA is the file submitted to the predictor."
    ))
  }

  if (anyNA(idx)) {
    stop(sprintf(
      paste("Could not match %d of the %d protein(s) in the prediction file to the %d FASTA",
            "record(s). Unmatched: %s. Is this the FASTA used for the prediction?"),
      sum(is.na(idx)), length(np), nrow(fasta), ee_truncate_list(pred_ids[is.na(idx)])
    ), call. = FALSE)
  }

  name <- fasta$Name[idx]
  len  <- fasta$Length[idx]

  bad <- which(!is.na(il) & il != len)
  if (length(bad)) {
    notes <- c(notes, sprintf(
      paste("%d protein(s) have a FASTA length that disagrees with the prediction file",
            "(e.g. %s: FASTA %d aa vs %d aa implied by the last peptide position).",
            "Lengths were taken from the FASTA."),
      length(bad), name[bad[1]], len[bad[1]], il[bad[1]]))
  }
  nonexact <- setdiff(unique(how), c("exact name match", "positional fallback"))
  if (length(nonexact)) {
    notes <- c(notes, paste0("Protein IDs resolved by ", paste(nonexact, collapse = " and "), "."))
  }

  list(name = name, length = len, notes = notes)
}

ee_truncate_list <- function(x, n = 5L) {
  x <- as.character(x)
  if (length(x) <= n) return(paste(x, collapse = ", "))
  paste0(paste(x[seq_len(n)], collapse = ", "), ", ... (", length(x) - n, " more)")
}

# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

#' Guess which predictor produced a file
#'
#' v1 required the user to pick the format and produced an unhelpful error when
#' they picked wrong. The UI now pre-selects the detected format and only warns
#' if the user overrides it.
#'
#' @return One of "NetMHC", "NetMHCpan", "NetMHCIIpan", "MHCFlurry",
#'   "IEDB Consensus", "Other", or NA_character_ if nothing matched.
ee_detect_format <- function(path) {
  head_lines <- tryCatch(readLines(path, n = 2L, warn = FALSE), error = function(e) character(0))
  if (!length(head_lines)) return(NA_character_)
  l1 <- head_lines[1]
  l2 <- if (length(head_lines) > 1L) head_lines[2] else ""

  if (stri_detect_fixed(l1, "sequence_name") && stri_detect_fixed(l1, "best_allele")) {
    return("MHCFlurry")
  }
  if (stri_detect_fixed(l1, "seq_num") && stri_detect_fixed(l1, "allele")) {
    return("IEDB Consensus")
  }

  f1 <- stri_split_fixed(l1, "\t")[[1]]
  f2 <- stri_split_fixed(l2, "\t")[[1]]
  # NetMHC-family: row 1 holds allele labels in otherwise-empty cells, row 2 the
  # per-column names starting with Pos / Peptide / ID.
  if (length(f2) > 3L && stri_trim_both(f2[1]) %in% c("Pos", "pos")) {
    if (any(stri_detect_regex(f2, "^(EL[-_]?score|EL_Rank|BA[-_]?score|BA_Rank)$"))) {
      # icore is emitted by NetMHCpan (class I) but not NetMHCIIpan.
      if (any(stri_detect_fixed(f2, "icore"))) return("NetMHCpan")
      return("NetMHCIIpan")
    }
    if (any(f2 == "Target") || any(stri_detect_regex(f1, "D[RQP][AB]"))) return("NetMHCIIpan")
    if (any(f2 == "nM")) return("NetMHC")
    return("NetMHCpan")
  }
  if (stri_detect_fixed(l1, "\t")) return("Other")
  NA_character_
}

# ---------------------------------------------------------------------------
# NetMHC / NetMHCpan / NetMHCIIpan
# ---------------------------------------------------------------------------

# Column names, within one allele block, that hold a percentile rank or a score.
# Ordered by preference: elution-ligand rank first, then binding-affinity rank.
.EE_RANK_PATTERNS  <- c("^EL_Rank$", "^%?Rank_EL$", "^%?Rank$", "^Rank$", "^BA_Rank$", "^%?Rank_BA$")
.EE_SCORE_PATTERNS <- c("^nM$", "^Aff\\(nM\\)$", "^Affinity$", "^BA[-_]?score$", "^Score$", "^EL[-_]?score$")

#' Parse any NetMHC-family "-xls" table
#'
#' Both header rows are used: row 1 gives the allele label at the first column of
#' each allele's block, row 2 gives the column names. Block boundaries therefore
#' come from the file itself and the parser is independent of how many columns a
#' given NetMHC version emits per allele (3, 4 or 6, with or without -BA).
ee_parse_netmhc_family <- function(prediction_file, fasta, score_type,
                                   flavour = c("NetMHCpan", "NetMHC", "NetMHCIIpan")) {
  flavour <- match.arg(flavour)

  hdr <- readLines(prediction_file, n = 2L, warn = FALSE)
  if (length(hdr) < 2L) {
    stop("A ", flavour, " table needs two header rows; this file has fewer than two lines.",
         call. = FALSE)
  }
  h1 <- stri_trim_both(stri_split_fixed(hdr[1], "\t")[[1]])
  h2 <- stri_trim_both(stri_split_fixed(hdr[2], "\t")[[1]])

  # Column layout is worked out from the header alone, so that fread can be
  # told to read ONLY the columns that are actually used. A NetMHCpan table has
  # four columns per allele but only one of them is a score; the other three are
  # repeated peptide strings. Reading all of them roughly quadrupled peak memory
  # for no benefit -- enough to matter on a 1 GB container.
  probe <- readLines(prediction_file, n = 3L, warn = FALSE)
  if (length(probe) < 3L) stop("The ", flavour, " table contains no data rows.", call. = FALSE)
  ncol_body <- length(stri_split_fixed(probe[3], "\t")[[1]])

  length(h1) <- ncol_body; length(h2) <- ncol_body   # pad trimmed trailing tabs
  h1[is.na(h1)] <- ""; h2[is.na(h2)] <- ""

  allele_at <- which(nzchar(h1))
  if (!length(allele_at)) {
    stop("No MHC allele labels were found in the first header row of this ", flavour,
         " file. Export the table with the allele header intact.", call. = FALSE)
  }
  alleles <- h1[allele_at]

  # Each allele block runs to the start of the next allele; the final block ends
  # before the trailing summary columns (Ave / NB / H_Avg_Ranks / N_binders).
  summary_at <- which(stri_detect_regex(h2, "^(Ave|NB|H_Avg_Ranks|N_binders|Sum)$"))
  tail_stop  <- if (length(summary_at)) min(summary_at[summary_at > max(allele_at)], ncol_body + 1L)
                else ncol_body + 1L
  block_end  <- c(allele_at[-1] - 1L, min(tail_stop, ncol_body + 1L) - 1L)

  patterns <- if (identical(score_type, "Rank")) .EE_RANK_PATTERNS else .EE_SCORE_PATTERNS
  pick_col <- function(from, to) {
    idx <- from:to
    nm  <- h2[idx]
    for (p in patterns) {
      m <- which(stri_detect_regex(nm, p, case_insensitive = TRUE))
      if (length(m)) return(idx[m[1]])
    }
    NA_integer_
  }
  score_cols <- mapply(pick_col, allele_at, block_end)

  if (anyNA(score_cols)) {
    stop(sprintf(
      paste0("No '%s' column was found inside the allele block for %s in this %s file. ",
             "Columns available per allele: %s. Try the other Score Type."),
      score_type, ee_truncate_list(alleles[is.na(score_cols)], 3L), flavour,
      paste(unique(h2[allele_at[1]:block_end[1]]), collapse = ", ")
    ), call. = FALSE)
  }

  # Leading metadata columns, located by name rather than by index.
  meta   <- h2[seq_len(allele_at[1] - 1L)]
  i_pos  <- match(TRUE, stri_detect_regex(meta, "^pos$",     case_insensitive = TRUE))
  i_pep  <- match(TRUE, stri_detect_regex(meta, "^peptide$", case_insensitive = TRUE))
  i_id   <- match(TRUE, stri_detect_regex(meta, "^(id|protein|seq_?id)$", case_insensitive = TRUE))
  if (anyNA(c(i_pos, i_pep, i_id))) {
    stop("Expected 'Pos', 'Peptide' and 'ID' columns in the ", flavour,
         " header; found: ", paste(meta, collapse = ", "), call. = FALSE)
  }

  # Read the four metadata columns and one score column per allele, nothing else.
  keep_idx <- c(i_pos, i_pep, i_id, score_cols)
  body <- data.table::fread(
    prediction_file, sep = "\t", header = FALSE, skip = 2L,
    select = keep_idx, fill = TRUE, showProgress = FALSE, data.table = TRUE,
    colClasses = list(character = c(i_pep, i_id), numeric = c(i_pos, score_cols)),
    strip.white = TRUE, na.strings = c("", "NA", "-")
  )
  if (!nrow(body)) stop("The ", flavour, " table contains no data rows.", call. = FALSE)
  # fread returns the selected columns in file order; re-map to that order.
  ord <- order(keep_idx)
  at <- integer(length(keep_idx)); at[ord] <- seq_along(keep_idx)

  pos_raw <- suppressWarnings(as.numeric(body[[at[1]]]))
  peptide <- as.character(body[[at[2]]])
  pred_id <- as.character(body[[at[3]]])
  score_at <- at[-(1:3)]

  keep <- !is.na(pos_raw) & !is.na(peptide) & nzchar(peptide)
  if (!all(keep)) {
    body <- body[keep]; pos_raw <- pos_raw[keep]; peptide <- peptide[keep]; pred_id <- pred_id[keep]
  }
  if (!length(peptide)) stop("No usable peptide rows in the ", flavour, " table.", call. = FALSE)

  # NetMHC/NetMHCpan emit 0-based Pos, NetMHCIIpan 1-based. Detect rather than
  # assume: if the smallest position is 0 the file is 0-based.
  one_based <- min(pos_raw, na.rm = TRUE) > 0
  pos1 <- if (one_based) pos_raw else pos_raw + 1

  scores <- ee_score_matrix(body, score_at, alleles)
  ee_assemble(peptide, pos1, pred_id, scores, fasta,
              predictor = flavour, score_type = score_type)
}

#' Coerce the selected columns into a numeric allele x peptide score matrix
ee_score_matrix <- function(body, cols, alleles) {
  m <- matrix(NA_real_, nrow = nrow(body), ncol = length(cols),
              dimnames = list(NULL, ee_clean_allele(alleles)))
  for (j in seq_along(cols)) {
    m[, j] <- suppressWarnings(as.numeric(body[[cols[j]]]))
  }
  allna <- colSums(!is.na(m)) == 0L
  if (all(allna)) {
    stop("Every selected score column parsed to NA. The chosen Score Type probably ",
         "does not exist in this file.", call. = FALSE)
  }
  if (any(allna)) {
    warning(sprintf("Dropping %d allele(s) with no numeric values: %s.",
                    sum(allna), ee_truncate_list(colnames(m)[allna])), call. = FALSE)
    m <- m[, !allna, drop = FALSE]
  }
  m
}

#' Tidy an allele label without destroying it
#'
#' v1 replaced ':' and '-' with '.' so that column names were syntactically valid
#' in data.frames, which turned "HLA-A01:01" into "HLA.A01.01" in every plot and
#' download. Storing scores in a matrix removes that constraint, so the original
#' allele name is preserved and only whitespace is trimmed.
ee_clean_allele <- function(x) {
  x <- stri_trim_both(as.character(x))
  x <- stri_replace_all_regex(x, "\\s+", " ")
  make.unique(ee_canonical_allele(x), sep = " #")
}

#' Undo the dot-substitution some tools apply to allele names
#'
#' v1 rewrote ':' and '-' as '.' so allele names were valid data.frame column
#' names, and every table it exported carries that form -- "HLA.A01.01" instead
#' of "HLA-A01:01". Those names then flow back in through the generic "Other"
#' reader (the paper's own class I dataset is one such file), where they sort
#' oddly, read badly in plots, and used to defeat MHC-class detection.
#'
#' Only unambiguous, fully-dotted patterns are rewritten; anything else is left
#' exactly as the file spelled it, because inventing punctuation in an allele
#' name is worse than displaying an unusual one.
ee_canonical_allele <- function(x) {
  out <- x
  # HLA.A01.01 / HLA.DRB1.03.01 -> HLA-A01:01 / HLA-DRB1*03:01
  i <- stri_detect_regex(out, "^HLA\\.[A-Z]+[0-9]*\\.[0-9]+$", case_insensitive = TRUE)
  out[i] <- stri_replace_first_regex(out[i], "^(HLA)\\.([A-Za-z]+[0-9]*)\\.([0-9]+)$", "$1-$2:$3")
  # DRB1.03.01 -> DRB1_0301   (NetMHCIIpan's own spelling)
  j <- stri_detect_regex(out, "^(DRB[0-9]|DQA[0-9]|DQB[0-9]|DPA[0-9]|DPB[0-9])\\.([0-9]+)\\.([0-9]+)$")
  out[j] <- stri_replace_first_regex(out[j], "^([A-Z]+[0-9])\\.([0-9]+)\\.([0-9]+)$", "$1_$2$3")
  out
}

# ---------------------------------------------------------------------------
# MHCFlurry
# ---------------------------------------------------------------------------

ee_parse_mhcflurry <- function(prediction_file, fasta, score_type) {
  # Read the header alone first so only the five columns actually used are
  # loaded; MHCflurry emits a dozen, several of them long flank sequences.
  hdr <- names(data.table::fread(prediction_file, nrows = 0L, showProgress = FALSE))
  need <- c("sequence_name", "pos", "peptide", "best_allele")
  miss <- setdiff(need, hdr)
  if (length(miss)) {
    stop("This does not look like an MHCflurry output: missing column(s) ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  val_col <- if (identical(score_type, "Rank")) {
    ee_first_present_nm(hdr, c("affinity_percentile", "presentation_percentile", "processing_percentile"))
  } else {
    ee_first_present_nm(hdr, c("affinity", "presentation_score", "processing_score"))
  }
  if (is.na(val_col)) {
    stop("No ", score_type, " column found in the MHCflurry output ",
         "(looked for affinity_percentile / affinity).", call. = FALSE)
  }

  dt <- data.table::fread(prediction_file, select = c(need, val_col),
                          showProgress = FALSE, data.table = TRUE,
                          na.strings = c("", "NA", "-"))
  dt <- dt[!is.na(peptide) & nzchar(peptide) & !is.na(pos)]
  if (!nrow(dt)) stop("The MHCflurry output contains no usable rows.", call. = FALSE)

  # MHCflurry reports one row per (peptide, sample); best_allele names the winner.
  # Keep the best value per (peptide, allele) so repeated samples collapse cleanly.
  best <- if (identical(score_type, "Rank") || val_col == "affinity") min else max
  wide <- data.table::dcast(
    dt, peptide + pos + sequence_name ~ best_allele,
    value.var = val_col, fun.aggregate = function(v) if (all(is.na(v))) NA_real_ else best(v, na.rm = TRUE)
  )
  allele_cols <- setdiff(names(wide), c("peptide", "pos", "sequence_name"))
  if (!length(allele_cols)) stop("No MHC alleles found in column 'best_allele'.", call. = FALSE)

  scores <- ee_score_matrix(wide, match(allele_cols, names(wide)), allele_cols)
  ee_assemble(wide$peptide, as.numeric(wide$pos) + 1, wide$sequence_name, scores, fasta,
              predictor = "MHCFlurry", score_type = score_type)
}

ee_first_present <- function(dt, candidates) ee_first_present_nm(names(dt), candidates)

#' First candidate name present in a header, in preference order
ee_first_present_nm <- function(hdr, candidates) {
  hit <- candidates[candidates %in% hdr]
  if (length(hit)) hit[1] else NA_character_
}

# ---------------------------------------------------------------------------
# IEDB consensus
# ---------------------------------------------------------------------------

ee_parse_iedb <- function(prediction_file, fasta, score_type) {
  hdr <- names(data.table::fread(prediction_file, sep = "\t", nrows = 0L, showProgress = FALSE))
  need <- c("allele", "seq_num", "start", "peptide")
  miss <- setdiff(need, hdr)
  if (length(miss)) {
    stop("This does not look like an IEDB output: missing column(s) ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  val_col <- if (identical(score_type, "Rank")) {
    ee_first_present_nm(hdr, c("consensus_percentile_rank", "percentile_rank", "ann_rank", "smm_rank"))
  } else {
    ee_first_present_nm(hdr, c("ann_ic50", "smm_ic50", "ic50", "score"))
  }
  if (is.na(val_col)) {
    stop("No ", score_type, " column found in the IEDB output.", call. = FALSE)
  }

  dt <- data.table::fread(prediction_file, sep = "\t", select = c(need, val_col),
                          showProgress = FALSE, data.table = TRUE,
                          na.strings = c("", "NA", "-"))
  dt <- dt[!is.na(peptide) & nzchar(peptide) & !is.na(start)]
  if (!nrow(dt)) stop("The IEDB output contains no usable rows.", call. = FALSE)
  dt[, `:=`(.val = suppressWarnings(as.numeric(get(val_col))))]

  wide <- data.table::dcast(
    dt, peptide + start + seq_num ~ allele, value.var = ".val",
    fun.aggregate = function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
  )
  allele_cols <- setdiff(names(wide), c("peptide", "start", "seq_num"))
  if (!length(allele_cols)) stop("No MHC alleles found in column 'allele'.", call. = FALSE)

  # seq_num is an explicit 1-based index into the submitted FASTA.
  sn <- as.integer(wide$seq_num)
  if (max(sn, na.rm = TRUE) > nrow(fasta)) {
    stop(sprintf(paste("The IEDB file references sequence number %d but the FASTA has only",
                       "%d record(s). The FASTA is not the one used for the prediction."),
                 max(sn, na.rm = TRUE), nrow(fasta)), call. = FALSE)
  }
  pred_id <- fasta$Name[sn]

  scores <- ee_score_matrix(wide, match(allele_cols, names(wide)), allele_cols)
  ee_assemble(wide$peptide, as.numeric(wide$start), pred_id, scores, fasta,
              predictor = "IEDB Consensus", score_type = score_type)
}

# ---------------------------------------------------------------------------
# Generic / "Other"
# ---------------------------------------------------------------------------

#' Parse a user-supplied table
#'
#' Documented layout: Peptide, Position, Protein ID, Protein length, then one
#' column per allele. v1 trusted that order blindly; here the columns are found
#' by name when possible and validated either way.
ee_parse_other <- function(prediction_file, fasta, score_type) {
  dt <- data.table::fread(prediction_file, showProgress = FALSE, data.table = TRUE,
                          na.strings = c("", "NA", "-"))
  if (ncol(dt) < 5L) {
    stop("A custom table needs at least 5 columns: peptide, position, protein ID, ",
         "protein length, and one column per MHC allele. This file has ", ncol(dt), ".",
         call. = FALSE)
  }
  nm  <- stri_trans_tolower(names(dt))
  fnd <- function(p, default) {
    i <- match(TRUE, stri_detect_regex(nm, p))
    if (is.na(i)) default else i
  }
  i_pep <- fnd("^(peptide|sequence|epitope)$", 1L)
  i_pos <- fnd("^(pos|position|start)$",       2L)
  i_id  <- fnd("^(id|protein|protein_?id)$",   3L)
  i_len <- fnd("^(length|prot_?length|protein_?length)$", 4L)

  meta_idx    <- unique(c(i_pep, i_pos, i_id, i_len))
  allele_idx  <- setdiff(seq_len(ncol(dt)), meta_idx)
  if (!length(allele_idx)) stop("No allele columns left after the four metadata columns.", call. = FALSE)

  scores <- ee_score_matrix(dt, allele_idx, names(dt)[allele_idx])
  ee_assemble(dt[[i_pep]], as.numeric(dt[[i_pos]]), dt[[i_id]], scores, fasta,
              predictor = "Other", score_type = score_type,
              fallback_length = suppressWarnings(as.numeric(dt[[i_len]])))
}

# ---------------------------------------------------------------------------
# Assembly + canonical object
# ---------------------------------------------------------------------------

#' Build the canonical dataset from parsed pieces
#'
#' @param peptide,pos1,pred_id Per-row peptide sequence, 1-based start, protein ID.
#' @param scores Numeric matrix (rows aligned to peptide).
#' @param fasta data.table from ee_read_fasta(), or NULL to use fallback_length.
ee_assemble <- function(peptide, pos1, pred_id, scores, fasta,
                        predictor, score_type, fallback_length = NULL) {
  peptide <- stri_trans_toupper(stri_trim_both(as.character(peptide)))
  pred_id <- as.character(pred_id)
  klen    <- stri_length(peptide)
  notes   <- character(0)

  ok <- !is.na(peptide) & nzchar(peptide) & !is.na(pos1) & !is.na(pred_id)
  if (!all(ok)) {
    notes <- c(notes, sprintf("Dropped %d incomplete row(s).", sum(!ok)))
    peptide <- peptide[ok]; pos1 <- pos1[ok]; pred_id <- pred_id[ok]
    klen <- klen[ok]; scores <- scores[ok, , drop = FALSE]
  }
  if (!length(peptide)) stop("No usable rows remained after parsing.", call. = FALSE)

  # Length implied by the prediction file itself: predictors emit every sliding
  # window, so the last start position plus the peptide length is the protein
  # length. Keep first-appearance order -- the positional fallback depends on it.
  pid_f   <- factor(pred_id, levels = unique(pred_id))
  uid     <- levels(pid_f)
  implied <- setNames(as.numeric(tapply(pos1 + klen - 1, pid_f, max, na.rm = TRUE))[seq_along(uid)], uid)

  if (!is.null(fasta) && nrow(fasta)) {
    m <- ee_match_proteins(uid, fasta, implied_len = implied)
    map_name <- setNames(m$name,   uid)
    map_len  <- setNames(m$length, uid)
    notes    <- c(notes, m$notes)
  } else if (!is.null(fallback_length)) {
    map_name <- setNames(uid, uid)
    fl <- fallback_length[seq_along(pred_id)]
    map_len <- setNames(as.numeric(tapply(fl, pid_f, function(v) {
      v <- v[!is.na(v)]; if (length(v)) v[1] else NA_real_
    }))[seq_along(uid)], uid)
    # Fall back to the implied length wherever the table did not supply one.
    map_len[is.na(map_len)] <- implied[is.na(map_len)]
    notes <- c(notes, "No FASTA supplied: protein lengths were taken from the input table.")
  } else {
    stop("A FASTA file is required to determine protein lengths.", call. = FALSE)
  }

  peptides <- data.table(
    Peptide    = peptide,
    Pos        = as.integer(pos1),
    End        = as.integer(pos1 + klen - 1L),
    PepLength  = as.integer(klen),
    ID         = unname(map_name[pred_id]),
    ProtLength = as.integer(unname(map_len[pred_id]))
  )

  # Peptides longer than their protein, or starting past its end, mean the FASTA
  # and the prediction disagree. Report instead of producing negative densities.
  bad <- peptides[!is.na(ProtLength) & End > ProtLength, .N]
  if (bad > 0L) {
    notes <- c(notes, sprintf(
      "%d peptide(s) end past the length of their protein in the FASTA; check that the FASTA matches the prediction.",
      bad))
  }

  ee_dataset(peptides, scores, predictor, score_type, notes)
}

#' Construct the canonical epitope dataset
#'
#' Precomputes, once, the two things every tool needs and v1 recomputed on every
#' render: the unique-peptide index and the per-peptide protein count.
ee_dataset <- function(peptides, scores, predictor, score_type, notes = character(0)) {
  stopifnot(nrow(peptides) == nrow(scores))

  alleles <- colnames(scores)
  cls     <- ee_mhc_class(alleles)

  # Predictions are sequence-based, so a peptide occurring in several proteins
  # carries identical scores. Index the first occurrence of each once here.
  #
  # The factor levels MUST be in first-appearance order, not the alphabetical
  # order factor() defaults to: uniq_idx (row of each first occurrence) and
  # uniq_pep (the distinct peptides) are used as parallel vectors throughout, so
  # alphabetical levels would label every unique-peptide result with the wrong
  # sequence. tests/test_core.R asserts the alignment.
  pep_f    <- factor(peptides$Peptide, levels = unique(peptides$Peptide))
  first_at <- !duplicated(pep_f)

  proteins <- unique(peptides[, .(ID, ProtLength)])[order(ID)]

  structure(list(
    peptides    = peptides,
    scores      = scores,
    alleles     = alleles,
    mhc_class   = cls,
    predictor   = predictor,
    score_type  = score_type,
    proteins    = proteins,
    uniq_idx    = which(first_at),          # rows holding each distinct peptide
    uniq_pep    = levels(pep_f),
    pep_index   = as.integer(pep_f),        # row -> index into uniq_pep
    n_rows      = nrow(peptides),
    n_peptides  = nlevels(pep_f),
    n_proteins  = nrow(proteins),
    notes       = notes
  ), class = "ee_dataset")
}

#' Class I vs class II from allele nomenclature
ee_mhc_class <- function(alleles) {
  a <- stri_trans_toupper(alleles)
  # The separator between locus and number varies by predictor and by whatever
  # sanitising a file has been through: "HLA-A01:01", "HLA_A0101", "HLA A01:01"
  # and "HLA.A01.01" all occur. v1 rewrote ':' and '-' as '.', so any table
  # exported from it arrives in the dotted form -- which is exactly the case
  # that used to fall through to "unknown".
  sep <- "[-_. ]?"
  ii <- stri_detect_regex(a, paste0("D[RQP][AB]?[0-9_.*]|^H2", sep, "I|^H-2", sep, "I|DRB|DQA|DQB|DPA|DPB"))
  i  <- stri_detect_regex(a, paste0("^HLA", sep, "[ABCEFG]|^H2", sep, "[KDL]|^H-2", sep,
                                    "[KDL]|^(BOLA|SLA|MAMU|PATR)"))
  if (sum(ii) > sum(i)) "II" else if (sum(i) > 0L) "I" else "unknown"
}

print.ee_dataset <- function(x, ...) {
  cat(sprintf("<ee_dataset> %s / %s\n", x$predictor, x$score_type))
  cat(sprintf("  %s rows, %s unique peptides, %d proteins, %d alleles (MHC class %s)\n",
              format(x$n_rows, big.mark = ","), format(x$n_peptides, big.mark = ","),
              x$n_proteins, length(x$alleles), x$mhc_class))
  if (length(x$notes)) cat(paste0("  note: ", x$notes, collapse = "\n"), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

#' Parse a prediction file into an ee_dataset
#'
#' Warnings raised anywhere in the pipeline are collected onto the object rather
#' than printed to the console, so the UI can show them to the user.
#'
#' @param prediction_file Path to the predictor output.
#' @param fasta_file Path to the FASTA used for the prediction (or NULL for "Other").
#' @param predictor One of the supported names, or "auto" to detect.
#' @param score_type "Rank" or "Score".
ee_parse <- function(prediction_file, fasta_file, predictor = "auto", score_type = "Rank") {
  if (is.null(prediction_file) || !nzchar(prediction_file)) {
    stop("No prediction file supplied.", call. = FALSE)
  }
  if (!file.exists(prediction_file)) stop("Prediction file not found.", call. = FALSE)
  if (!identical(score_type, "Rank") && !identical(score_type, "Score")) {
    stop("score_type must be \"Rank\" or \"Score\".", call. = FALSE)
  }

  if (identical(predictor, "auto")) {
    predictor <- ee_detect_format(prediction_file)
    if (is.na(predictor)) {
      stop("Could not recognise this file format. Choose the predictor explicitly, ",
           "or use \"Other\" with the documented column layout.", call. = FALSE)
    }
  }

  fasta <- if (!is.null(fasta_file) && nzchar(fasta_file)) ee_read_fasta(fasta_file) else NULL
  if (is.null(fasta) && !identical(predictor, "Other")) {
    stop("A FASTA file is required for ", predictor, " output.", call. = FALSE)
  }

  warns <- character(0)
  ds <- withCallingHandlers(
    switch(predictor,
      "NetMHC"         = ee_parse_netmhc_family(prediction_file, fasta, score_type, "NetMHC"),
      "NetMHCpan"      = ee_parse_netmhc_family(prediction_file, fasta, score_type, "NetMHCpan"),
      "NetMHCIIpan"    = ee_parse_netmhc_family(prediction_file, fasta, score_type, "NetMHCIIpan"),
      "MHCFlurry"      = ee_parse_mhcflurry(prediction_file, fasta, score_type),
      "IEDB Consensus" = ee_parse_iedb(prediction_file, fasta, score_type),
      "Other"          = ee_parse_other(prediction_file, fasta, score_type),
      stop("Unknown predictor: ", predictor, call. = FALSE)
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  ds$notes <- c(ds$notes, warns)
  ds
}
