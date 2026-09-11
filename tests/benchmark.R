# tests/benchmark.R -- scale check for Epitope-Evaluator 2
#
#   Rscript tests/benchmark.R [n_peptides] [n_alleles]
#
# Generates a synthetic proteome + prediction file at the target scale
# (~10x the bundled SARS-CoV-2 example), then times parsing and every tool.
# Where v1's algorithm can still be run on current R, it is timed alongside.

suppressPackageStartupMessages({
  library(data.table); library(stringi); library(matrixStats)
})
if (!file.exists("global.R")) setwd("..")
source("utils/constants.R"); source("utils/parsing_functions.R"); source("utils/core_functions.R")

args <- commandArgs(trailingOnly = TRUE)
N_PEP  <- if (length(args) >= 1) as.integer(args[1]) else 1000000L
N_ALL  <- if (length(args) >= 2) as.integer(args[2]) else 20L
N_PROT <- 400L
K      <- 9L

tmp <- file.path(tempdir(), "ee_bench"); dir.create(tmp, showWarnings = FALSE)
fa_path   <- file.path(tmp, "bench.fasta")
pred_path <- file.path(tmp, "bench_netmhcpan.xls")

timing <- list()
mb <- function(g, col) sum(g[, col] * c(56, 8)) / 1024^2
time_it <- function(label, expr) {
  before <- mb(gc(reset = TRUE, full = TRUE), "used")
  t <- system.time(v <- force(expr))[["elapsed"]]
  peak <- mb(gc(), "max used")
  timing[[length(timing) + 1L]] <<- data.table(step = label, seconds = t,
                                               peak_mb = peak, extra_mb = peak - before)
  cat(sprintf("  %-52s %7.2f s   +%5.0f MB (peak %.0f MB)\n", label, t, peak - before, peak))
  invisible(v)
}

# ---------------------------------------------------------------------------
cat(sprintf("\nGenerating a synthetic dataset: %s peptides x %d alleles, %d proteins\n",
            format(N_PEP, big.mark = ","), N_ALL, N_PROT))
# ---------------------------------------------------------------------------
set.seed(42)
per_prot <- rep(ceiling(N_PEP / N_PROT), N_PROT)
per_prot[N_PROT] <- N_PEP - sum(per_prot[-N_PROT])
per_prot <- pmax(per_prot, 1L)

AA <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

time_it("generate FASTA", {
  con <- file(fa_path, "w")
  for (i in seq_len(N_PROT)) {
    L <- per_prot[i] + K - 1L
    s <- stri_flatten(sample(AA, L, TRUE))
    writeLines(c(sprintf(">sp|BENCH%05d|PROT%05d_TEST Synthetic protein %d", i, i, i),
                 substring(s, seq(1, L, 60), pmin(seq(1, L, 60) + 59, L))), con)
  }
  close(con)
})

alleles <- sprintf("HLA-%s%02d:%02d", rep(c("A", "B", "C"), length.out = N_ALL),
                   seq_len(N_ALL), seq_len(N_ALL))

time_it("write a NetMHCpan-format prediction file", {
  # Two header rows, then core/icore/EL-score/EL_Rank per allele, as NetMHCpan emits.
  h1 <- c("", "", "", as.vector(rbind(alleles, "", "", "")), "")
  h2 <- c("Pos", "Peptide", "ID", rep(c("core", "icore", "EL-score", "EL_Rank"), N_ALL), "Ave", "NB")

  pos <- unlist(lapply(per_prot, function(n) seq_len(n) - 1L), use.names = FALSE)
  pid <- rep(sprintf("BENCH%05d_PROT%05d", seq_len(N_PROT), seq_len(N_PROT)), per_prot)
  pep <- stri_flatten(sample(AA, N_PEP * K, TRUE))
  pep <- stri_sub(pep, seq(1, N_PEP * K, K), length = K)

  dt <- data.table(Pos = pos, Peptide = pep, ID = pid)
  for (a in alleles) {
    # A realistic rank distribution: mostly non-binders, a heavy tail of good ones.
    r <- round(100 * stats::rbeta(N_PEP, 0.6, 2.2), 4)
    dt[, (paste0(a, "_core"))  := pep]
    dt[, (paste0(a, "_icore")) := pep]
    dt[, (paste0(a, "_el"))    := round(stats::runif(N_PEP), 6)]
    dt[, (paste0(a, "_rank"))  := r]
  }
  dt[, Ave := 0][, NB := 0L]

  con <- file(pred_path, "w")
  writeLines(c(paste(h1, collapse = "\t"), paste(h2, collapse = "\t")), con)
  close(con)
  fwrite(dt, pred_path, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE)
})

sz <- file.size(pred_path) / 1024^2
cat(sprintf("  prediction file: %.0f MB   FASTA: %.1f MB\n", sz, file.size(fa_path) / 1024^2))

# Drop everything the generator allocated, then reset the GC high-water mark so
# the memory figure at the end reflects the application, not this script.
rm(list = intersect(ls(), c("dt", "pep", "pos", "pid", "con", "h1", "h2")))
invisible(gc(reset = TRUE, full = TRUE))
mem_baseline <- sum(gc()[, "used"] * c(56, 8)) / 1024^2

# ---------------------------------------------------------------------------
cat("\nParsing\n")
# ---------------------------------------------------------------------------
ds <- time_it(sprintf("ee_parse (%.0f MB NetMHCpan)", sz),
              ee_parse(pred_path, fa_path, "NetMHCpan", "Rank"))
cat(sprintf("  -> %s rows, %s unique peptides, %d proteins, %d alleles\n",
            format(ds$n_rows, big.mark = ","), format(ds$n_peptides, big.mark = ","),
            ds$n_proteins, length(ds$alleles)))
cat(sprintf("  -> resident size of the parsed object: %.0f MB\n",
            as.numeric(utils::object.size(ds)) / 1024^2))

sel  <- ds$alleles
sel8 <- ds$alleles[seq_len(min(8L, N_ALL))]
CUT  <- 2

# ---------------------------------------------------------------------------
cat("\nTool 1 -- Distribution\n")
# ---------------------------------------------------------------------------
v <- time_it("ee_combine, all alleles, Union", ee_combine(ds, sel, "Union", unique_only = TRUE))
time_it("ee_histogram, 200 bins", ee_histogram(v, 0, 20, 0.1))
time_it("ee_cutoff_grid, all alleles x 40 cutoffs",
        ee_cutoff_grid(ds, sel, seq(0.05, 2, length.out = 40)))
time_it("  [v1 equivalent] per-cutoff full column scan", {
  out <- vector("list", 40); cs <- seq(0.05, 2, length.out = 40)
  for (i in seq_along(cs)) out[[i]] <- colSums(ds$scores[ds$uniq_idx, sel, drop = FALSE] <= cs[i],
                                               na.rm = TRUE)
  out })
time_it("ee_distribution_table", ee_distribution_table(ds, sel, "Union", 0, CUT))

# ---------------------------------------------------------------------------
cat("\nTools 2 & 6 -- set logic\n")
# ---------------------------------------------------------------------------
mem <- time_it(sprintf("ee_membership, %d alleles", length(sel)),
               ee_membership(ds, sel, CUT, by = "allele"))
cmb <- time_it("ee_combinations (no 2^k enumeration)",
               ee_combinations(mem, collect_peptides = FALSE))
cat(sprintf("  -> %s peptides in %s distinct combinations (2^%d - 1 = %s possible)\n",
            format(nrow(mem$m), big.mark = ","), format(nrow(cmb), big.mark = ","),
            length(sel), format(2^length(sel) - 1, big.mark = ",", scientific = FALSE)))
time_it("ee_combinations with peptide lists", ee_combinations(mem))

prots <- ds$proteins$ID[1:8]
time_it("ee_membership over 8 proteins (Conservation)",
        ee_membership(ds, prots, CUT, by = "protein", alleles = sel8, mode = "Union"))

# ---------------------------------------------------------------------------
cat("\nTool 3 -- Density\n")
# ---------------------------------------------------------------------------
time_it("ee_density_table", ee_density_table(ds, sel, "Union", CUT))
time_it(sprintf("ee_protein_allele_counts (%d x %d grid)", N_PROT, length(sel)),
        ee_protein_allele_counts(ds, sel, CUT))
time_it("  [v1 equivalent] melt to long form then group", {
  long <- data.table::melt(
    data.table(ID = ds$peptides$ID, as.data.table(ds$scores[, sel, drop = FALSE])),
    id.vars = "ID", variable.name = "allele")
  long[value <= CUT, .N, by = .(ID, allele)] })

# ---------------------------------------------------------------------------
cat("\nTool 4 -- Viewer\n")
# ---------------------------------------------------------------------------
big_prot <- ds$proteins[which.max(ProtLength), ID]
lay <- time_it(sprintf("ee_viewer_layout on the longest protein (%s aa)",
                       format(max(ds$proteins$ProtLength), big.mark = ",")),
               ee_viewer_layout(ds, big_prot, sel, "Union", 100))
cat(sprintf("  -> %s epitopes packed into %d rows\n",
            format(nrow(lay$epitopes), big.mark = ","), lay$n_lanes))

# v1's lane assignment was a nested pairwise overlap search. Timed on a small
# slice, because at full size it does not finish in a reasonable interval.
v1_pack <- function(start, end) {
  n <- length(start); lev <- rep(0L, n); id <- seq_len(n); level <- 0L
  while (level %in% lev) {
    idx <- which(lev == level)
    changed <- TRUE
    while (changed && length(idx) > 1L) {
      moved <- FALSE
      for (i in seq_along(idx)) {
        a <- idx[i]; vs <- seq(start[a], end[a])
        for (b in idx[-i]) {
          if (length(intersect(vs, seq(start[b], end[b]))) != 0L) { lev[b] <- lev[b] - 1L; moved <- TRUE }
        }
        idx <- which(lev == level)
        if (moved) break
      }
      changed <- moved
    }
    level <- level - 1L
    if (level < -60L) break
  }
  lev
}
ep <- lay$epitopes
for (n in c(200L, 800L)) {
  if (nrow(ep) < n) next
  s <- ep$Pos[seq_len(n)]; e <- ep$End[seq_len(n)]
  time_it(sprintf("  ee_pack_lanes, %d intervals", n),        ee_pack_lanes(s, e))
  time_it(sprintf("  [v1 equivalent] pairwise, %d intervals", n), v1_pack(s, e))
}

# ---------------------------------------------------------------------------
cat("\nTool 5 -- Promiscuity\n")
# ---------------------------------------------------------------------------
time_it("ee_promiscuity", ee_promiscuity(ds, sel, 0.5, 2, max(2L, length(sel) - 4L)))

# ---------------------------------------------------------------------------
res <- rbindlist(timing)
cat("\n", strrep("-", 64), "\n", sep = "")
cat(sprintf("Total analysis time after parsing: %.2f s\n",
            sum(res$seconds[!grepl("^generate|^write|ee_parse|v1 equiv", res$step)])))
cat(sprintf("Slowest analysis step: %s (%.2f s)\n",
            res$step[which.max(res$seconds * !grepl("^generate|^write|v1 equiv", res$step))],
            max(res$seconds * !grepl("^generate|^write|v1 equiv", res$step))))
app <- res[!grepl("^generate|^write|v1 equivalent", step)]
cat(sprintf("Baseline after parsing: %.0f MB (parsed object %.0f MB)\n",
            mem_baseline, as.numeric(utils::object.size(ds)) / 1024^2))
cat(sprintf("Highest peak in any application step: %.0f MB  (%s)\n",
            max(app$peak_mb), app$step[which.max(app$peak_mb)]))
unlink(tmp, recursive = TRUE)
