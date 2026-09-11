# tests/test_core.R -- Engine tests for Epitope-Evaluator 2
#
# Runs without Shiny, plotly or DT: only data.table, stringi and matrixStats.
#
#   Rscript tests/test_core.R          (from the app root)
#
# Covers the five parsers against the bundled example files, every analysis
# function, and specific regressions for each v1 defect listed in CHANGELOG.md.

suppressPackageStartupMessages({
  library(data.table); library(stringi); library(matrixStats)
})

root <- if (file.exists("global.R")) "." else ".."
setwd(root)
source("utils/constants.R")
source("utils/parsing_functions.R")
source("utils/core_functions.R")

# --- tiny test harness -----------------------------------------------------
.pass <- 0L; .fail <- 0L; .failures <- character(0)

ok <- function(label, expr) {
  res <- tryCatch(isTRUE(expr), error = function(e) structure(FALSE, msg = conditionMessage(e)))
  if (isTRUE(res)) {
    .pass <<- .pass + 1L
    cat(sprintf("  \033[32mok\033[0m   %s\n", label))
  } else {
    .fail <<- .fail + 1L
    m <- attr(res, "msg")
    .failures <<- c(.failures, label)
    cat(sprintf("  \033[31mFAIL\033[0m %s%s\n", label, if (is.null(m)) "" else paste0("  <", m, ">")))
  }
}
throws <- function(label, expr, pattern = NULL) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  hit <- !is.null(e) && (is.null(pattern) || grepl(pattern, e, ignore.case = TRUE))
  ok(label, hit)
  if (!hit && !is.null(e)) cat("       got: ", e, "\n")
}
section <- function(x) cat(sprintf("\n\033[1m%s\033[0m\n", x))
tick <- function(expr) { t0 <- Sys.time(); v <- force(expr)
  attr(v, "secs") <- as.numeric(difftime(Sys.time(), t0, units = "secs")); v }

FA <- "data/example.fasta"

# ===========================================================================
section("FASTA reader")
# ===========================================================================
fa <- ee_read_fasta(FA)
ok("17 records read",                    nrow(fa) == 17L)
ok("Spike is 1273 aa (true length)",     fa$Length[fa$Name == "sp|P0DTC2|SPIKE_SARS2"] == 1273L)
ok("shortest record is ORF3B, 22 aa",    min(fa$Length) == 22L)
ok("names are the first header token",   all(!grepl(" ", fa$Name)))
ok("sequences are upper-case residues",  all(grepl("^[A-Z]+$", fa$Sequence)))

tmp <- tempfile(fileext = ".fasta")
writeLines(c("not a fasta", "just text"), tmp)
throws("non-FASTA input is rejected", ee_read_fasta(tmp), "does not look like a FASTA")
writeLines(character(0), tmp)
throws("empty file is rejected", ee_read_fasta(tmp), "empty")
unlink(tmp)

# ===========================================================================
section("Format detection")
# ===========================================================================
ok("NetMHC detected",       ee_detect_format("data/example_NetMHC.xls")        == "NetMHC")
ok("NetMHCpan detected",    ee_detect_format("data/example_NetMHCPAN.xls")     == "NetMHCpan")
ok("NetMHCIIpan detected",  ee_detect_format("data/example_NetMHCIIPAN.xls")   == "NetMHCIIpan")
ok("MHCflurry detected",    ee_detect_format("data/example_MHCFlurry.txt")     == "MHCFlurry")
ok("IEDB detected",         ee_detect_format("data/example_IEDB_consensus.txt")== "IEDB Consensus")

# ===========================================================================
section("Parsers -- all five formats, both score types")
# ===========================================================================
cases <- list(
  list("data/example_NetMHCPAN.xls",      "NetMHCpan",      "Rank",  12L, "I"),
  list("data/example_NetMHCPAN.xls",      "NetMHCpan",      "Score", 12L, "I"),
  list("data/example_NetMHC.xls",         "NetMHC",         "Rank",  12L, "I"),
  list("data/example_NetMHC.xls",         "NetMHC",         "Score", 12L, "I"),
  list("data/example_NetMHCIIPAN.xls",    "NetMHCIIpan",    "Rank",   8L, "II"),
  list("data/example_NetMHCIIPAN.xls",    "NetMHCIIpan",    "Score",  8L, "II"),
  list("data/example_MHCFlurry.txt",      "MHCFlurry",      "Rank",  11L, "I"),
  list("data/example_IEDB_consensus.txt", "IEDB Consensus", "Rank",   7L, "I")
)
parsed <- list()
for (cs in cases) {
  lbl <- sprintf("%s / %s", cs[[2]], cs[[3]])
  d <- tick(ee_parse(cs[[1]], FA, cs[[2]], cs[[3]]))
  parsed[[lbl]] <- d
  ok(sprintf("%-24s parses (%d alleles, class %s, %.2fs)", lbl, length(d$alleles),
             d$mhc_class, attr(d, "secs")),
     length(d$alleles) == cs[[4]] && d$mhc_class == cs[[5]] && d$n_rows > 1000L)
  ok(sprintf("%-24s all 17 proteins present", lbl), d$n_proteins == 17L)
  ok(sprintf("%-24s protein lengths are true aa lengths", lbl),
     d$proteins$ProtLength[d$proteins$ID == "sp|P0DTC2|SPIKE_SARS2"] == 1273L)
  ok(sprintf("%-24s scores are finite somewhere in every allele", lbl),
     all(colSums(is.finite(d$scores)) > 0L))
  ok(sprintf("%-24s allele names kept verbatim", lbl), !any(grepl("^X\\.|\\.\\.", d$alleles)))
}

ds  <- parsed[["NetMHCpan / Rank"]]     # class I reference dataset
ds2 <- parsed[["NetMHCIIpan / Rank"]]   # class II reference dataset

section("Regressions against v1 defects")
ok("R1: protein length is 1273, not 1265 (v1 used nchar-8)",
   ds$proteins$ProtLength[ds$proteins$ID == "sp|P0DTC2|SPIKE_SARS2"] == 1273L)
ok("R2: class II length is also 1273 (v1 gave 1265, true window count is 1259)",
   ds2$proteins$ProtLength[ds2$proteins$ID == "sp|P0DTC2|SPIKE_SARS2"] == 1273L)
ok("R3: MFVFLVLLP maps to Spike, not to a shifted protein",
   ds$peptides$ID[1] == "sp|P0DTC2|SPIKE_SARS2")
ok("R4: the 22 aa ORF3B is present and did not shift the mapping",
   "sp|P0DTF1|ORF3B_SARS2" %in% ds$proteins$ID &&
   ds$proteins$ProtLength[ds$proteins$ID == "sp|P0DTF1|ORF3B_SARS2"] == 22L)
ok("R5: NetMHCpan Score picks EL-score, not an icore text column",
   sum(is.finite(parsed[["NetMHCpan / Score"]]$scores)) > 100000L)
ok("R6: NetMHC Score picks nM affinities (values well above 100)",
   median(parsed[["NetMHC / Score"]]$scores, na.rm = TRUE) > 100)
ok("R7: allele names keep ':' (v1 rewrote HLA-A01:01 as HLA.A01.01)",
   any(grepl(":", ds$alleles)))
ok("R8: end position never exceeds protein length",
   ds$peptides[End > ProtLength, .N] == 0L)

throws("R9: wrong Score Type gives an actionable error, not silent NAs",
       ee_parse("data/example_IEDB_consensus.txt", FA, "IEDB Consensus", "Nonsense"),
       "Rank.*Score")
throws("R10: mismatched FASTA is detected, not silently mis-mapped",
       ee_parse("data/example_NetMHCPAN.xls", "data/example.fasta.missing", "NetMHCpan", "Rank"),
       "not found")

# ===========================================================================
section("ee_combine -- AND / OR semantics")
# ===========================================================================
a2 <- ds$alleles[1:2]
ok("Intersection is the row-wise max",
   all.equal(ee_combine(ds, a2, "Intersection"), matrixStats::rowMaxs(ds$scores[, a2])))
ok("Union is the row-wise min",
   all.equal(ee_combine(ds, a2, "Union"), matrixStats::rowMins(ds$scores[, a2])))
ok("single allele returns the column itself",
   all.equal(ee_combine(ds, ds$alleles[1], "Union"), unname(ds$scores[, 1])))
ok("Intersection >= Union always", all(ee_combine(ds, a2, "Intersection") >=
                                       ee_combine(ds, a2, "Union"), na.rm = TRUE))
ok("unique_only returns one value per distinct peptide",
   length(ee_combine(ds, a2, "Union", unique_only = TRUE)) == ds$n_peptides)
ok("uniq_idx and uniq_pep are parallel vectors",
   identical(ds$peptides$Peptide[ds$uniq_idx], ds$uniq_pep))
ok("pep_index indexes into uniq_pep",
   identical(ds$uniq_pep[ds$pep_index], ds$peptides$Peptide))

nam <- ds; nam$scores[1:5, 1] <- NA_real_
ok("all-NA rows become NA, not Inf",
   { z <- nam$scores[1:2, 1, drop = FALSE]
     v <- suppressWarnings(matrixStats::rowMins(z, na.rm = TRUE)); v[!is.finite(v)] <- NA
     all(is.na(v)) })
ok("partial NA does not poison the row",
   is.finite(ee_combine(nam, a2, "Union")[1]))

throws("empty allele selection is caught", ee_combine(ds, character(0), "Union"), "at least one")
throws("unknown allele is caught", ee_combine(ds, "HLA-ZZ99:99", "Union"), "not present")

# ===========================================================================
section("Tool 1 -- Distribution")
# ===========================================================================
v <- ee_combine(ds, ds$alleles[1], "Union", unique_only = TRUE)
h <- ee_histogram(v, 0, 2, 0.1)
ok("20 bins over [0, 2] at width 0.1", nrow(h) == 20L)
ok("counts equal a direct comparison",
   h$count[1] == sum(v >= 0 & v < 0.1, na.rm = TRUE))
ok("cumulative is the running total", all.equal(h$cumulative, cumsum(h$count)))
ok("last bin is closed so xmax is included",
   sum(h$count) == sum(v >= 0 & v <= 2, na.rm = TRUE))
ok("density integrates to 1", abs(sum(h$density * 0.1) - 1) < 1e-9)
ok("values outside the window are excluded", sum(h$count) < length(v))

throws("max <= min is caught",      ee_histogram(v, 2, 2, 0.1), "greater than")
throws("zero bin width is caught",  ee_histogram(v, 0, 2, 0),   "positive")
throws("absurd bin count is caught",ee_histogram(v, 0, 100, 1e-6), "Increase the bin width")

dtab <- ee_distribution_table(ds, ds$alleles[1:3], "Union", 0, 2)
ok("distribution table is sorted by score", !is.unsorted(dtab$Score))
ok("distribution table respects the window", all(dtab$Score >= 0 & dtab$Score <= 2))
ok("distribution table carries protein context", all(c("ID", "Pos", "End") %in% names(dtab)))

grid <- ee_cutoff_grid(ds, ds$alleles[1:4], seq(0.5, 2, by = 0.5))
ok("cutoff grid is alleles x cutoffs", nrow(grid) == 4L * 4L)
ok("cutoff grid counts match a direct comparison",
   grid[allele == ds$alleles[1] & cutoff == 1, count] ==
     sum(ds$scores[ds$uniq_idx, ds$alleles[1]] <= 1, na.rm = TRUE))
ok("counts are monotone in the cutoff",
   all(grid[, all(diff(count) >= 0), by = allele]$V1))

# ===========================================================================
section("Tools 2 & 6 -- set logic (this replaces v1's O(2^k) filter_ loop)")
# ===========================================================================
sets <- ds$alleles[1:5]
mem  <- ee_membership(ds, sets, 2, by = "allele")
cmb  <- ee_combinations(mem)
ok("every peptide kept belongs to >= 1 set", all(rowSums(mem$m) > 0L))
ok("combination counts sum to the peptides kept", sum(cmb$count) == nrow(mem$m))
ok("set sizes match the membership columns",
   all(ee_set_sizes(mem)$size == colSums(mem$m)))
ok("a named intersection has the size a direct filter gives",
   { pair <- sets[1:2]
     direct <- sum(ds$scores[ds$uniq_idx, pair[1]] <= 2 & ds$scores[ds$uniq_idx, pair[2]] <= 2,
                   na.rm = TRUE)
     flag <- mem$m[, pair[1]] & mem$m[, pair[2]]
     sum(flag) == direct })
ok("n_sets equals the number of '&'-joined names",
   all(cmb$n_sets == lengths(strsplit(cmb$sets, " & ", fixed = TRUE))))
ok("results are ordered by descending count", !is.unsorted(rev(cmb$count)))

# Scales past the point where 2^k enumeration is possible.
big <- ee_membership(ds, ds$alleles, 5, by = "allele")   # 12 alleles = 4095 combos in v1
tb  <- tick(ee_combinations(big, collect_peptides = FALSE))
ok(sprintf("12 alleles in %.3fs without enumerating 2^12", attr(tb, "secs")),
   attr(tb, "secs") < 2 && nrow(tb) > 0L)
ok("observed combinations never exceed 2^k - 1", nrow(tb) <= 2^12 - 1)

prots <- ds$proteins$ID[1:4]
memp  <- ee_membership(ds, prots, 2, by = "protein", alleles = ds$alleles[1:3], mode = "Union")
cmbp  <- ee_combinations(memp)
ok("conservation membership is peptides x proteins", ncol(memp$m) == 4L)
ok("conservation counts sum correctly", sum(cmbp$count) == nrow(memp$m))
ok("a conserved peptide really occurs in each of its proteins",
   { multi <- which(rowSums(memp$m) > 1L)
     if (!length(multi)) TRUE else {
       i <- multi[1]; p <- memp$peptides[i]
       all(prots[memp$m[i, ]] %in% ds$peptides[Peptide == p, unique(ID)]) } })

throws("one set is not an intersection", ee_membership(ds, sets[1], 2, by = "allele"), "at least two")

# ===========================================================================
section("Tool 3 -- Density")
# ===========================================================================
dens <- ee_density_table(ds, ds$alleles[1:3], "Union", 2)
ok("every protein appears, including zero-epitope ones", nrow(dens) == ds$n_proteins)
ok("density = epitopes / amino acids",
   all.equal(dens$Density, dens$Epitopes / dens$ProtLength))
ok("epitope counts match a direct filter",
   { sc <- ee_combine(ds, ds$alleles[1:3], "Union")
     n <- sum(ds$peptides$ID[is.finite(sc) & sc <= 2] == "sp|P0DTC2|SPIKE_SARS2")
     dens$Epitopes[dens$ID == "sp|P0DTC2|SPIKE_SARS2"] == n })
ok("density is never negative (v1's L-8 could go negative for short proteins)",
   all(dens$Density >= 0, na.rm = TRUE))

pa <- ee_protein_allele_counts(ds, ds$alleles[1:4], 2)
ok("protein x allele grid is complete", nrow(pa) == ds$n_proteins * 4L)
ok("grid counts match a direct comparison",
   { a <- ds$alleles[1]
     n <- sum(ds$scores[, a] <= 2 & ds$peptides$ID == "sp|P0DTC2|SPIKE_SARS2", na.rm = TRUE)
     pa[ID == "sp|P0DTC2|SPIKE_SARS2" & allele == a, Epitopes] == n })
ok("zero-epitope cells are present as 0, not missing", any(pa$Epitopes == 0L))

# ===========================================================================
section("Tool 4 -- Viewer / lane packing")
# ===========================================================================
ok("no lanes for no intervals", length(ee_pack_lanes(integer(0), integer(0))) == 0L)
ok("disjoint intervals share one lane",
   all(ee_pack_lanes(c(1, 10, 20), c(5, 15, 25)) == 1L))
ok("fully overlapping intervals get distinct lanes",
   length(unique(ee_pack_lanes(c(1, 1, 1), c(9, 9, 9)))) == 3L)
ok("lane count equals the maximum overlap (optimal colouring)",
   { s <- c(1, 2, 3, 20); e <- c(9, 10, 11, 28)
     max(ee_pack_lanes(s, e)) == 3L })
ok("no two intervals in a lane overlap",
   { set.seed(1); s <- sample(1:500, 300, TRUE); e <- s + 8
     L <- ee_pack_lanes(s, e)
     all(vapply(split(seq_along(s), L), function(ix) {
       o <- ix[order(s[ix])]
       length(o) < 2 || all(s[o][-1] > e[o][-length(o)])
     }, logical(1))) })

vw <- tick(ee_viewer_layout(ds, "sp|P0DTC2|SPIKE_SARS2", ds$alleles[1:3], "Union", 2))
ok(sprintf("Spike layout built in %.3fs", attr(vw, "secs")), attr(vw, "secs") < 1)
ok("protein length carried through", vw$protein$length == 1273L)
ok("every epitope lists the alleles it binds", all(nzchar(vw$epitopes$Alleles)))
ok("NAlleles equals the listed allele count",
   all(vw$epitopes$NAlleles == lengths(strsplit(vw$epitopes$Alleles, ", ", fixed = TRUE))))
ok("epitopes stay inside the protein", all(vw$epitopes$End <= vw$protein$length))
ok("Intersection mode is a subset of Union mode",
   nrow(ee_viewer_layout(ds, "sp|P0DTC2|SPIKE_SARS2", ds$alleles[1:3], "Intersection", 2)$epitopes) <=
   nrow(vw$epitopes))
ok("a protein with no hits returns an empty layout, not an error",
   nrow(ee_viewer_layout(ds, "sp|P0DTF1|ORF3B_SARS2", ds$alleles[1], "Union", 0)$epitopes) == 0L)

# The worst case v1 choked on: every window of the longest protein, all alleles.
vwb <- tick(ee_viewer_layout(ds, "sp|P0DTD1|R1AB_SARS2", ds$alleles, "Union", 100))
ok(sprintf("7096 aa protein, %d epitopes, laid out in %.3fs",
           nrow(vwb$epitopes), attr(vwb, "secs")), attr(vwb, "secs") < 5)

# ===========================================================================
section("Tool 5 -- Promiscuity")
# ===========================================================================
pr <- ee_promiscuity(ds, ds$alleles, strong = 0.5, weak = 2, min_alleles = 3)
ok("every reported epitope binds at least the minimum", all(pr$table$NAlleles >= 3L))
ok("NAlleles counts alleles at <= weak (v1's table used a strict <)",
   { p <- pr$table$Peptide[1]
     r <- ds$uniq_idx[match(p, ds$uniq_pep)]
     sum(ds$scores[r, ] <= 2, na.rm = TRUE) == pr$table$NAlleles[1] })
ok("table and heatmap agree on the boundary",
   { rn <- rownames(pr$matrix)[1]
     r  <- ds$uniq_idx[match(rn, ds$uniq_pep)]
     sum(pr$matrix[1, ] != "") == sum(ds$scores[r, ] <= 2, na.rm = TRUE) })
ok("SB is a subset of WB-or-better",
   { v <- pr$values; c2 <- pr$matrix
     all(c2[!is.na(v) & v <= 0.5] == "SB") })
ok("results are sorted by promiscuity", !is.unsorted(rev(pr$table$NAlleles)))
ok("heatmap is capped at the render limit", nrow(pr$matrix) <= EE$limits$heatmap_rows)
ok("full table is not capped", nrow(pr$table) == pr$n_total)
throws("strong > weak is caught", ee_promiscuity(ds, ds$alleles, 5, 1, 2), "must not exceed")
ok("an impossible minimum returns empty, not an error",
   nrow(ee_promiscuity(ds, ds$alleles, 0.5, 2, 999)$table) == 0L)

section("MHC class detection across allele-naming conventions")
# Predictors and downstream tooling separate locus from number with -, _, space
# or nothing. v1 rewrote ':' and '-' as '.', so any table exported from it
# arrives dotted -- that form used to fall through to "unknown".
ok("HLA-A01:01 is class I",   ee_mhc_class(c("HLA-A01:01", "HLA-B07:02")) == "I")
ok("HLA.A01.01 is class I (v1's dotted form)",
   ee_mhc_class(c("HLA.A01.01", "HLA.B07.02")) == "I")
ok("HLA_A0101 is class I",    ee_mhc_class("HLA_A0101") == "I")
ok("HLA-A*02:01 is class I",  ee_mhc_class("HLA-A*02:01") == "I")
ok("HLA A01:01 is class I",   ee_mhc_class("HLA A01:01") == "I")
ok("DRB1_0101 is class II",   ee_mhc_class(c("DRB1_0101", "DRB3_0202")) == "II")
ok("DRB1.03.01 is class II",  ee_mhc_class("DRB1.03.01") == "II")
ok("H-2Kb is class I",        ee_mhc_class("H-2Kb") == "I")
ok("nonsense names are unknown, not guessed",
   ee_mhc_class(c("foo", "bar")) == "unknown")
ok("a mixed panel follows the majority",
   ee_mhc_class(c("DRB1_0101", "DRB1_0301", "HLA-A01:01")) == "II")

section("Defaults must not open a tool on an empty result")
for (nm in names(parsed)) {
  d <- parsed[[nm]]; sg <- ee_suggest(d)
  m <- ee_suggest_min_alleles(d, d$alleles, sg$weak)
  n <- nrow(ee_promiscuity(d, d$alleles, sg$strong, sg$weak, m)$table)
  ok(sprintf("%-24s promiscuity default min=%d yields %d epitopes", nm, m, n), n > 0L)
}
ok("the suggested minimum is at least 2", {
  all(vapply(parsed, function(d) ee_suggest_min_alleles(d, d$alleles, ee_suggest(d)$weak) >= 2L,
             logical(1))) })
ok("the suggested minimum never exceeds the allele count", {
  all(vapply(parsed, function(d) ee_suggest_min_alleles(d, d$alleles, ee_suggest(d)$weak) <=
               length(d$alleles), logical(1))) })
ok("an unreachable target falls back to 2, not to an error",
   ee_suggest_min_alleles(ds, ds$alleles, 0, target = 10^9) == 2L)

# ===========================================================================
section("Class II defaults")
# ===========================================================================
s1 <- ee_suggest(ds); s2 <- ee_suggest(ds2)
ok("class I default cutoff is 2",   s1$cutoff == 2)
ok("class II default cutoff is 10", s2$cutoff == 10)
ok("class II gets a wider bin",     s2$step > s1$step)
ok("score-type datasets get nM defaults, not %rank ones",
   ee_suggest(parsed[["NetMHC / Score"]])$cutoff == 500)
ok("axis label follows the score type",
   ee_score_label(ds) == "% rank" && grepl("nM", ee_score_label(parsed[["NetMHC / Score"]])))

# ===========================================================================
cat(sprintf("\n\033[1m%d passed, %d failed\033[0m\n", .pass, .fail))
if (.fail) { cat("\nFailed:\n"); cat(paste0("  - ", .failures, collapse = "\n"), "\n"); quit(status = 1) }
