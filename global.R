# global.R -- Global configuration for Epitope-Evaluator 2
#
# Dependency policy: CRAN-only, no system libraries (GDAL/GEOS/PROJ/etc.), so the
# app deploys unchanged to shinyapps.io. Seven packages total, down from 22 in v1:
#
#   shiny bslib plotly DT data.table matrixStats stringi
#
# Deliberately removed relative to v1:
#   seqinr, phylotools  -> replaced by a vectorised built-in FASTA reader
#   reshape, reshape2   -> archived on CRAN; replaced by data.table
#   ggVennDiagram, sf   -> sf needs GDAL/GEOS/PROJ; Venn is now drawn natively in plotly
#   ggplot2, gridExtra  -> every plot is now a native plotly trace (no ggplotly round-trip)
#   shinyBS             -> archived on CRAN; collapsibles are now bslib accordions
#   shinythemes, shinydashboard, shinyWidgets, shinyjs, rlist, tidyselect, rsconnect

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(DT)
  library(data.table)
  library(matrixStats)
  library(stringi)
})

# ---------------------------------------------------------------------------
# Runtime options
# ---------------------------------------------------------------------------

# v1 capped uploads at 30 MB, which silently rejected the 25 MB MHCFlurry example
# plus its FASTA. Sized here for the ~100-300 MB prediction files the app targets.
options(shiny.maxRequestSize = 400 * 1024^2)

# data.table threading: shinyapps.io containers are 1-2 cores. Using all of them
# inside a single Shiny process starves the event loop, so cap at 2.
#
# detectCores() returns NA in many containers (shinyapps.io included), and
# setDTthreads(NA) is an error -- which would stop the app from booting there.
local({
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1L) n <- 1L
  data.table::setDTthreads(max(1L, min(2L, as.integer(n))))
})

options(
  stringsAsFactors = FALSE,
  # Show the real error text in the UI instead of the generic "an error occurred";
  # every output is wrapped by ee_safe() so nothing leaks a stack trace to the user.
  shiny.sanitize.errors = FALSE
)

# ---------------------------------------------------------------------------
# Analysis-wide constants (palette, default cutoffs, rendering guards)
# ---------------------------------------------------------------------------

source("utils/constants.R", local = FALSE)

# ---------------------------------------------------------------------------
# Source order matters: utils before modules.
# ---------------------------------------------------------------------------

source("utils/parsing_functions.R", local = FALSE)
source("utils/core_functions.R",    local = FALSE)
source("utils/plot_functions.R",    local = FALSE)
source("utils/ui_helpers.R",        local = FALSE)
