# constants.R -- Analysis constants shared by the engine and the Shiny layer.
#
# Kept separate from global.R so that utils/core_functions.R can be sourced and
# tested without Shiny installed (see tests/test_core.R).

EE <- local({
  e <- new.env(parent = emptyenv())

  e$version <- "2.0.0"

  # Default cutoffs by MHC class (Reynisson 2020 for class I, Jensen 2018 for class II).
  e$defaults <- list(
    I  = list(cutoff = 2,  strong = 0.5, weak = 2,  step = 0.1),
    II = list(cutoff = 10, strong = 2,   weak = 10, step = 0.5)
  )

  # Colourblind-safe palette (Okabe-Ito derived), applied consistently everywhere.
  e$col <- list(
    accent      = "#0072B2",
    accent_soft = "#7FB8DC",
    strong      = "#D55E00",   # strong binder
    weak        = "#E69F00",   # weak binder
    none        = "#F4F4F4",   # non-binder
    protein     = "#9CC3E5",   # protein backbone in the Viewer
    grid        = "#E4E4E4",
    text        = "#333333",
    highlight   = "#CC3311"
  )
  e$ramp    <- list(c(0, "#FFFFFF"), c(0.5, "#8CBFDE"), c(1, "#08519C"))
  # Viewer: number of alleles bound, low to high.
  e$ramp_ep <- list(c(0, "#FFE08A"), c(0.5, "#F08C3C"), c(1, "#B2182B"))

  # Rendering guards: past these sizes we aggregate or downsample rather than
  # shipping hundreds of thousands of SVG elements to the browser.
  e$limits <- list(
    heatmap_rows   = 300L,
    heatmap_cols   = 120L,
    viewer_rects   = 20000L,
    venn_max_sets  = 4L,
    upset_max_bars = 40L,
    table_rows     = 50000L,
    scatter_webgl  = 3000L    # switch the density scatter to WebGL past this
  )

  e
})
