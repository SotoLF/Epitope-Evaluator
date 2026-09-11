# plot_functions.R -- Native plotly figures for Epitope-Evaluator 2
#
# v1 built most figures in ggplot2 and pushed them through ggplotly(). That
# round-trip re-walks the whole grob tree on every redraw and is the single
# slowest step once a heatmap has more than a few thousand tiles. Everything
# here is a plotly trace built directly from a data.table, which also removes
# the ggplot2, gridExtra, sf and ggVennDiagram dependencies.
#
# Two rules every figure follows:
#   * aggregate server-side; never ship one browser object per data point
#   * always return a figure, never NULL -- an empty state is a message, not a
#     blank panel or a disconnected session

# ---------------------------------------------------------------------------
# Shared styling
# ---------------------------------------------------------------------------

ee_font <- list(family = "system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif",
                size = 13, color = EE$col$text)

#' Common layout applied to every figure
ee_layout <- function(p, xlab = "", ylab = "", title = NULL, legend = FALSE, ...) {
  plotly::layout(
    p,
    title  = if (is.null(title)) NULL else list(text = title, font = c(ee_font, list(size = 15)), x = 0),
    xaxis  = list(title = list(text = xlab, font = ee_font), tickfont = ee_font,
                  gridcolor = EE$col$grid, zerolinecolor = EE$col$grid, automargin = TRUE),
    yaxis  = list(title = list(text = ylab, font = ee_font), tickfont = ee_font,
                  gridcolor = EE$col$grid, zerolinecolor = EE$col$grid, automargin = TRUE),
    font        = ee_font,
    showlegend  = legend,
    hoverlabel  = list(font = ee_font, bgcolor = "white", bordercolor = EE$col$grid),
    plot_bgcolor  = "white",
    paper_bgcolor = "white",
    margin = list(l = 60, r = 20, t = if (is.null(title)) 20 else 46, b = 50),
    ...
  )
}

#' Strip the plotly toolbar down to what is actually useful, and make the PNG
#' export high-resolution (v1 exported at screen resolution).
ee_config <- function(p, filename = "epitope-evaluator") {
  plotly::config(
    p,
    displaylogo = FALSE,
    modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverClosestCartesian",
                               "hoverCompareCartesian", "toggleSpikelines"),
    toImageButtonOptions = list(format = "png", filename = filename,
                                width = 1600, height = 1000, scale = 2)
  )
}

#' A figure that says why it is empty
#'
#' Returned instead of NULL wherever a filter yields nothing. v1 returned NULL
#' from several renderPlotly() blocks, which leaves a grey panel with no
#' explanation, or errors outright in the Viewer.
ee_empty_plot <- function(message = "No data to display.", hint = NULL) {
  txt <- if (is.null(hint)) message else paste0(message, "<br><span style='font-size:12px'>", hint, "</span>")
  plotly::plot_ly(type = "scatter", mode = "markers", x = numeric(0), y = numeric(0)) |>
    plotly::layout(
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
      plot_bgcolor = "white", paper_bgcolor = "white",
      annotations = list(list(text = txt, showarrow = FALSE, xref = "paper", yref = "paper",
                              x = 0.5, y = 0.5, font = c(ee_font, list(size = 14, color = "#777"))))
    ) |>
    ee_config()
}

#' Shorten long labels for axis ticks while keeping them unique
ee_short <- function(x, n = 28L) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)
}

#' A height that grows with the number of rows, within sane bounds
ee_auto_height <- function(n_rows, per_row = 18, min = 320, max = 1400) {
  max(min, min(max, 120 + per_row * n_rows))
}

# ---------------------------------------------------------------------------
# Tool 1 -- Distribution
# ---------------------------------------------------------------------------

#' Histogram or cumulative histogram of combined scores
#'
#' Bars come from ee_histogram(), so the browser receives a few hundred numbers
#' regardless of how many peptides were analysed.
ee_plot_histogram <- function(h, cumulative = FALSE, xlab = "% rank", show_density = TRUE) {
  if (!nrow(h) || sum(h$count) == 0L) {
    return(ee_empty_plot("No epitopes in this range.",
                         "Widen the %rank window or select more alleles."))
  }
  y     <- if (cumulative) h$cumulative else h$count
  ylab  <- if (cumulative) "Cumulative number of epitopes" else "Number of epitopes"
  width <- h$bin_end - h$bin_start

  p <- plotly::plot_ly() |>
    plotly::add_bars(
      x = h$bin_mid, y = y, width = width * 0.96,
      marker = list(color = EE$col$accent_soft,
                    line = list(color = EE$col$accent, width = 1)),
      hovertext = sprintf("<b>%.3g – %.3g</b><br>%s: %s",
                          h$bin_start, h$bin_end, ylab, format(y, big.mark = ",")),
      hovertemplate = "%{hovertext}<extra></extra>",
      name = "Epitopes"
    )

  # Density overlay only makes sense on the non-cumulative view.
  if (show_density && !cumulative && sum(h$count) > 0) {
    p <- plotly::add_lines(
      p, x = h$bin_mid, y = h$density * max(y) / max(h$density),
      line = list(color = EE$col$highlight, width = 2, shape = "spline"),
      yaxis = "y2", hoverinfo = "skip", name = "Density"
    ) |>
      plotly::layout(yaxis2 = list(overlaying = "y", side = "right", showgrid = FALSE,
                                   showticklabels = FALSE, rangemode = "tozero"))
  }
  p |> ee_layout(xlab = xlab, ylab = ylab) |> ee_config("epitope-distribution")
}

#' Allele x cutoff heatmap of epitope counts
ee_plot_cutoff_heatmap <- function(grid, xlab = "% rank cutoff") {
  if (!nrow(grid)) return(ee_empty_plot("No alleles selected."))

  alleles <- levels(grid$allele)
  cutoffs <- sort(unique(grid$cutoff))
  z <- matrix(grid$count[order(grid$allele, grid$cutoff)],
              nrow = length(alleles), ncol = length(cutoffs), byrow = TRUE)

  # Annotate cell values only when the grid is small enough to read them.
  annotate <- length(cutoffs) * length(alleles) <= 300

  p <- plotly::plot_ly(
    x = cutoffs, y = ee_short(alleles), z = z, type = "heatmap",
    colorscale = EE$ramp, xgap = 1, ygap = 1,
    colorbar = list(title = list(text = "Epitopes", font = ee_font), tickfont = ee_font,
                    thickness = 12, len = 0.7),
    hovertemplate = "Allele: <b>%{y}</b><br>Cutoff ≤ %{x}<br>Epitopes: %{z:,}<extra></extra>"
  )
  if (annotate) {
    p <- plotly::add_annotations(
      p, x = rep(cutoffs, each = length(alleles)), y = rep(ee_short(alleles), length(cutoffs)),
      text = format(as.vector(z), big.mark = ","), showarrow = FALSE,
      font = list(size = 10, family = ee_font$family,
                  color = ifelse(as.vector(z) > max(z, na.rm = TRUE) * 0.6, "white", EE$col$text))
    )
  }
  p |> ee_layout(xlab = xlab, ylab = "") |> ee_config("epitopes-by-cutoff")
}

# ---------------------------------------------------------------------------
# Tools 2 & 6 -- UpSet and Venn
# ---------------------------------------------------------------------------

#' UpSet plot, built natively in plotly
#'
#' Three panels sharing coordinates: intersection sizes on top, the dot matrix
#' below it, and set sizes to the left. Unlike v1 this never enumerates 2^k
#' combinations -- it renders exactly the combinations ee_combinations() found.
ee_plot_upset <- function(cmb, sizes, max_bars = EE$limits$upset_max_bars,
                          set_label = "Set") {
  if (!nrow(cmb)) return(ee_empty_plot("No shared epitopes at this cutoff."))

  cmb <- utils::head(cmb, max_bars)
  sets <- sizes[order(sizes$size), ]           # smallest at the top of the matrix
  set_levels <- sets$set
  ns <- length(set_levels)
  ni <- nrow(cmb)

  members <- strsplit(cmb$sets, " & ", fixed = TRUE)
  dot_x <- rep(seq_len(ni), lengths(members))
  dot_y <- match(unlist(members), set_levels)

  # --- top: intersection sizes
  top <- plotly::plot_ly() |>
    plotly::add_bars(
      x = seq_len(ni), y = cmb$count,
      marker = list(color = EE$col$accent),
      text = format(cmb$count, big.mark = ","), textposition = "outside",
      textfont = list(size = 10, family = ee_font$family),
      hovertext = sprintf("<b>%s</b><br>Epitopes: %s",
                          cmb$sets, format(cmb$count, big.mark = ",")),
      hovertemplate = "%{hovertext}<extra></extra>"
    ) |>
    plotly::layout(
      xaxis = list(range = c(0.4, ni + 0.6), showticklabels = FALSE, showgrid = FALSE,
                   zeroline = FALSE, fixedrange = TRUE),
      yaxis = list(title = list(text = "Intersection size", font = ee_font),
                   tickfont = ee_font, gridcolor = EE$col$grid, rangemode = "tozero")
    )

  # --- bottom right: the dot matrix, with connecting lines per intersection
  seg_x <- unlist(lapply(seq_len(ni), function(i) c(i, i, NA)))
  seg_y <- unlist(lapply(members, function(mm) {
    r <- range(match(mm, set_levels)); c(r[1], r[2], NA)
  }))

  matrixp <- plotly::plot_ly() |>
    # background dots for every cell
    plotly::add_markers(
      x = rep(seq_len(ni), each = ns), y = rep(seq_len(ns), ni),
      marker = list(color = "#E8E8E8", size = 9), hoverinfo = "skip", showlegend = FALSE
    ) |>
    plotly::add_segments(
      x = seg_x[c(TRUE, FALSE, FALSE)], xend = seg_x[c(FALSE, TRUE, FALSE)],
      y = seg_y[c(TRUE, FALSE, FALSE)], yend = seg_y[c(FALSE, TRUE, FALSE)],
      line = list(color = EE$col$accent, width = 2), hoverinfo = "skip", showlegend = FALSE
    ) |>
    plotly::add_markers(
      x = dot_x, y = dot_y,
      marker = list(color = EE$col$accent, size = 10),
      hovertext = cmb$sets[dot_x],
      hovertemplate = "<b>%{hovertext}</b><extra></extra>", showlegend = FALSE
    ) |>
    plotly::layout(
      xaxis = list(range = c(0.4, ni + 0.6), showticklabels = FALSE, showgrid = FALSE,
                   zeroline = FALSE, fixedrange = TRUE),
      yaxis = list(range = c(0.4, ns + 0.6), tickvals = seq_len(ns),
                   ticktext = ee_short(set_levels), tickfont = ee_font,
                   showgrid = FALSE, zeroline = FALSE, fixedrange = TRUE)
    )

  # --- bottom left: set sizes, pointing left so they read into the matrix
  size_ticks <- pretty(c(0, max(sets$size)), 4)
  size_ticks <- size_ticks[size_ticks >= 0 & size_ticks <= max(sets$size) * 1.05]
  left <- plotly::plot_ly() |>
    plotly::add_bars(
      x = -sets$size, y = seq_len(ns), orientation = "h",
      marker = list(color = "#7F7F7F"),
      hovertext = sprintf("<b>%s</b><br>Total epitopes: %s",
                          sets$set, format(sets$size, big.mark = ",")),
      hovertemplate = "%{hovertext}<extra></extra>"
    ) |>
    plotly::layout(
      # The bars are drawn at negative x so they grow leftwards into the dot
      # matrix; the tick labels must still read as positive set sizes.
      xaxis = list(title = list(text = paste(set_label, "size"), font = ee_font),
                   tickfont = ee_font, gridcolor = EE$col$grid,
                   tickmode = "array", tickvals = -size_ticks,
                   ticktext = format(size_ticks, big.mark = ",", trim = TRUE)),
      yaxis = list(range = c(0.4, ns + 0.6), showticklabels = FALSE,
                   showgrid = FALSE, zeroline = FALSE)
    )

  blank <- plotly::plotly_empty(type = "scatter", mode = "markers")

  plotly::subplot(
    blank, top, left, matrixp,
    nrows = 2, widths = c(0.26, 0.74), heights = c(0.55, 0.45),
    shareX = FALSE, shareY = FALSE, margin = c(0.02, 0.01, 0.03, 0.02), titleY = TRUE, titleX = TRUE
  ) |>
    plotly::layout(showlegend = FALSE, font = ee_font,
                   plot_bgcolor = "white", paper_bgcolor = "white",
                   margin = list(l = 40, r = 20, t = 24, b = 50)) |>
    ee_config("upset-plot")
}

#' Geometry for a 2-, 3- or 4-set Venn diagram
#'
#' Circles for 2 and 3 sets, the classic four-ellipse layout for 4. Region label
#' anchors are found by rasterising the shapes on a grid and taking the centroid
#' of each region, which works for any layout without hard-coded coordinates.
ee_venn_geometry <- function(n) {
  if (n == 2L) {
    list(cx = c(0.36, 0.64), cy = c(0.5, 0.5), rx = c(0.30, 0.30), ry = c(0.30, 0.30),
         rot = c(0, 0), lx = c(0.14, 0.86), ly = c(0.86, 0.86))
  } else if (n == 3L) {
    list(cx = c(0.38, 0.62, 0.50), cy = c(0.58, 0.58, 0.36),
         rx = rep(0.29, 3), ry = rep(0.29, 3), rot = c(0, 0, 0),
         lx = c(0.15, 0.85, 0.50), ly = c(0.90, 0.90, 0.03))
  } else {
    list(cx = c(0.36, 0.45, 0.55, 0.64), cy = c(0.44, 0.53, 0.53, 0.44),
         rx = rep(0.34, 4), ry = rep(0.19, 4), rot = c(-40, -40, 40, 40),
         lx = c(0.06, 0.24, 0.76, 0.94), ly = c(0.60, 0.92, 0.92, 0.60))
  }
}

#' Venn / Euler diagram drawn with plotly shapes
#'
#' Replaces ggVennDiagram, whose sf dependency needs GDAL, GEOS and PROJ on the
#' host -- the main reason a v1 deploy could fail on a clean server.
ee_plot_venn <- function(mem, cmb, max_sets = EE$limits$venn_max_sets) {
  sets <- colnames(mem$m)
  n <- length(sets)
  if (n < 2L || n > max_sets) {
    return(ee_empty_plot(sprintf("A Venn diagram is only drawn for 2–%d sets.", max_sets),
                         "Switch the plot type to UpSet, which handles any number of sets."))
  }
  if (!nrow(mem$m)) return(ee_empty_plot("No epitopes pass the cutoff."))

  g <- ee_venn_geometry(n)
  th <- seq(0, 2 * pi, length.out = 121)

  # SVG path for each (possibly rotated) ellipse.
  paths <- vapply(seq_len(n), function(i) {
    a <- g$rot[i] * pi / 180
    x <- g$cx[i] + g$rx[i] * cos(th) * cos(a) - g$ry[i] * sin(th) * sin(a)
    y <- g$cy[i] + g$rx[i] * cos(th) * sin(a) + g$ry[i] * sin(th) * cos(a)
    paste0("M ", paste(sprintf("%.4f %.4f", x, y), collapse = " L "), " Z")
  }, character(1))

  pal <- c(EE$col$accent, EE$col$strong, "#009E73", "#CC79A7")
  shapes <- lapply(seq_len(n), function(i) list(
    type = "path", path = paths[i], layer = "below",
    fillcolor = paste0(pal[i], "33"),
    line = list(color = pal[i], width = 2)
  ))

  # Rasterise to locate each region's centroid.
  gr <- 260L
  gx <- rep(seq(0, 1, length.out = gr), times = gr)
  gy <- rep(seq(0, 1, length.out = gr), each  = gr)
  inside <- vapply(seq_len(n), function(i) {
    a <- -g$rot[i] * pi / 180
    dx <- gx - g$cx[i]; dy <- gy - g$cy[i]
    xr <- dx * cos(a) - dy * sin(a); yr <- dx * sin(a) + dy * cos(a)
    (xr / g$rx[i])^2 + (yr / g$ry[i])^2 <= 1
  }, logical(gr * gr))
  dim(inside) <- c(gr * gr, n)

  key <- apply(inside, 1L, function(r) paste(sets[r], collapse = " & "))
  counts <- setNames(cmb$count, cmb$sets)

  ann <- list()
  for (k in unique(key[nzchar(key)])) {
    sel <- which(key == k)
    ann[[length(ann) + 1L]] <- list(
      x = mean(gx[sel]), y = mean(gy[sel]),
      text = format(if (is.na(counts[k])) 0L else counts[[k]], big.mark = ","),
      showarrow = FALSE, font = c(ee_font, list(size = 13))
    )
  }
  for (i in seq_len(n)) {
    ann[[length(ann) + 1L]] <- list(
      x = g$lx[i], y = g$ly[i], text = paste0("<b>", ee_short(sets[i], 22), "</b>"),
      showarrow = FALSE, font = c(ee_font, list(size = 12, color = pal[i]))
    )
  }

  plotly::plot_ly(type = "scatter", mode = "markers", x = numeric(0), y = numeric(0)) |>
    plotly::layout(
      shapes = shapes, annotations = ann,
      xaxis = list(visible = FALSE, range = c(-0.05, 1.05), fixedrange = TRUE),
      yaxis = list(visible = FALSE, range = c(-0.05, 1.05), fixedrange = TRUE,
                   scaleanchor = "x", scaleratio = 1),
      plot_bgcolor = "white", paper_bgcolor = "white", showlegend = FALSE,
      margin = list(l = 10, r = 10, t = 10, b = 10)
    ) |>
    ee_config("venn-diagram")
}

# ---------------------------------------------------------------------------
# Tool 3 -- Density
# ---------------------------------------------------------------------------

#' Epitope count vs protein length, with table-driven highlighting
ee_plot_density_scatter <- function(dens, selected = character(0), log_axes = FALSE) {
  if (!nrow(dens)) return(ee_empty_plot("No proteins to plot."))

  # v1 drew the highlighted points on a categorical axis and the rest on a
  # numeric one, which silently broke the x scale whenever a row was selected.
  # One trace, one scale; selection is expressed through marker styling.
  is_sel <- dens$ID %in% selected
  txt <- sprintf(
    "<b>%s</b><br>Length: %s aa<br>Epitopes: %s<br>Density: %.4f per aa",
    dens$ID, format(dens$ProtLength, big.mark = ","),
    format(dens$Epitopes, big.mark = ","), dens$Density
  )
  mode <- if (nrow(dens) > EE$limits$scatter_webgl) "scattergl" else "scatter"

  p <- plotly::plot_ly(type = mode, mode = "markers") |>
    plotly::add_trace(
      x = dens$ProtLength, y = dens$Epitopes,
      marker = list(
        size = ifelse(is_sel, 13, 9),
        color = ifelse(is_sel, EE$col$highlight, EE$col$accent_soft),
        opacity = ifelse(is_sel, 1, 0.8),
        line = list(color = ifelse(is_sel, "#000000", EE$col$accent), width = ifelse(is_sel, 2, 1))
      ),
      hovertext = txt, hovertemplate = "%{hovertext}<extra></extra>", showlegend = FALSE
    )

  ax <- if (log_axes) "log" else "linear"
  p |>
    ee_layout(xlab = "Protein length (amino acids)", ylab = "Number of epitopes") |>
    plotly::layout(xaxis = list(type = ax, rangemode = "tozero"),
                   yaxis = list(type = ax, rangemode = "tozero")) |>
    ee_config("epitope-density")
}

#' Protein x allele heatmap or grouped bar plot
ee_plot_protein_allele <- function(pa, fill = c("Number of epitopes", "Density"),
                                   plot_type = c("Heatmap", "Bar plot"),
                                   sort_by = c("Total", "Input order"),
                                   log_scale = TRUE,
                                   max_rows = EE$limits$heatmap_rows) {
  fill <- match.arg(fill); plot_type <- match.arg(plot_type); sort_by <- match.arg(sort_by)
  if (!nrow(pa)) return(ee_empty_plot("No proteins or alleles selected."))

  pa <- data.table::copy(pa)
  pa[, value := if (fill == "Density") Density else as.numeric(Epitopes)]

  prot_order <- if (sort_by == "Total") {
    pa[, .(t = sum(value, na.rm = TRUE)), by = ID][order(-t), ID]
  } else unique(pa$ID)

  truncated <- length(prot_order) > max_rows
  if (truncated) prot_order <- prot_order[seq_len(max_rows)]
  pa <- pa[ID %chin% prot_order]

  alleles <- levels(pa$allele)
  # One pre-formatted string per point. Indexed customdata (%{customdata[0]}...)
  # does not survive plotly's R->JSON serialisation: the field is dropped and
  # the template renders literally. tests/test_app.R asserts every hover
  # reference resolves to a field that is actually present.
  hov <- function(d) sprintf(
    "Protein: <b>%s</b><br>Allele: <b>%s</b><br>Epitopes: %s<br>Length: %s aa<br>Density: %.5f",
    d$ID, as.character(d$allele), format(d$Epitopes, big.mark = ","),
    format(d$ProtLength, big.mark = ","), d$Density)

  if (plot_type == "Bar plot") {
    p <- plotly::plot_ly()
    for (pr in prot_order) {
      d <- pa[ID == pr][order(match(allele, alleles))]
      p <- plotly::add_bars(p, x = as.character(d$allele), y = d$value, name = ee_short(pr, 24),
                            hovertext = hov(d), hovertemplate = "%{hovertext}<extra></extra>")
    }
    return(
      p |>
        plotly::layout(barmode = "group",
                       xaxis = list(categoryorder = "array", categoryarray = alleles, tickangle = -45)) |>
        ee_layout(xlab = "", ylab = fill, legend = TRUE) |>
        plotly::layout(legend = list(font = ee_font, orientation = "v", x = 1.01, y = 1)) |>
        ee_config("epitopes-by-protein-allele")
    )
  }

  # Heatmap. log10 keeps a handful of huge proteins from flattening everything
  # else, exactly as in v1, but the colourbar is labelled with real counts.
  z <- matrix(NA_real_, length(prot_order), length(alleles),
              dimnames = list(prot_order, alleles))
  z[cbind(match(pa$ID, prot_order), match(as.character(pa$allele), alleles))] <- pa$value
  zplot <- if (log_scale) log10(pmax(z, 0) + if (fill == "Density") 1e-6 else 1) else z

  rng <- range(zplot, na.rm = TRUE)
  ticks <- pretty(rng, 4)
  ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
  lab <- if (log_scale) {
    # if/else, not ifelse(): the condition is a scalar, so ifelse() would
    # collapse the whole label vector to a single element and plotly would
    # silently fall back to showing the raw log10 tick values.
    v <- 10^ticks - if (fill == "Density") 1e-6 else 1
    if (fill == "Density") formatC(pmax(v, 0), format = "g", digits = 2)
    else format(round(pmax(v, 0)), big.mark = ",", trim = TRUE)
  } else format(signif(ticks, 3), trim = TRUE)

  # Heatmap hover text is a matrix aligned to z.
  htxt <- matrix("", length(prot_order), length(alleles))
  htxt[cbind(match(pa$ID, prot_order), match(as.character(pa$allele), alleles))] <- hov(pa)

  plotly::plot_ly(
    x = alleles, y = ee_short(prot_order), z = zplot, type = "heatmap",
    colorscale = EE$ramp, xgap = 1, ygap = 1,
    text = htxt, hovertemplate = "%{text}<extra></extra>",
    # tickmode = "array" is required, otherwise plotly ignores ticktext and the
    # colourbar shows the raw log10 values under a "Number of epitopes" title.
    colorbar = list(title = list(text = fill, font = ee_font), tickfont = ee_font,
                    tickmode = "array", tickvals = ticks, ticktext = lab,
                    thickness = 12, len = 0.7)
  ) |>
    plotly::layout(xaxis = list(tickangle = -45),
                   yaxis = list(categoryorder = "array", categoryarray = rev(ee_short(prot_order)))) |>
    ee_layout(xlab = "", ylab = "") |>
    ee_config("epitopes-by-protein-allele")
}

# ---------------------------------------------------------------------------
# Tool 4 -- Viewer
# ---------------------------------------------------------------------------

#' Epitope map along one protein
#'
#' Epitopes are grouped by the number of alleles they bind and each group is
#' drawn as ONE filled trace containing many NA-separated quadrilaterals. That
#' keeps the browser object count at (number of distinct allele counts + 2)
#' instead of one shape per epitope, which is what let v1's Viewer lock the tab
#' on proteins with thousands of epitopes.
ee_plot_viewer <- function(layout, show_labels = TRUE) {
  ep   <- layout$epitopes
  plen <- layout$protein$length
  if (!nrow(ep)) {
    return(ee_empty_plot("No epitopes in this protein at the current settings.",
                         "Raise the cutoff, or switch the allele condition to Union."))
  }

  h <- 0.34
  quad <- function(x0, x1, y) list(
    x = as.vector(rbind(x0, x1, x1, x0, NA_real_)),
    y = as.vector(rbind(y - h, y - h, y + h, y + h, NA_real_))
  )

  # Protein backbone at lane 0.
  bb <- quad(1, plen, 0)
  p <- plotly::plot_ly() |>
    plotly::add_polygons(
      x = bb$x, y = bb$y, fill = "toself",
      fillcolor = EE$col$protein, line = list(color = "#5B8DB8", width = 1),
      hoverinfo = "skip", showlegend = FALSE
    )

  ncols <- max(ep$NAlleles)
  ramp <- grDevices::colorRampPalette(c("#FFE08A", "#F08C3C", "#B2182B"))(max(ncols, 2L))

  for (k in sort(unique(ep$NAlleles))) {
    d <- ep[NAlleles == k]
    y <- -d$lane
    q <- quad(d$Pos, d$End + 1, y)
    p <- plotly::add_polygons(
      p, x = q$x, y = q$y, fill = "toself",
      fillcolor = ramp[k], line = list(color = "#00000055", width = 0.5),
      hoverinfo = "skip", showlegend = FALSE
    )
    # A marker layer carries the tooltip: polygons cannot hold per-shape hover.
    p <- plotly::add_markers(
      p, x = (d$Pos + d$End + 1) / 2, y = y,
      marker = list(size = 1, color = "rgba(0,0,0,0)"),
      hovertext = sprintf("<b>%s</b><br>Position: %d–%d<br>Binds %d allele(s)<br>%s",
                          d$Peptide, d$Pos, d$End, d$NAlleles, ee_short(d$Alleles, 90)),
      hovertemplate = "%{hovertext}<extra></extra>",
      showlegend = FALSE
    )
  }

  ann <- if (show_labels) list(list(
    x = plen / 2, y = 0, text = paste0("<b>", ee_short(layout$protein$id, 40), "</b>"),
    showarrow = FALSE, font = c(ee_font, list(size = 12))
  )) else list()

  # A discrete colourbar for "number of alleles bound". plotly only draws a
  # colourbar for a trace that has data, so this one sits on a real coordinate
  # and is hidden with opacity rather than with NA (which plotly warns about).
  p <- plotly::add_markers(
    p, x = plen / 2, y = 0, opacity = 0,
    marker = list(color = 1, colorscale = EE$ramp_ep, cmin = 1, cmax = max(ncols, 1),
                  showscale = TRUE, size = 0.1,
                  colorbar = list(title = list(text = "Alleles<br>bound", font = ee_font),
                                  tickfont = ee_font, thickness = 12, len = 0.6,
                                  dtick = max(1, floor(ncols / 6)))),
    hoverinfo = "skip", showlegend = FALSE
  )

  p |>
    plotly::layout(
      annotations = ann,
      xaxis = list(title = list(text = "Amino acid position", font = ee_font),
                   tickfont = ee_font, range = c(-plen * 0.01, plen * 1.01),
                   gridcolor = EE$col$grid, zeroline = FALSE),
      yaxis = list(visible = FALSE, range = c(-layout$n_lanes - 0.8, 0.8), fixedrange = TRUE),
      dragmode = "zoom", plot_bgcolor = "white", paper_bgcolor = "white",
      font = ee_font, showlegend = FALSE,
      margin = list(l = 20, r = 20, t = 16, b = 46)
    ) |>
    ee_config("epitope-viewer")
}

# ---------------------------------------------------------------------------
# Tool 5 -- Promiscuity
# ---------------------------------------------------------------------------

#' Peptide x allele strong/weak binder heatmap
ee_plot_promiscuity <- function(pr, score_label = "% rank") {
  cls <- pr$matrix
  if (!nrow(cls)) {
    return(ee_empty_plot("No epitope binds that many alleles.",
                         "Lower the minimum number of alleles, or raise the weak-binding cutoff."))
  }
  alleles <- colnames(cls)
  peps <- rownames(cls)

  # Three discrete states mapped onto a stepped colourscale.
  z <- matrix(0, nrow(cls), ncol(cls))
  z[cls == "WB"] <- 1; z[cls == "SB"] <- 2
  scale <- list(
    list(0,     EE$col$none),  list(0.333, EE$col$none),
    list(0.333, EE$col$weak),  list(0.667, EE$col$weak),
    list(0.667, EE$col$strong),list(1,     EE$col$strong)
  )

  state <- ifelse(cls == "", "not a binder", ifelse(cls == "SB", "strong binder", "weak binder"))
  shown <- ifelse(is.na(pr$values), "n/a", format(signif(pr$values, 4)))
  htxt <- matrix(sprintf("<b>%s</b><br>%s<br>%s: %s<br>%s",
                         rep(peps, ncol(cls)), rep(alleles, each = nrow(cls)),
                         score_label, shown, state),
                 nrow(cls), ncol(cls))

  plotly::plot_ly(
    x = alleles, y = peps, z = z, type = "heatmap",
    colorscale = scale, zmin = 0, zmax = 2, showscale = FALSE, xgap = 1, ygap = 1,
    text = htxt, hovertemplate = "%{text}<extra></extra>"
  ) |>
    plotly::layout(
      xaxis = list(tickangle = -45),
      # Keep the promiscuity ordering from ee_promiscuity() instead of letting
      # plotly sort the categories alphabetically.
      yaxis = list(categoryorder = "array", categoryarray = rev(peps),
                   tickfont = c(ee_font, list(size = 10, family = "ui-monospace, monospace"))),
      annotations = list(
        list(x = 0, y = 1.04, xref = "paper", yref = "paper", showarrow = FALSE,
             xanchor = "left", font = c(ee_font, list(size = 12)),
             text = paste0("<span style='color:", EE$col$strong, "'>■</span> strong binder   ",
                           "<span style='color:", EE$col$weak,   "'>■</span> weak binder   ",
                           "<span style='color:#BBBBBB'>■</span> not a binder"))
      )
    ) |>
    ee_layout(xlab = "", ylab = "") |>
    ee_config("epitope-promiscuity")
}
