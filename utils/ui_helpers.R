# ui_helpers.R -- Shared UI scaffolding and error handling for Epitope-Evaluator 2
#
# The single most common failure mode of v1 was a grey panel or a
# "Disconnected from the server" banner: an error inside renderPlotly() or
# renderDataTable() with nothing to tell the user what went wrong. Every output
# in v2 goes through ee_safe(), which turns any error into a readable message
# inside the panel and keeps the rest of the session alive.

#' NULL/NA-coalescing operator, used across the modules for unset inputs
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

#' Run an expression, converting failure into a visible in-panel message
#'
#' @param expr The rendering expression.
#' @param as "plot", "table" or "ui" -- how to express a failure.
#' @param context Short label naming the tool, prefixed to the message.
ee_safe <- function(expr, as = c("plot", "table", "ui"), context = NULL) {
  as <- match.arg(as)
  tryCatch(
    force(expr),
    # shiny::validate()/req() signal their own condition; let Shiny handle those.
    shiny.silent.error = function(e) NULL,
    validation = function(e) NULL,
    error = function(e) {
      msg <- conditionMessage(e)
      if (!is.null(context)) message(sprintf("[%s] %s", context, msg))
      switch(as,
        plot  = ee_empty_plot("Could not build this figure.", htmltools::htmlEscape(msg)),
        table = DT::datatable(
          data.frame(Problem = msg, check.names = FALSE),
          rownames = FALSE, options = list(dom = "t", ordering = FALSE),
          class = "compact"
        ),
        ui    = ee_alert(msg, "danger")
      )
    }
  )
}

#' A coloured message box
ee_alert <- function(text, type = c("info", "warning", "danger", "success"), icon = TRUE) {
  type <- match.arg(type)
  sym <- c(info = "ℹ", warning = "▲", danger = "✕", success = "✓")[[type]]
  div(class = paste0("alert alert-", type, " ee-alert"), role = "alert",
      if (icon) span(class = "ee-alert-icon", sym), span(HTML(text)))
}

#' Show a transient toast without stacking duplicates
ee_notify <- function(text, type = "warning", id = NULL, duration = 8) {
  showNotification(HTML(text), type = type, id = id, duration = duration)
}

# ---------------------------------------------------------------------------
# Layout building blocks
# ---------------------------------------------------------------------------

#' The collapsible "how this tool works" panel
#'
#' Replaces shinyBS::bsCollapsePanel. shinyBS has been archived on CRAN, so a
#' v1 deploy to a fresh server can no longer resolve it.
ee_help_panel <- function(content, title = "How this tool works", open = FALSE) {
  bslib::accordion(
    open = if (open) title else FALSE, class = "ee-help",
    bslib::accordion_panel(title = title, value = title, content)
  )
}

#' The standard tool page: parameters on the left, results on the right
#'
#' fillable/fill are FALSE deliberately. bslib's default puts the main area
#' inside an `html-fill-container` flexbox that stretches to the viewport and
#' scrolls internally, so the page does not scroll as one document and content
#' above the fold is lost inside a nested scroller. With them off the results
#' column is ordinary document flow, the window scrollbar is the only one, and
#' the sidebar is made sticky in styles.css so the parameters stay in view.
ee_tool_page <- function(title, subtitle, sidebar_content, help_content, ...) {
  bslib::layout_sidebar(
    fillable = FALSE, fill = FALSE,
    sidebar = bslib::sidebar(
      width = 300, title = "Parameters", open = "desktop",
      class = "ee-sidebar", sidebar_content
    ),
    div(
      class = "ee-tool-head",
      h4(title, class = "ee-tool-title"),
      p(subtitle, class = "ee-tool-sub")
    ),
    ee_help_panel(help_content),
    ...
  )
}

#' A results card with an optional download button in its header
ee_card <- function(title, ..., download_id = NULL, download_label = "Download table",
                    height = NULL, full_screen = TRUE) {
  bslib::card(
    full_screen = full_screen, height = height, fill = FALSE, class = "ee-card",
    bslib::card_header(
      class = "ee-card-header",
      span(title),
      if (!is.null(download_id))
        downloadButton(download_id, download_label, class = "btn-sm btn-outline-primary ee-dl")
    ),
    bslib::card_body(..., padding = 10)
  )
}

#' Run-analysis button, styled consistently across the six tools
ee_run_button <- function(id, label = "Run analysis") {
  actionButton(id, label, icon = icon("play"), class = "btn-primary w-100 ee-run")
}

#' A numeric input that states its unit and range in the label
ee_numeric <- function(id, label, value, min = NA, max = NA, step = NA, help = NULL) {
  tagList(
    numericInput(id, label, value = value, min = min, max = max, step = step),
    if (!is.null(help)) div(class = "ee-hint", help)
  )
}

#' Spinner wrapper -- one place to change the busy indicator
ee_busy <- function(output_tag) {
  div(class = "ee-busy", output_tag)
}

#' Standard DT options: server-side paging, horizontal scroll, CSV/TSV buttons off
#' (downloads go through downloadHandler so the full, uncapped table is written).
ee_datatable <- function(df, ..., page = 15L, scrollY = NULL) {
  n <- nrow(df)
  capped <- n > EE$limits$table_rows
  if (capped) df <- utils::head(df, EE$limits$table_rows)
  DT::datatable(
    df,
    rownames = FALSE, selection = "multiple", class = "compact stripe hover",
    caption = if (capped) htmltools::tags$caption(
      style = "caption-side:top; font-size:12px; color:#a05000;",
      sprintf("Showing the first %s of %s rows. The download button gives every row.",
              format(EE$limits$table_rows, big.mark = ","), format(n, big.mark = ","))
    ) else NULL,
    options = list(
      pageLength = page, lengthMenu = c(10, 15, 25, 50, 100),
      scrollX = TRUE, scrollY = scrollY, autoWidth = FALSE,
      dom = "<'row'<'col-sm-6'l><'col-sm-6'f>>tip",
      search = list(regex = FALSE, caseInsensitive = TRUE)
    ),
    ...
  )
}

#' A downloadHandler that writes a tab-separated table and never fails silently
#'
#' Downloads always contain the complete result, even where the on-screen table
#' or figure was capped for rendering. A failing query writes a one-row
#' explanation rather than a zero-byte file, which is what v1 produced whenever
#' the underlying reactive threw.
#'
#' The two closures are also attached as attributes so tests can exercise them
#' without reaching into Shiny's render-function internals.
ee_download_table <- function(name, data_fun) {
  fname <- function() sprintf("%s_%s.tsv", name, format(Sys.Date(), "%Y%m%d"))
  writer <- function(file) {
    d <- tryCatch(data_fun(), error = function(e)
      # Flatten the message: a tab or newline inside it would corrupt the TSV.
      data.frame(Error = gsub("[\t\r\n]+", " ", conditionMessage(e))))
    if (is.null(d) || !NROW(d)) d <- data.frame(Note = "No rows matched the current settings.")
    data.table::fwrite(d, file, sep = "\t", quote = FALSE)
  }
  structure(downloadHandler(filename = fname, content = writer),
            ee_filename = fname, ee_content = writer)
}

# ---------------------------------------------------------------------------
# Parameter synchronisation
# ---------------------------------------------------------------------------

#' Push dataset-derived defaults into a tool's inputs
#'
#' This is the fix for one of v1's quietly broken features: data_input.R called
#' updateNumericInput() with the *data-input module's* session, so the IDs
#' resolved to "data_input-cutoff" and never reached the tools. Class II users
#' therefore always got class I cutoffs. Here each module updates its own inputs
#' from its own session, which is the only session that can namespace them
#' correctly.
#'
#' @param session The calling module's session.
#' @param ds An ee_dataset.
#' @param numeric_map Named list: input id -> field name in ee_suggest(ds).
ee_apply_defaults <- function(session, ds, numeric_map) {
  s <- ee_suggest(ds)
  hi <- max(s$score_max, s$cutoff, na.rm = TRUE)
  for (id in names(numeric_map)) {
    field <- numeric_map[[id]]
    v <- s[[field]]
    if (is.null(v) || !is.finite(v)) next
    updateNumericInput(session, id, value = v, min = 0,
                       max = if (identical(field, "step")) hi else hi,
                       step = if (identical(field, "step")) NULL else s$step)
  }
}

#' Resolve a tool's analysis parameters
#'
#' One consistent rule across all six tools:
#'
#'   * ANALYSIS parameters (alleles, proteins, cutoffs, AND/OR) are read here,
#'     isolated, so they only take effect when Run is pressed -- changing them
#'     never silently recomputes a large analysis under the user.
#'   * DISPLAY parameters (plot type, sort order, log scale, labels) are read
#'     directly inside the render functions, so they apply immediately.
#'
#' The block also re-runs whenever a new dataset is loaded, so a tool tab is
#' never blank on arrival: it shows the analysis at the dataset's own defaults.
#' At that moment the client has not yet echoed back the updated selectors, so
#' every value falls back to something derived from the dataset, and stale
#' selections left over from a previous dataset are filtered out.
#'
#' @param ds Reactive returning the ee_dataset.
#' @param go The Run button's input value.
#' @param fn function(d, s) returning the parameter list; `d` is the dataset and
#'   `s` is ee_suggest(d). Called inside isolate().
ee_analysis_params <- function(ds, go, fn) {
  reactive({
    d <- req(ds())
    force(go())                      # re-run on Run
    s <- ee_suggest(d)
    isolate(fn(d, s))
  })
}

#' Keep only selections that exist in this dataset, else fall back
ee_pick <- function(chosen, available, n_default = 1L) {
  keep <- intersect(chosen, available)
  if (length(keep)) keep else utils::head(available, n_default)
}

#' A numeric input value, falling back when the client has not sent one yet
ee_val <- function(x, fallback) {
  if (is.null(x) || length(x) != 1L || !is.finite(x)) fallback else x
}

#' Render the first rows of a real prediction file as an HTML table
#'
#' v1's Documentation showed a live DT of each bundled example, so you could
#' see what the file actually looks like rather than read a description of it.
#' That is the more useful reference, so it is back -- but built statically at
#' UI time (the files never change) instead of through a renderDataTable round
#' trip on every session.
#'
#' The NetMHC family's two header rows are shown as header rows, because that
#' two-row structure is the single most confusing thing about the format and
#' the thing users most often destroy by re-saving through Excel.
#'
#' @param path File to preview.
#' @param sep Field separator.
#' @param header_rows How many leading lines are header (2 for the NetMHC family).
#' @param n_rows,n_cols How much to show.
ee_file_preview <- function(path, sep = "\t", header_rows = 1L, n_rows = 4L, n_cols = 14L) {
  if (!file.exists(path)) {
    return(div(class = "ee-hint", sprintf("[%s not bundled]", basename(path))))
  }
  lines <- tryCatch(readLines(path, n = header_rows + n_rows, warn = FALSE),
                    error = function(e) character(0))
  if (!length(lines)) return(div(class = "ee-hint", "[could not read the file]"))

  cells <- lapply(strsplit(lines, sep, fixed = TRUE), function(v) {
    v <- if (length(v) > n_cols) c(v[seq_len(n_cols)], "\u2026") else v
    ifelse(is.na(v) | !nzchar(trimws(v)), "", trimws(v))
  })
  width <- max(lengths(cells))
  cells <- lapply(cells, function(v) c(v, rep("", width - length(v))))

  hdr <- cells[seq_len(min(header_rows, length(cells)))]
  bdy <- if (length(cells) > header_rows) cells[-seq_len(header_rows)] else list()

  div(
    class = "ee-preview-wrap",
    tags$table(
      class = "ee-preview",
      tags$thead(lapply(hdr, function(r) tags$tr(lapply(r, tags$th)))),
      tags$tbody(lapply(bdy, function(r) tags$tr(lapply(r, tags$td))))
    ),
    div(class = "ee-figcaption",
        sprintf("%s \u2014 first %d data row(s)%s.", basename(path), length(bdy),
                if (width > n_cols) sprintf(", first %d of many columns", n_cols) else ""))
  )
}

#' A captioned screenshot, degrading gracefully when the file is absent
#'
#' Every figure on the static pages is a real capture of THIS build, produced by
#' tools/screenshots.py. If one is missing (a source checkout that skipped the
#' images, say) the caption is still shown rather than a broken image icon --
#' which is what v1 did for eight tutorial videos that were never committed.
#'
#' @param file Filename inside www/images.
#' @param caption Text shown beneath the image.
#' @param href Optional deep link that opens the live view being illustrated.
ee_figure <- function(file, caption, href = NULL) {
  path <- file.path("www", "images", file)
  if (!file.exists(path)) {
    return(div(class = "ee-figure ee-figure-missing",
               div(class = "ee-hint", sprintf("[screenshot %s not bundled]", file)),
               div(class = "ee-figcaption", caption)))
  }
  img <- tags$img(src = file.path("images", file), alt = caption,
                  class = "ee-shot", loading = "lazy")
  tags$figure(
    class = "ee-figure",
    if (is.null(href)) img else tags$a(href = href, title = "Open this view", img),
    tags$figcaption(class = "ee-figcaption", caption)
  )
}

#' Sidebar heading marking the controls that apply without pressing Run
ee_display_head <- function() {
  tagList(hr(), div(class = "ee-sub-head", "Display"),
          div(class = "ee-hint", "Applies immediately \u2014 no need to press Run."))
}

#' Refresh an allele or protein selector, preserving the user's choice when possible
ee_update_choices <- function(session, id, choices, current, default_n = 1L,
                              type = c("checkbox", "select", "picker")) {
  type <- match.arg(type)
  keep <- intersect(current, choices)
  sel  <- if (length(keep)) keep else utils::head(choices, default_n)
  switch(type,
    checkbox = updateCheckboxGroupInput(session, id, choices = choices, selected = sel),
    select   = updateSelectizeInput(session, id, choices = choices, selected = sel, server = TRUE),
    picker   = updateSelectInput(session, id, choices = choices, selected = sel)
  )
}

# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

ee_num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

#' Small "N epitopes / N proteins" summary strip shown above results
ee_stat_strip <- function(...) {
  items <- list(...)
  div(class = "ee-stats",
      lapply(items, function(it) div(class = "ee-stat",
        div(class = "ee-stat-value", it$value),
        div(class = "ee-stat-label", it$label))))
}
ee_stat <- function(label, value) list(label = label, value = value)

#' Placeholder for the summary strip rendered by a module's server
ee_stat_strip_output <- function(id) uiOutput(id, class = "ee-stats-slot")

#' A plotly output whose height is computed server-side
#'
#' Heatmaps need to grow with their row count. Rather than shipping custom
#' JavaScript to resize a fixed container, the module renders the output tag
#' itself with the height it needs.
ee_dynamic_plot <- function(id) uiOutput(id)

ee_render_dynamic_plot <- function(id, height) {
  plotlyOutput(id, height = paste0(round(height), "px"))
}
