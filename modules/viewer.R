# viewer.R -- Tool 4: where the epitopes sit along one protein

viewer_ui <- function(id) {
  ns <- NS(id)
  ee_tool_page(
    title = "Epitope Viewer",
    subtitle = "Positions of predicted epitopes along a single protein, coloured by how many alleles bind them.",
    sidebar_content = tagList(
      selectizeInput(ns("protein"), "Protein", choices = NULL,
                     options = list(placeholder = "Select a protein")),
      selectizeInput(ns("alleles"), "MHC alleles", choices = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button"), placeholder = "Select alleles")),
      radioButtons(ns("mode"), "Condition over alleles",
                   choices = c("Intersection", "Union"), selected = "Union"),
      numericInput(ns("cutoff"), "Score cutoff", value = 2, min = 0, step = 0.1),
      hr(),
      ee_run_button(ns("go")),
      ee_display_head(),
      checkboxInput(ns("label"), "Label the protein bar", TRUE)
    ),
    help_content = HTML(
      "<p>The blue bar is the protein, drawn to scale in amino acids. Each block underneath
       is a predicted epitope, placed at its real position and stacked into rows so that
       overlapping epitopes stay readable. Colour runs yellow to red with the number of
       alleles the epitope is predicted to bind, so dense red clusters are the promiscuous
       hotspots worth following up. Hover any block for its sequence, coordinates and
       allele list; drag to zoom into a region.</p>
       <ul>
         <li><b>Protein</b> &mdash; one at a time, from the loaded FASTA.</li>
         <li><b>Condition</b> &mdash; <b>Union</b> shows epitopes binding any selected allele;
             <b>Intersection</b> only those binding all of them.</li>
         <li><b>Score cutoff</b> &mdash; defaults to 2 %rank for class I, 10 for class II.</li>
       </ul>
       <p class='ee-note'><b>Changed in v2:</b> the row-stacking is a first-fit interval sweep
       instead of v1's pairwise overlap search, so a 7,000-residue protein with every window
       selected lays out in well under a second rather than minutes.</p>"),

    ee_stat_strip_output(ns("stats")),
    ee_card("Epitope map", ee_busy(ee_dynamic_plot(ns("plot_box")))),
    ee_card("Epitopes in this protein", download_id = ns("dl"),
            ee_busy(DT::dataTableOutput(ns("table"))))
  )
}

viewer_server <- function(id, ds) {
  moduleServer(id, function(input, output, session) {

    observeEvent(ds(), {
      d <- ds()
      ee_update_choices(session, "protein", d$proteins$ID, isolate(input$protein), 1L, "select")
      ee_update_choices(session, "alleles", d$alleles, isolate(input$alleles),
                        length(d$alleles), "select")
      ee_apply_defaults(session, d, list(cutoff = "cutoff"))
    }, ignoreNULL = TRUE)

    p <- ee_analysis_params(ds, reactive(input$go), function(d, s) list(
      protein = ee_pick(input$protein, d$proteins$ID, 1L)[1],
      alleles = ee_pick(input$alleles, d$alleles, length(d$alleles)),
      mode    = input$mode %||% "Union",
      cutoff  = ee_val(input$cutoff, s$cutoff)
    ))

    lay <- reactive({
      d <- req(ds()); q <- p()
      req(nzchar(q$protein %||% ""))
      ee_viewer_layout(d, q$protein, q$alleles, q$mode, q$cutoff)
    })

    output$stats <- renderUI({
      ee_safe({
        L <- lay(); ep <- L$epitopes
        cov <- if (nrow(ep)) {
          # Residues covered by at least one epitope, via an interval union.
          o <- ep[order(Pos)]
          hi <- cummax(o$End); starts <- o$Pos
          new <- c(TRUE, starts[-1] > hi[-length(hi)])
          grp <- cumsum(new)
          sum(vapply(split(seq_len(nrow(o)), grp),
                     function(ix) max(o$End[ix]) - min(o$Pos[ix]) + 1L, numeric(1)))
        } else 0
        tagList(
          ee_stat_strip(
            ee_stat("Protein length", paste(ee_num(L$protein$length), "aa")),
            ee_stat("Epitopes",       ee_num(nrow(ep))),
            ee_stat("Stacked rows",   ee_num(L$n_lanes)),
            ee_stat("Residues covered", sprintf("%s (%.0f%%)", ee_num(cov),
                                                100 * cov / max(L$protein$length, 1))),
            ee_stat("Max alleles bound", if (nrow(ep)) ee_num(max(ep$NAlleles)) else "0")
          ),
          if (isTRUE(L$truncated))
            ee_alert(sprintf(paste("Showing the %s most promiscuous epitopes.",
                                   "The table and the download contain every one."),
                             ee_num(EE$limits$viewer_rects)), "warning")
        )
      }, as = "ui", context = "viewer stats")
    })

    output$plot_box <- renderUI({
      n <- tryCatch(lay()$n_lanes, error = function(e) 8L)
      ee_render_dynamic_plot(session$ns("plot"),
                             ee_auto_height(n, per_row = 22, min = 260, max = 900))
    })

    output$plot <- renderPlotly({
      # Display-only: the label toggle applies without pressing Run.
      ee_safe(ee_plot_viewer(lay(), show_labels = isTRUE(input$label)), context = "viewer plot")
    })

    tbl <- reactive({
      ep <- data.table::copy(lay()$epitopes)
      if (!nrow(ep)) return(ep)
      data.table::setnames(ep[, .(Peptide, Pos, End, NAlleles, Alleles)],
                           c("Peptide", "Start", "End", "N alleles", "Alleles"))
    })

    output$table <- DT::renderDataTable({
      ee_safe(ee_datatable(tbl(), page = 10), as = "table", context = "viewer table")
    })

    output$dl <- ee_download_table("epitope_viewer", tbl)
  })
}
