# density.R -- Tool 3: epitope count vs protein length, and the protein x allele grid

density_ui <- function(id) {
  ns <- NS(id)
  ee_tool_page(
    title = "Epitope Density",
    subtitle = "How epitope count scales with protein length, and which proteins are unusually epitope-rich.",
    sidebar_content = tagList(
      selectizeInput(ns("alleles"), "MHC alleles", choices = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button"), placeholder = "Select alleles")),
      radioButtons(ns("mode"), "Condition over alleles",
                   choices = c("Intersection", "Union"), selected = "Union"),
      numericInput(ns("cutoff"), "Score cutoff", value = 2, min = 0, step = 0.1),
      hr(),
      ee_run_button(ns("go")),
      ee_display_head(),
      checkboxInput(ns("log_axes"), "Log scale on both axes", FALSE),
      div(class = "ee-sub-head mt-3", "Protein × allele grid"),
      radioButtons(ns("fill"), "Colour by",
                   choices = c("Number of epitopes", "Density"), selected = "Number of epitopes"),
      radioButtons(ns("plot_type"), "Plot type",
                   choices = c("Heatmap", "Bar plot"), selected = "Heatmap"),
      radioButtons(ns("sort_by"), "Order proteins by",
                   choices = c("Total", "Input order"), selected = "Total"),
      checkboxInput(ns("log_fill"), "log10 colour scale", TRUE)
    ),
    help_content = HTML(
      "<p>The scatter plot puts protein length on the x-axis and epitope count on the y-axis.
       Long proteins carry more epitopes simply because they contain more peptides, so the
       interesting proteins are the ones far off the diagonal: high count for their length.
       <b>Density</b> makes that explicit &mdash; epitopes per amino acid. Selecting rows in
       the table highlights them in the plot.</p>
       <p>The grid below breaks the same counts down by allele, so you can see whether a
       protein is broadly immunogenic or only presented by a few HLA types. A log<sub>10</sub>
       colour scale is on by default because one very long protein would otherwise flatten
       the whole scale; the colourbar is labelled with real counts.</p>
       <ul>
         <li><b>Score cutoff</b> &mdash; defaults to 2 %rank for class I, 10 for class II.</li>
         <li><b>Plot type</b> &mdash; the heatmap suits many proteins, the bar plot a handful.</li>
       </ul>
       <p class='ee-note'><b>Changed in v2:</b> protein length is the true amino-acid length
       from the FASTA. v1 used <code>nchar(sequence) - 8</code> for every predictor, so
       lengths were 8 short and every density value was slightly inflated (and, for class II
       15-mers, wrong in the other direction as a window count too).</p>"),

    ee_stat_strip_output(ns("stats")),
    bslib::layout_columns(
      col_widths = c(7, 5),
      ee_card("Epitopes vs protein length", height = "520px",
              ee_busy(plotlyOutput(ns("scatter"), height = "450px"))),
      ee_card("Per protein", download_id = ns("dl"), height = "520px",
              ee_busy(DT::dataTableOutput(ns("table"))))
    ),
    ee_card("Epitopes per protein and allele", download_id = ns("dl_grid"),
            download_label = "Download grid",
            ee_busy(ee_dynamic_plot(ns("grid_box"))))
  )
}

density_server <- function(id, ds) {
  moduleServer(id, function(input, output, session) {

    observeEvent(ds(), {
      d <- ds()
      ee_update_choices(session, "alleles", d$alleles, isolate(input$alleles),
                        length(d$alleles), "select")
      ee_apply_defaults(session, d, list(cutoff = "cutoff"))
    }, ignoreNULL = TRUE)

    p <- ee_analysis_params(ds, reactive(input$go), function(d, s) list(
      alleles = ee_pick(input$alleles, d$alleles, length(d$alleles)),
      mode    = input$mode %||% "Union",
      cutoff  = ee_val(input$cutoff, s$cutoff)
    ))

    dens <- reactive({
      d <- req(ds()); q <- p()
      ee_density_table(d, q$alleles, q$mode, q$cutoff)
    })

    grid <- reactive({
      d <- req(ds()); q <- p()
      ee_protein_allele_counts(d, q$alleles, q$cutoff)
    })

    output$stats <- renderUI({
      ee_safe({
        x <- dens()
        top <- x[which.max(Density)]
        ee_stat_strip(
          ee_stat("Proteins",              ee_num(nrow(x))),
          ee_stat("Total epitopes",        ee_num(sum(x$Epitopes))),
          ee_stat("Proteins with none",    ee_num(sum(x$Epitopes == 0L))),
          ee_stat("Median density",        sprintf("%.4f / aa", stats::median(x$Density, na.rm = TRUE))),
          ee_stat("Densest protein",       ee_short(top$ID, 22))
        )
      }, as = "ui", context = "density stats")
    })

    output$scatter <- renderPlotly({
      ee_safe({
        x <- dens()
        # Display-only: selection and log scale apply without pressing Run.
        ee_plot_density_scatter(x, selected = x$ID[input$table_rows_selected],
                                log_axes = isTRUE(input$log_axes))
      }, context = "density scatter")
    })

    tbl <- reactive({
      x <- data.table::copy(dens())
      x[, Density := round(Density, 5)]
      data.table::setnames(x, c("ID", "ProtLength", "Epitopes", "Density"),
                           c("Protein", "Length (aa)", "Epitopes", "Density (per aa)"))
      x[]
    })

    output$table <- DT::renderDataTable({
      ee_safe(ee_datatable(tbl(), page = 10), as = "table", context = "density table")
    })

    output$grid_box <- renderUI({
      n <- length(unique(grid()$ID))
      ee_render_dynamic_plot(
        session$ns("grid"),
        if (identical(input$plot_type, "Bar plot")) 480
        else ee_auto_height(min(n, EE$limits$heatmap_rows), per_row = 20, min = 320, max = 1400))
    })

    output$grid <- renderPlotly({
      ee_safe(
        # Display-only: applies without pressing Run.
        ee_plot_protein_allele(grid(),
                               fill      = input$fill %||% "Number of epitopes",
                               plot_type = input$plot_type %||% "Heatmap",
                               sort_by   = input$sort_by %||% "Total",
                               log_scale = isTRUE(input$log_fill)),
        context = "density grid")
    })

    output$dl      <- ee_download_table("epitope_density_per_protein", tbl)
    output$dl_grid <- ee_download_table("epitope_density_protein_allele", function() {
      g <- data.table::copy(grid())
      g[, Density := round(Density, 6)]
      data.table::setnames(g, c("ID", "allele", "ProtLength"),
                           c("Protein", "Allele", "Length (aa)"))
      g[]
    })
  })
}
