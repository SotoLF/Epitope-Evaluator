# distribution.R -- Tool 1: distribution of epitopes across the score range

distribution_ui <- function(id) {
  ns <- NS(id)
  ee_tool_page(
    title = "Epitope Distribution",
    subtitle = "How predicted epitopes are distributed across the score range, per allele or across allele combinations.",
    sidebar_content = tagList(
      selectizeInput(ns("alleles"), "MHC alleles", choices = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button"), placeholder = "Select alleles")),
      radioButtons(ns("mode"), "Condition over alleles",
                   choices = c("Intersection", "Union"), selected = "Union"),
      div(class = "ee-hint",
          strong("Intersection"), " = binds every selected allele (worst rank governs). ",
          strong("Union"), " = binds at least one (best rank governs)."),
      hr(),
      numericInput(ns("xmin"), "Min score", value = 0, min = 0, step = 0.1),
      numericInput(ns("xmax"), "Max score", value = 2, min = 0, step = 0.1),
      numericInput(ns("step"), "Bin width", value = 0.1, min = 1e-6, step = 0.05),
      hr(),
      ee_run_button(ns("go")),
      ee_display_head(),
      radioButtons(ns("plot_type"), "Plot type",
                   choices = c("Histogram", "Cumulative histogram"), selected = "Histogram"),
      checkboxInput(ns("density"), "Overlay density curve", TRUE)
    ),
    help_content = HTML(
      "<p>The histogram counts <b>unique peptides</b> whose combined score falls in each bin.
       When several alleles are selected they are combined first: <b>Intersection</b> takes the
       worst (highest) rank, so a peptide is only counted if it binds them all;
       <b>Union</b> takes the best (lowest) rank, so binding any one allele is enough.
       This is the number of epitopes a heterozygous individual, or a population carrying
       those alleles, could present.</p>
      <p>The heatmap below shows, for each allele separately, how many epitopes fall below a
       ladder of increasing cutoffs, so you can see how sensitive your epitope set is to the
       threshold you picked.</p>
      <ul>
        <li><b>MHC alleles</b> &mdash; from the loaded file. Choose one or more.</li>
        <li><b>Min / Max score</b> &mdash; the window analysed. Defaults follow the MHC class
            (class I: 0&ndash;2 %rank, class II: 0&ndash;10 %rank).</li>
        <li><b>Bin width</b> &mdash; histogram bin size, and the cutoff spacing in the heatmap.</li>
        <li><b>Plot type</b> &mdash; per-bin counts, or the running total up to each cutoff.</li>
      </ul>"),

    ee_stat_strip_output(ns("stats")),
    bslib::layout_columns(
      col_widths = c(6, 6),
      ee_card("Score distribution", height = "480px",
              ee_busy(plotlyOutput(ns("hist"), height = "410px"))),
      ee_card("Epitopes in range", download_id = ns("dl"), height = "480px",
              ee_busy(DT::dataTableOutput(ns("table"))))
    ),
    ee_card("Epitopes per allele across cutoffs",
            ee_busy(ee_dynamic_plot(ns("heatmap_box"))))
  )
}

distribution_server <- function(id, ds) {
  moduleServer(id, function(input, output, session) {

    # Refresh selectors and class-appropriate defaults whenever a dataset loads.
    # v1 tried to do this from the data-input module's session, where the
    # namespaced ids could never resolve, so it silently did nothing.
    observeEvent(ds(), {
      d <- ds()
      ee_update_choices(session, "alleles", d$alleles, isolate(input$alleles), 1L, "select")
      ee_apply_defaults(session, d, list(xmin = "xmin", xmax = "xmax", step = "step"))
    }, ignoreNULL = TRUE)

    p <- ee_analysis_params(ds, reactive(input$go), function(d, s) list(
      alleles = ee_pick(input$alleles, d$alleles, 1L),
      mode    = input$mode %||% "Union",
      xmin    = ee_val(input$xmin, s$xmin),
      xmax    = ee_val(input$xmax, s$xmax),
      step    = ee_val(input$step, s$step)
    ))

    scores <- reactive({
      d <- req(ds()); q <- p()
      ee_combine(d, q$alleles, q$mode, unique_only = TRUE)
    })

    hist_data <- reactive({
      q <- p()
      ee_histogram(scores(), q$xmin, q$xmax, q$step)
    })

    tbl <- reactive({
      d <- req(ds()); q <- p()
      ee_distribution_table(d, q$alleles, q$mode, q$xmin, q$xmax)
    })

    output$stats <- renderUI({
      ee_safe({
        d <- req(ds()); q <- p(); h <- hist_data()
        ee_stat_strip(
          ee_stat("Epitopes in range", ee_num(sum(h$count))),
          ee_stat("Unique peptides scored", ee_num(d$n_peptides)),
          ee_stat("Alleles", paste(length(q$alleles), "selected")),
          ee_stat("Condition", q$mode),
          ee_stat("Bins", ee_num(nrow(h)))
        )
      }, as = "ui", context = "distribution stats")
    })

    output$hist <- renderPlotly({
      ee_safe(
        # Display-only: applies without pressing Run.
        ee_plot_histogram(hist_data(),
                          cumulative = identical(input$plot_type, "Cumulative histogram"),
                          xlab = ee_score_label(req(ds())),
                          show_density = isTRUE(input$density)),
        context = "distribution histogram")
    })

    output$table <- DT::renderDataTable({
      ee_safe(ee_datatable(tbl(), page = 10), as = "table", context = "distribution table")
    })

    output$heatmap <- renderPlotly({
      ee_safe({
        d <- req(ds()); q <- p()
        ee_check_range(q$xmin, q$xmax, q$step)
        # Cap the cutoff ladder so a tiny bin width cannot request 5,000 columns.
        n <- max(1L, min(ceiling((q$xmax - q$xmin) / q$step), 40L))
        cutoffs <- round(seq(q$xmin, q$xmax, length.out = n + 1L)[-1L], 6)
        ee_plot_cutoff_heatmap(ee_cutoff_grid(d, q$alleles, cutoffs),
                               xlab = paste(ee_score_label(d), "cutoff"))
      }, context = "distribution heatmap")
    })

    # Height follows the number of alleles so rows never squash together.
    output$heatmap_box <- renderUI({
      n <- max(length(p()$alleles), 1L)
      ee_render_dynamic_plot(session$ns("heatmap"),
                             ee_auto_height(n, per_row = 34, min = 260, max = 900))
    })

    output$dl <- ee_download_table("epitope_distribution", tbl)
  })
}
