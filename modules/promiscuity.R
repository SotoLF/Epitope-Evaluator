# promiscuity.R -- Tool 5: epitopes that bind many alleles at once

promiscuity_ui <- function(id) {
  ns <- NS(id)
  ee_tool_page(
    title = "Epitope Promiscuity",
    subtitle = "Epitopes predicted to bind many MHC alleles, split into strong and weak binders.",
    sidebar_content = tagList(
      selectizeInput(ns("alleles"), "MHC alleles", choices = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button"), placeholder = "Select alleles")),
      numericInput(ns("min_alleles"), "Minimum number of alleles", value = 3, min = 1, step = 1),
      hr(),
      numericInput(ns("strong"), "Strong-binding cutoff", value = 0.5, min = 0, step = 0.1),
      numericInput(ns("weak"),   "Weak-binding cutoff",   value = 2,   min = 0, step = 0.1),
      div(class = "ee-hint",
          "Strong binder: score ≤ strong cutoff. Weak binder: above the strong cutoff but ≤ the weak cutoff."),
      hr(),
      ee_run_button(ns("go"))
    ),
    help_content = HTML(
      "<p>Each row of the heatmap is an epitope, each column an MHC allele. Red marks a
       <b>strong binder</b> (score at or below the strong cutoff), orange a <b>weak binder</b>
       (above the strong cutoff but at or below the weak one), pale grey a non-binder.
       Rows are sorted by promiscuity, so the broadest binders are at the top.</p>
       <p>A peptide near the top of this plot is presented by many HLA types at once, which
       is the property that makes an epitope useful in a vaccine intended for a genetically
       diverse population rather than a single patient.</p>
       <ul>
         <li><b>Minimum number of alleles</b> &mdash; only epitopes binding at least this many
             alleles (at the weak cutoff or better) are shown.</li>
         <li><b>Cutoffs</b> &mdash; default to 0.5 / 2 %rank for class I and 2 / 10 for class II.</li>
       </ul>
       <p class='ee-note'><b>Changed in v2:</b> the allele count uses <code>≤</code> the weak
       cutoff, matching the documentation and the heatmap colours. v1's table used a strict
       <code>&lt;</code> while its heatmap used <code>≤</code>, so the two disagreed for every
       epitope sitting exactly on the cutoff.</p>"),

    ee_stat_strip_output(ns("stats")),
    ee_card("Binding heatmap", ee_busy(ee_dynamic_plot(ns("plot_box")))),
    ee_card("Promiscuous epitopes", download_id = ns("dl"),
            ee_busy(DT::dataTableOutput(ns("table"))))
  )
}

promiscuity_server <- function(id, ds) {
  moduleServer(id, function(input, output, session) {

    observeEvent(ds(), {
      d <- ds()
      ee_update_choices(session, "alleles", d$alleles, isolate(input$alleles),
                        length(d$alleles), "select")
      ee_apply_defaults(session, d, list(strong = "strong", weak = "weak"))
      # Derived from the data so the tool never opens on an empty heatmap.
      updateNumericInput(session, "min_alleles",
                         value = ee_suggest_min_alleles(d, d$alleles, ee_suggest(d)$weak),
                         min = 1, max = length(d$alleles))
    }, ignoreNULL = TRUE)

    p <- ee_analysis_params(ds, reactive(input$go), function(d, s) list(
      alleles     = ee_pick(input$alleles, d$alleles, length(d$alleles)),
      strong      = ee_val(input$strong, s$strong),
      weak        = ee_val(input$weak, s$weak),
      min_alleles = ee_val(input$min_alleles,
                           ee_suggest_min_alleles(d, d$alleles, s$weak))
    ))

    pr <- reactive({
      d <- req(ds()); q <- p()
      ee_promiscuity(d, q$alleles, q$strong, q$weak, q$min_alleles)
    })

    output$stats <- renderUI({
      ee_safe({
        x <- pr(); q <- p()
        tagList(
          ee_stat_strip(
            ee_stat("Epitopes found",  ee_num(x$n_total)),
            ee_stat("Alleles",         ee_num(length(q$alleles))),
            ee_stat("Minimum alleles", ee_num(q$min_alleles)),
            ee_stat("Most promiscuous",
                    if (nrow(x$table)) sprintf("%d alleles", max(x$table$NAlleles)) else "—"),
            ee_stat("Strong-binder cells", ee_num(sum(x$matrix == "SB")))
          ),
          if (isTRUE(x$truncated))
            ee_alert(sprintf(paste("The heatmap shows the %s most promiscuous of %s epitopes.",
                                   "Raise the minimum allele count to narrow it, or use the",
                                   "download for the complete list."),
                             ee_num(nrow(x$matrix)), ee_num(x$n_total)), "warning")
        )
      }, as = "ui", context = "promiscuity stats")
    })

    output$plot_box <- renderUI({
      n <- tryCatch(nrow(pr()$matrix), error = function(e) 20L)
      ee_render_dynamic_plot(session$ns("plot"),
                             ee_auto_height(n, per_row = 16, min = 320, max = 1400))
    })

    output$plot <- renderPlotly({
      ee_safe(ee_plot_promiscuity(pr(), score_label = ee_score_label(req(ds()))),
              context = "promiscuity plot")
    })

    tbl <- reactive({
      x <- data.table::copy(pr()$table)
      if (!nrow(x)) return(x)
      data.table::setnames(x, c("Peptide", "Pos", "End", "ID", "NAlleles"),
                           c("Peptide", "Start", "End", "Protein", "N alleles"))
      x[]
    })

    output$table <- DT::renderDataTable({
      ee_safe(ee_datatable(tbl(), page = 10), as = "table", context = "promiscuity table")
    })

    output$dl <- ee_download_table("epitope_promiscuity", tbl)
  })
}
