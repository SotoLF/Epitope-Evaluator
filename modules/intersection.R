# intersection.R -- Tool 2: epitopes shared between MHC allele combinations
#
# Intersection (over alleles) and Conservation (over proteins) are the same
# question with a different definition of "set", so both modules call the same
# ee_membership() / ee_combinations() pair and the same UpSet and Venn
# renderers. v1 implemented the set logic twice, and both copies used the
# O(2^k) dplyr::filter_() loop that no longer runs on current dplyr.

intersection_ui <- function(id) {
  ns <- NS(id)
  ee_tool_page(
    title = "Epitope Intersection",
    subtitle = "Which epitopes are shared between combinations of MHC alleles, and which are private to one.",
    sidebar_content = tagList(
      selectizeInput(ns("alleles"), "MHC alleles", choices = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button"),
                                    placeholder = "Select two or more alleles")),
      numericInput(ns("cutoff"), "Score cutoff", value = 2, min = 0, step = 0.1),
      div(class = "ee-hint", "A peptide joins an allele's set when its score is at or below this cutoff."),
      hr(),
      ee_run_button(ns("go")),
      ee_display_head(),
      radioButtons(ns("plot_type"), "Plot type",
                   choices = c("UpSet", "Venn diagram"), selected = "UpSet"),
      div(class = "ee-hint",
          sprintf("A Venn diagram is drawn for 2–%d alleles; UpSet handles any number.",
                  EE$limits$venn_max_sets)),
      numericInput(ns("max_bars"), "Max intersections shown",
                   value = EE$limits$upset_max_bars, min = 5, max = 200, step = 5)
    ),
    help_content = HTML(
      "<p>Each selected allele defines a set: the unique peptides predicted to bind it at or
       below the cutoff. The <b>UpSet</b> plot shows every combination that actually occurs —
       the bars on top are intersection sizes, the dot matrix says which alleles each bar
       refers to, and the grey bars on the left are the total size of each allele's set.
       A single dot means an epitope private to that allele; connected dots mean an epitope
       shared by exactly those alleles and no others.</p>
       <p>Promiscuous epitopes (large, highly-connected bars) are the attractive vaccine
       candidates: one peptide covering many HLA types covers more of a population.</p>
       <ul>
         <li><b>MHC alleles</b> &mdash; two or more. UpSet stays readable well past the point
             where a Venn diagram does not.</li>
         <li><b>Score cutoff</b> &mdash; defaults to 2 %rank for class I, 10 for class II.</li>
         <li><b>Max intersections shown</b> &mdash; only affects the figure; the table and the
             download always contain every combination.</li>
       </ul>"),

    ee_stat_strip_output(ns("stats")),
    ee_card("Set overlap", ee_busy(plotlyOutput(ns("plot"), height = "560px"))),
    ee_card("Epitopes per combination", download_id = ns("dl"),
            ee_busy(DT::dataTableOutput(ns("table"))))
  )
}

intersection_server <- function(id, ds) {
  moduleServer(id, function(input, output, session) {

    observeEvent(ds(), {
      d <- ds()
      ee_update_choices(session, "alleles", d$alleles, isolate(input$alleles),
                        min(4L, length(d$alleles)), "select")
      ee_apply_defaults(session, d, list(cutoff = "cutoff"))
    }, ignoreNULL = TRUE)

    p <- ee_analysis_params(ds, reactive(input$go), function(d, s) list(
      alleles = ee_pick(input$alleles, d$alleles, min(4L, length(d$alleles))),
      cutoff  = ee_val(input$cutoff, s$cutoff)
    ))

    mem <- reactive({
      d <- req(ds()); q <- p()
      ee_membership(d, q$alleles, q$cutoff, by = "allele")
    })
    cmb <- reactive(ee_combinations(mem()))

    output$stats <- renderUI({
      ee_safe({
        m <- mem(); c2 <- cmb(); q <- p()
        shared <- sum(c2$count[c2$n_sets > 1L])
        ee_stat_strip(
          ee_stat("Alleles compared",      ee_num(ncol(m$m))),
          ee_stat("Epitopes (any allele)", ee_num(nrow(m$m))),
          ee_stat("Shared by 2+ alleles",  ee_num(shared)),
          ee_stat("Distinct combinations", ee_num(nrow(c2))),
          ee_stat("Cutoff",                format(q$cutoff))
        )
      }, as = "ui", context = "intersection stats")
    })

    output$plot <- renderPlotly({
      ee_safe({
        m <- mem(); c2 <- cmb()
        # Display-only: applies without pressing Run.
        if (identical(input$plot_type, "Venn diagram")) ee_plot_venn(m, c2)
        else ee_plot_upset(c2, ee_set_sizes(m),
                           max_bars = max(5L, as.integer(ee_val(input$max_bars,
                                                                EE$limits$upset_max_bars))),
                           set_label = "Allele")
      }, context = "intersection plot")
    })

    tbl <- reactive({
      c2 <- cmb()
      out <- data.table::copy(c2)[, .(n_sets, sets, count, peptides)]
      out[, peptides := ee_peptide_preview(peptides)]
      data.table::setnames(out, c("N alleles", "Alleles", "N epitopes", "Epitopes"))
    })

    output$table <- DT::renderDataTable({
      ee_safe(ee_datatable(tbl(), page = 10), as = "table", context = "intersection table")
    })

    # The download keeps the complete peptide list, not the shortened preview.
    output$dl <- ee_download_table("epitope_intersection", function() {
      out <- data.table::copy(cmb())[, .(n_sets, sets, count, peptides)]
      data.table::setnames(out, c("N alleles", "Alleles", "N epitopes", "Epitopes"))
    })
  })
}
