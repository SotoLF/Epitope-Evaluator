# conservation.R -- Tool 6: epitopes shared between proteins / strains / variants
#
# Same set machinery as Tool 2, with proteins as the sets instead of alleles.

conservation_ui <- function(id) {
  ns <- NS(id)
  ee_tool_page(
    title = "Epitope Conservation",
    subtitle = "Which epitopes are conserved across proteins, strains or variants, and which are unique to one.",
    sidebar_content = tagList(
      selectizeInput(ns("proteins"), "Proteins", choices = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button"),
                                    placeholder = "Select two or more proteins")),
      hr(),
      selectizeInput(ns("alleles"), "MHC alleles", choices = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button"), placeholder = "Select alleles")),
      radioButtons(ns("mode"), "Condition over alleles",
                   choices = c("Intersection", "Union"), selected = "Union"),
      numericInput(ns("cutoff"), "Score cutoff", value = 2, min = 0, step = 0.1),
      hr(),
      ee_run_button(ns("go")),
      ee_display_head(),
      radioButtons(ns("plot_type"), "Plot type",
                   choices = c("UpSet", "Venn diagram"), selected = "UpSet"),
      numericInput(ns("max_bars"), "Max intersections shown",
                   value = EE$limits$upset_max_bars, min = 5, max = 200, step = 5)
    ),
    help_content = HTML(
      "<p>An epitope is first selected by the allele filter (a peptide must bind the chosen
       alleles at or below the cutoff, under the Intersection or Union condition), then
       assigned to every protein whose sequence contains it. The plot shows which
       combinations of proteins share which epitopes.</p>
       <p>Run the prediction on a multi-FASTA of the <i>same</i> protein from several strains
       or variants and this becomes a conservation analysis: epitopes in the all-proteins
       intersection are the ones a vaccine could target without losing coverage as the
       pathogen drifts. Epitopes private to one variant mark the positions where escape has
       already happened.</p>
       <ul>
         <li><b>Proteins</b> &mdash; two or more, from the loaded FASTA.</li>
         <li><b>MHC alleles / condition</b> &mdash; which epitopes qualify before the
             cross-protein comparison is made.</li>
         <li><b>Score cutoff</b> &mdash; defaults to 2 %rank for class I, 10 for class II.</li>
       </ul>"),

    ee_stat_strip_output(ns("stats")),
    ee_card("Protein overlap", ee_busy(plotlyOutput(ns("plot"), height = "560px"))),
    ee_card("Epitopes per protein combination", download_id = ns("dl"),
            ee_busy(DT::dataTableOutput(ns("table"))))
  )
}

conservation_server <- function(id, ds) {
  moduleServer(id, function(input, output, session) {

    observeEvent(ds(), {
      d <- ds()
      ee_update_choices(session, "proteins", d$proteins$ID, isolate(input$proteins),
                        min(3L, d$n_proteins), "select")
      ee_update_choices(session, "alleles", d$alleles, isolate(input$alleles),
                        min(3L, length(d$alleles)), "select")
      ee_apply_defaults(session, d, list(cutoff = "cutoff"))
    }, ignoreNULL = TRUE)

    p <- ee_analysis_params(ds, reactive(input$go), function(d, s) list(
      proteins = ee_pick(input$proteins, d$proteins$ID, min(3L, d$n_proteins)),
      alleles  = ee_pick(input$alleles, d$alleles, min(3L, length(d$alleles))),
      mode     = input$mode %||% "Union",
      cutoff   = ee_val(input$cutoff, s$cutoff)
    ))

    mem <- reactive({
      d <- req(ds()); q <- p()
      ee_membership(d, q$proteins, q$cutoff, by = "protein",
                    alleles = q$alleles, mode = q$mode)
    })
    cmb <- reactive(ee_combinations(mem()))

    output$stats <- renderUI({
      ee_safe({
        m <- mem(); c2 <- cmb(); q <- p()
        core <- sum(c2$count[c2$n_sets == ncol(m$m)])
        uniq <- sum(c2$count[c2$n_sets == 1L])
        ee_stat_strip(
          ee_stat("Proteins compared",   ee_num(ncol(m$m))),
          ee_stat("Epitopes considered", ee_num(nrow(m$m))),
          ee_stat("In every protein",    ee_num(core)),
          ee_stat("Unique to one",       ee_num(uniq)),
          ee_stat("Cutoff",              format(q$cutoff))
        )
      }, as = "ui", context = "conservation stats")
    })

    output$plot <- renderPlotly({
      ee_safe({
        m <- mem(); c2 <- cmb()
        # Display-only: applies without pressing Run.
        if (identical(input$plot_type, "Venn diagram")) ee_plot_venn(m, c2)
        else ee_plot_upset(c2, ee_set_sizes(m),
                           max_bars = max(5L, as.integer(ee_val(input$max_bars,
                                                                EE$limits$upset_max_bars))),
                           set_label = "Protein")
      }, context = "conservation plot")
    })

    tbl <- reactive({
      out <- data.table::copy(cmb())[, .(n_sets, sets, count, peptides)]
      out[, peptides := ee_peptide_preview(peptides)]
      data.table::setnames(out, c("N proteins", "Proteins", "N epitopes", "Epitopes"))
    })

    output$table <- DT::renderDataTable({
      ee_safe(ee_datatable(tbl(), page = 10), as = "table", context = "conservation table")
    })

    # The download keeps the complete peptide list, not the shortened preview.
    output$dl <- ee_download_table("epitope_conservation", function() {
      out <- data.table::copy(cmb())[, .(n_sets, sets, count, peptides)]
      data.table::setnames(out, c("N proteins", "Proteins", "N epitopes", "Epitopes"))
    })
  })
}
