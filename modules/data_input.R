# data_input.R -- Upload, parse and summarise the input files
#
# One module serves both the Home tab (file upload) and the Run Example tab
# (bundled SARS-CoV-2 files), selected by the `preset` argument. v1 duplicated
# all seven modules into modules/example/ -- 2,149 lines of copy-paste with an
# "e_" prefix on every input id -- even though Shiny namespaces already make one
# module reusable.

data_input_ui <- function(id, preset = NULL) {
  ns <- NS(id)
  is_example <- !is.null(preset)

  sidebar <- if (is_example) {
    tagList(
      ee_alert(paste0(
        "Predicted MHC binding for the 17-protein SARS-CoV-2 proteome. ",
        "Pick a predictor output to explore, then use the six tools."), "info"),
      selectInput(ns("preset"), "Example dataset",
                  choices = ee_example_choices(), selected = names(EE_EXAMPLES)[1]),
      uiOutput(ns("preset_blurb")),
      radioButtons(ns("method"), "Score type", choices = c("Rank", "Score"), selected = "Rank"),
      div(class = "ee-hint",
          "\"Rank\" is the percentile rank against random natural peptides; ",
          "\"Score\" is the raw binding-affinity output of the predictor."),
      ee_run_button(ns("go"), "Load example")
    )
  } else {
    tagList(
      fileInput(ns("pred_file"), "Prediction file",
                accept = c(".txt", ".xls", ".tsv", ".csv", ".out"),
                placeholder = "No file selected"),
      fileInput(ns("fasta_file"), "FASTA file",
                accept = c(".fasta", ".fa", ".faa", ".fas", ".txt"),
                placeholder = "No file selected"),
      uiOutput(ns("detected")),
      selectInput(ns("predictor"), "Predictor",
                  choices = c("Auto-detect", "NetMHC", "NetMHCpan", "NetMHCIIpan",
                              "MHCFlurry", "IEDB Consensus", "Other"),
                  selected = "Auto-detect"),
      radioButtons(ns("method"), "Score type", choices = c("Rank", "Score"), selected = "Rank"),
      div(class = "ee-hint",
          "\"Rank\" is the percentile rank against random natural peptides; ",
          "\"Score\" is the raw binding-affinity output (nM or elution score)."),
      ee_run_button(ns("go"), "Load data"),
      div(class = "ee-hint mt-2",
          sprintf("Maximum upload: %d MB per file.",
                  round(getOption("shiny.maxRequestSize") / 1024^2)))
    )
  }

  bslib::layout_sidebar(
    fillable = FALSE, fill = FALSE,
    sidebar = bslib::sidebar(width = 320, title = "Input", open = "desktop",
                             class = "ee-sidebar", sidebar),
    div(class = "ee-tool-head",
        h4(if (is_example) "Example data" else "Input data", class = "ee-tool-title"),
        p(if (is_example)
            "Load a bundled prediction and explore every tool without uploading anything."
          else
            "Upload a predictor output and the FASTA it was run on, then load the data.",
          class = "ee-tool-sub")),
    uiOutput(ns("status")),
    uiOutput(ns("summary")),
    ee_card("Parsed data", download_id = ns("dl"), download_label = "Download parsed table",
            height = "560px",
            ee_busy(DT::dataTableOutput(ns("preview"))))
  )
}

#' @param id Module id.
#' @param preset NULL for the upload tab, or a label from EE_EXAMPLES.
#' @return A reactive returning an ee_dataset, or NULL before anything is loaded.
data_input_server <- function(id, preset = NULL) {
  moduleServer(id, function(input, output, session) {
    is_example <- !is.null(preset)
    state <- reactiveValues(ds = NULL, error = NULL, notes = NULL, elapsed = NA_real_)

    output$preset_blurb <- renderUI({
      ex <- EE_EXAMPLES[[input$preset %||% names(EE_EXAMPLES)[1]]]
      if (is.null(ex$blurb)) return(NULL)
      div(class = "ee-hint", ex$blurb)
    })

    # --- live format detection, before the user commits -------------------
    detected <- reactive({
      f <- input$pred_file
      if (is.null(f)) return(NULL)
      tryCatch(ee_detect_format(f$datapath), error = function(e) NULL)
    })

    output$detected <- renderUI({
      d <- detected()
      if (is.null(d)) return(NULL)
      if (is.na(d)) return(ee_alert("Format not recognised. Choose the predictor manually, or use \"Other\".", "warning"))
      ee_alert(sprintf("Detected format: <b>%s</b>", d), "success")
    })

    # Keep the dropdown in step with what was detected, without overriding a
    # deliberate manual choice.
    observeEvent(detected(), {
      d <- detected()
      if (!is.null(d) && !is.na(d) && identical(input$predictor, "Auto-detect")) {
        updateSelectInput(session, "predictor", selected = d)
      }
    })

    # --- load -------------------------------------------------------------
    # Factored out of the observer so the example tab can load once at startup
    # without faking a button press.
    load_data <- function() {
      state$error <- NULL; state$notes <- NULL

      # At module init the client has not sent its inputs yet, so fall back to
      # the same defaults the UI declares.
      method <- input$method %||% "Rank"

      args <- if (is_example) {
        nm <- input$preset %||% names(EE_EXAMPLES)[1]
        ex <- EE_EXAMPLES[[nm]]
        list(pred = ex$pred, fasta = ex$fasta, predictor = ex$predictor, label = nm)
      } else {
        list(pred = input$pred_file$datapath, fasta = input$fasta_file$datapath,
             predictor = if (identical(input$predictor %||% "Auto-detect", "Auto-detect")) "auto"
                         else input$predictor,
             label = if (is.null(input$pred_file)) NULL else input$pred_file$name)
      }

      if (is.null(args$pred)) {
        state$ds <- NULL
        state$error <- "Select a prediction file before loading."
        return()
      }
      if (is.null(args$fasta) && !identical(args$predictor, "Other")) {
        state$ds <- NULL
        state$error <- "Select the FASTA file that was submitted to the predictor."
        return()
      }

      withProgress(message = "Reading prediction file", value = 0.15, {
        t0 <- Sys.time()
        res <- tryCatch({
          incProgress(0.35, detail = "parsing")
          d <- ee_parse(args$pred, args$fasta, args$predictor, method)
          incProgress(0.4, detail = "indexing")
          d
        }, error = function(e) e)

        if (inherits(res, "error")) {
          state$ds <- NULL
          state$error <- conditionMessage(res)
        } else {
          res$source <- args$label
          state$ds <- res
          state$notes <- res$notes
          state$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        }
      })
      invisible(state$ds)
    }

    observeEvent(input$go, load_data(), ignoreInit = TRUE)

    # The example tab loads once on arrival so it is never empty. isolate()
    # keeps this from re-running when the dropdown changes -- that path goes
    # through the button like everywhere else.
    if (is_example) isolate(load_data())

    # --- status and summary ----------------------------------------------
    output$status <- renderUI({
      tags <- list()
      if (!is.null(state$error)) {
        tags <- c(tags, list(ee_alert(paste0(
          "<b>Could not load the data.</b><br>", htmltools::htmlEscape(state$error)), "danger")))
      }
      if (length(state$notes)) {
        tags <- c(tags, list(ee_alert(paste0(
          "<b>Loaded with notes:</b><ul><li>",
          paste(vapply(state$notes, htmltools::htmlEscape, ""), collapse = "</li><li>"),
          "</li></ul>"), "warning")))
      }
      if (!length(tags)) return(NULL)
      do.call(tagList, tags)
    })

    output$summary <- renderUI({
      d <- state$ds
      if (is.null(d)) {
        return(ee_alert(if (is_example) "Loading the example…" else
          "No data loaded yet. Choose your files and press <b>Load data</b>.", "info"))
      }
      ee_stat_strip(
        ee_stat("Peptide/protein rows", ee_num(d$n_rows)),
        ee_stat("Unique peptides",      ee_num(d$n_peptides)),
        ee_stat("Proteins",             ee_num(d$n_proteins)),
        ee_stat("MHC alleles",          ee_num(length(d$alleles))),
        ee_stat("MHC class",            d$mhc_class),
        ee_stat("Predictor",            d$predictor),
        ee_stat("Loaded in",            sprintf("%.1f s", state$elapsed))
      )
    })

    preview <- reactive({
      d <- req(state$ds)
      # Round for display only; downloads keep full precision.
      cbind(d$peptides[, .(Peptide, Pos, End, Length = PepLength, Protein = ID,
                           `Protein length` = ProtLength)],
            as.data.table(round(d$scores, 4)))
    })

    output$preview <- DT::renderDataTable({
      ee_safe(ee_datatable(preview(), page = 12), as = "table", context = "input preview")
    })

    output$dl <- ee_download_table("epitope_evaluator_parsed", preview)

    reactive(state$ds)
  })
}

# Bundled example datasets, declared once and shared by the Example tab.
#
# Two groups. "One predictor per format" is the same SARS-CoV-2 proteome run
# through each supported tool, so the formats can be compared like for like.
# "Biological applications" are the datasets behind the 2022 paper -- notably
# the Spike variants, which is the case the Conservation tool was built for.
EE_EXAMPLES <- list(
  "NetMHCpan 4.1 — MHC class I" = list(
    group = "One predictor per format",
    pred = "data/example_NetMHCPAN.xls", fasta = "data/example.fasta",
    predictor = "NetMHCpan",
    blurb = "SARS-CoV-2 proteome, 12 class I supertype alleles."),
  "NetMHCIIpan 4.0 — MHC class II" = list(
    group = "One predictor per format",
    pred = "data/example_NetMHCIIPAN.xls", fasta = "data/example.fasta",
    predictor = "NetMHCIIpan",
    blurb = "The same proteome, 8 DRB1 alleles. Note the wider class II cutoffs."),
  "NetMHC 4.0 — MHC class I" = list(
    group = "One predictor per format",
    pred = "data/example_NetMHC.xls", fasta = "data/example.fasta",
    predictor = "NetMHC",
    blurb = "Allele-specific predictor; its Score column is a binding affinity in nM."),
  "MHCflurry 2.0 — MHC class I" = list(
    group = "One predictor per format",
    pred = "data/example_MHCFlurry.txt", fasta = "data/example.fasta",
    predictor = "MHCFlurry",
    blurb = "Long format, one row per peptide and sample; not every allele is scored."),
  "IEDB consensus — MHC class I" = list(
    group = "One predictor per format",
    pred = "data/example_IEDB_consensus.txt", fasta = "data/example.fasta",
    predictor = "IEDB Consensus",
    blurb = "Long format keyed on seq_num, a 1-based index into the FASTA."),

  "SARS-CoV-2 proteome — class I (paper)" = list(
    group = "Biological applications (Soto et al. 2022)",
    pred = "Biological_Applications_Data/SARS_Class_I.txt", fasta = "Biological_Applications_Data/SARS_Cov2.fasta",
    predictor = "Other",
    blurb = "The class I dataset from the paper, as a generic table. Allele names arrive in the dotted form and are canonicalised on load."),
  "SARS-CoV-2 proteome — class II (paper)" = list(
    group = "Biological applications (Soto et al. 2022)",
    pred = "Biological_Applications_Data/SARS_Class_II.xls", fasta = "Biological_Applications_Data/SARS_Cov2.fasta",
    predictor = "NetMHCIIpan",
    blurb = "The class II dataset from the paper: a DRB1 + DRB3 panel, wider than the demo above."),
  "Spike variants — conservation (paper)" = list(
    group = "Biological applications (Soto et al. 2022)",
    pred = "Biological_Applications_Data/Spike_Class_II.xls", fasta = "Biological_Applications_Data/Spikes.fasta",
    predictor = "NetMHCIIpan",
    blurb = "Spike from Alpha, Beta, Delta, Gamma, Omicron and Wuhan. Go to the Conservation tool: 277 epitopes are shared by all six, and Omicron carries by far the most private ones.")
)

#' Example names grouped for a <select> with <optgroup>s
ee_example_choices <- function() {
  g <- vapply(EE_EXAMPLES, function(x) x$group %||% "Examples", character(1))
  split(names(EE_EXAMPLES), factor(g, levels = unique(g)))
}
