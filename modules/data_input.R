# data_input.R - Módulo para entrada y procesamiento de datos
#
# Este módulo maneja la carga y análisis de archivos de predicción de epitopos.

#' UI para el módulo de entrada de datos
#'
#' @param id ID único para el módulo
#' @return UI para la pestaña de entrada de datos
data_input_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Input Data', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          fileInput(
                            inputId = ns('input_file'),
                            label = 'Upload Prediction file'
                          ),
                          fileInput(
                            inputId = ns('input_fasta'),
                            label = 'Upload Fasta file',
                          ),
                          radioButtons(
                            inputId = ns('cond_pred'),
                            label = 'Predictor',
                            choices = list('NetMHC', 'NetMHCpan', 'NetMHCIIpan', 'MHCFlurry', 'IEDB Consensus', 'Other'),
                            selected = 'NetMHCpan'
                          ),
                          radioButtons(
                            inputId = ns('method'),
                            label = 'Score Type',
                            choices = list('Score', 'Rank'),
                            selected = 'Rank'
                          ),
                          actionButton(
                            inputId = ns('go_1'),
                            label = 'Run Analysis')
             ),
             mainPanel(width = 10, wellPanel(DT::dataTableOutput(ns('contents_1'))))
           ))
}

#' Servidor para el módulo de entrada de datos
#'
#' @param id ID único para el módulo
#' @param session Sesión actual
#' @return Reactive con datos analizados
data_input_server <- function(id, session) {
  moduleServer(id, function(input, output, session) {
    
    # Valores reactivos para tipo de predictor y método
    cond_pred <- eventReactive(input$go_1, {input$cond_pred})
    method <- eventReactive(input$go_1, {input$method})
    
    # Analiza los datos de entrada
    data <- eventReactive(input$go_1, {
      # Valida archivos de entrada
      validate(need(!is.null(input$input_file), message = "Need Output Prediction File."))
      validate(need(!is.null(input$input_fasta), message = "Need Fasta File."))
      validate(need(!is.na(tryCatch(read.fasta(input$input_fasta$datapath), error = function(e) return(NA))),
                    message = 'File submitted is not in Fasta format (See Documentation)'))
      
      # Analiza según tipo de predictor
      if (cond_pred() == "NetMHC") {
        validate(need(!is.na(tryCatch(Parse_NetMHC(input$input_file$datapath, input$input_fasta$datapath, method()), 
                                      error = function(e) return(NA))), 
                      message = 'Prediction File does not correspond to a *NetMHC* Output format (See Documentation)'))
        data <- Parse_NetMHC(input$input_file$datapath, input$input_fasta$datapath, method())
        
      } else if(cond_pred() == "NetMHCpan") {
        validate(need(!is.na(tryCatch(Parse_NetMHCPAN(input$input_file$datapath, input$input_fasta$datapath, method()), 
                                      error = function(e) return(NA))), 
                      message = 'Prediction File does not correspond to a *NetMHCPan* Output format (See Documentation)'))
        data <- Parse_NetMHCPAN(input$input_file$datapath, input$input_fasta$datapath, method())
        
      } else if(cond_pred() == "NetMHCIIpan") {
        validate(need(!is.na(tryCatch(Parse_NetMHCIIPAN(input$input_file$datapath, input$input_fasta$datapath, method()), 
                                      error = function(e) return(NA))), 
                      message = 'Prediction File does not correspond to a *NetMHCIIPan* Output format (See Documentation)'))
        data <- Parse_NetMHCIIPAN(input$input_file$datapath, input$input_fasta$datapath, method())
        
      } else if(cond_pred() == "MHCFlurry") {
        validate(need(!is.na(tryCatch(Parse_MHCFlurry(input$input_file$datapath, input$input_fasta$datapath, method()), 
                                      error = function(e) return(NA))), 
                      message = 'Prediction File does not correspond to a *MHCFlurry* Output format (See Documentation)'))
        data <- Parse_MHCFlurry(input$input_file$datapath, input$input_fasta$datapath, method())
        
      } else if(cond_pred() == "IEDB Consensus") {
        validate(need(!is.na(tryCatch(Parse_IEDB(input$input_file$datapath, input$input_fasta$datapath, method()), 
                                      error = function(e) return(NA))), 
                      message = 'Prediction File does not correspond to a *IEDB-Consensus* Output format (See Documentation)'))
        data <- Parse_IEDB(input$input_file$datapath, input$input_fasta$datapath, method())
        
      } else if (cond_pred() == "Other") {
        validate(need(!is.na(tryCatch(read.delim(input$input_file$datapath, sep = '\t') %>% as.data.frame(), 
                                      error = function(e) return(NA))), 
                      message = 'Prediction File does not correspond to the *specified* Output format (See Documentation)'))
        data <- read.delim(input$input_file$datapath, sep = "\t") %>% as.data.frame()
      }
      
      # Procesa columnas numéricas
      data[c(5:ncol(data))] <- sapply(data[c(5:ncol(data))], as.numeric)  
      data[c(5:ncol(data))] <- round(data[c(5:ncol(data))], 2)
      
      return(data)
    })
    
    # Muestra los datos analizados
    output$contents_1 <- DT::renderDataTable({
      data()
    }, options = list(scrollX = "300px"))
    
    # Procesamiento adicional de datos para uso en otros módulos
    parsed_data <- reactive({
      # Verifica si hay datos disponibles
      if (is.null(data())) {
        return(NULL)
      }
      
      # Estandariza nombres de columnas
      parsed_data <- data()
      colnames(parsed_data)[1:4] = c("Peptide", "Pos", "Length", "ID")
      
      # Actualiza elementos de la UI basados en los datos
      allele_example <- colnames(parsed_data)[5]
      if (startsWith(toupper(allele_example), "HLA")) {
        updateNumericInput(session = session, inputId = "xmin", min = 0, max = 2, value = 0)
        updateNumericInput(session = session, inputId = "xmax", min = 0, max = 2, value = 2)
        updateNumericInput(session = session, inputId = "cutoff", min = 0, max = 2, value = 2)
        updateNumericInput(session = session, inputId = "cutoff_4", min = 0, max = 2, value = 2)
        updateNumericInput(session = session, inputId = "sb_ct", min = 0, max = 2, value = 0.5)
        updateNumericInput(session = session, inputId = "wb_ct", min = 0, max = 10, value = 2)
        updateNumericInput(session = session, inputId = "cutoff_6", min = 0, max = 10, value = 2)
        updateNumericInput(session = session, inputId = "cutoff_7", min = 0, max = 10, value = 2)
      } else {
        updateNumericInput(session = session, inputId = "xmin", min = 0, max = 10, value = 0)
        updateNumericInput(session = session, inputId = "xmax", min = 0, max = 10, value = 10)
        updateNumericInput(session = session, inputId = "cutoff", min = 0, max = 10, value = 10)
        updateNumericInput(session = session, inputId = "cutoff_4", min = 0, max = 10, value = 10)
        updateNumericInput(session = session, inputId = "sb_ct", min = 0, max = 2, value = 2)
        updateNumericInput(session = session, inputId = "wb_ct", min = 0, max = 10, value = 10)
        updateNumericInput(session = session, inputId = "cutoff_6", min = 0, max = 10, value = 10)
        updateNumericInput(session = session, inputId = "cutoff_7", min = 0, max = 10, value = 10)
      }
      
      return(parsed_data)
    })
    
    # Devuelve datos analizados para uso en otros módulos
    return(parsed_data)
  })
}
