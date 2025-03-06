# example_data.R - Module for displaying example data
#
# This module provides information about the example data and displays it

#' UI for the example data module
#'
#' @param id Unique ID for the module
#' @return UI for the example data tab
example_data_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Example Data', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(
               width = 4, 
               h3('Input Data'),
               wellPanel(uiOutput(ns('e_text_1')))
             ),
             mainPanel(
               width = 8, 
               downloadButton(ns('e_downloadData_example'), 'Download Table'),
               wellPanel(DT::dataTableOutput(ns('e_contents_1')))
             )
           ))
}

#' Server for the example data module
#'
#' @param id Unique ID for the module
#' @param e_data Reactive expression containing example data
#' @return NULL
example_data_server <- function(id, e_data) {
  moduleServer(id, function(input, output, session) {
    
    # Prepare explanatory text
    output$e_text_1 <- renderUI({
      fasta_file <- htmltools::a("(Download file)", href = "https://github.com/SotoLF/Epitope-Evaluator/blob/main/Examples/example.fasta")
      prediction_file <- htmltools::a("(Download file)", href = "https://github.com/SotoLF/Epitope-Evaluator/blob/main/Examples/example.xls")
      list(p(HTML(paste("The example data shows the <b> MHC class I epitopes</b> from the SARS-CoV-2 proteome using the NetMHCPan4.1.
      <br>The input files (below) and the parsed table (right) can be downloaded:
      <ol> 
      <li>Fasta File:", fasta_file, "</li>
      <li>Output Prediction File:", prediction_file, "</li>
      </ol>
                        
      The parameters for the prediction file were the following:
                        <ul>
                        <li>PEPTIDE LENGTH: 9mer peptides</li>
                        <li>ALLELES: HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01 </li>
                        <li>Strong Binder: 0.5 </li>
                        <li>Weak Binder: 2 </li>
                        <li>Save prediction to XLS file: Yes </li>
                        </ul>
                        "))))
    })
    
    # Display example data
    output$e_contents_1 <- DT::renderDataTable({
      e_data()
    }, options = list(scrollX = '300px'))
    
    # Download button for example data
    output$e_downloadData_example <- downloadHandler(
      filename = function() {
        paste("Ouput_example.txt")
      },
      content = function(file) {
        write.table(e_data(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
  })
}