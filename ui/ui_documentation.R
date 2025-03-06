# ui_documentation.R - UI para la pestaña Documentation
#
# Este archivo contiene la definición de la pestaña "Documentation" que proporciona
# información sobre los formatos de archivo soportados y los predictores.

#' UI para la pestaña Documentation
#'
#' @return UI para la pestaña Documentation
documentation_ui <- function() {
  tabPanel(title = span("Documentation", style = 'font-size:130%'), 
           sidebarLayout(
             sidebarPanel(
               width = 4, 
               h2('Github Documentation'),
               wellPanel(uiOutput('d_text_1'))
             ),
             mainPanel(
               width = 8, 
               h2('Supported predictors of T-cell epitopes'),
               h3('MHC Class-I epitopes'),
               bsCollapse(
                 bsCollapsePanel(
                   'NetMHCpan',
                   wellPanel(
                     style = 'overflow-x:auto;',
                     uiOutput('d_text_netmhcpan')
                   ),
                   wellPanel(
                     DT::dataTableOutput('d_table_netmhcpan')
                   ),
                   style = 'primary'
                 )
               ),
               
               bsCollapse(
                 bsCollapsePanel(
                   'NetMHC',
                   wellPanel(
                     style = 'overflow-x:auto;', 
                     uiOutput('d_text_netmhc')
                   ),
                   wellPanel(
                     DT::dataTableOutput('d_table_netmhc')
                   ),
                   style = 'primary'
                 )
               ),
               
               bsCollapse(
                 bsCollapsePanel(
                   'MHCFlurry',
                   wellPanel(
                     style = 'overflow-x:auto', 
                     uiOutput('d_text_mhcflurry')
                   ),
                   wellPanel(
                     DT::dataTableOutput('d_table_mhcflurry')
                   ),
                   style = 'primary'
                 )
               ),
               
               bsCollapse(
                 bsCollapsePanel(
                   'IEDB Consensus',
                   wellPanel(
                     style = 'overflow-x:auto;', 
                     uiOutput('d_text_iedb')
                   ),
                   wellPanel(
                     DT::dataTableOutput('d_table_iedb')
                   ),
                   style = 'primary'
                 )
               ),
               
               h3('MHC Class-II epitopes'),
               
               bsCollapse(
                 bsCollapsePanel(
                   'NetMHCIIpan',
                   wellPanel(
                     style = 'overflow-x:auto;', 
                     uiOutput('d_text_netmhciipan')
                   ),
                   wellPanel(
                     DT::dataTableOutput('d_table_netmhciipan')
                   ),
                   style = 'primary'
                 )
               )
             )
           ))
}

#' Servidor para la pestaña Documentation
#'
#' @param output Objeto output de Shiny
#' @return NULL
documentation_server <- function(output) {
  # Texto principal
  output$d_text_1 <- renderUI({
    github_code <- htmltools::a("GitHub", href = "https://github.com/SotoLF/Epitope-Evaluator")
    list(p(HTML(paste("<b>The documentation can be found on</b>", github_code))))
  })
  
  # NetMHCpan
  output$d_text_netmhcpan <- renderUI({
    predictor_url <- htmltools::a("NetMHCpan", href = "https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1")
    list(p(HTML(paste("<ul> 
                      <li>Paper: NetMHCpan-4.1 and NetMHCIIpan-4.0: Improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data</li>
                      <li>Version: 4.1</li>
                      <li>URL: ", predictor_url, "</li>
                      </ul>"))))
  })
  
  output$d_table_netmhcpan <- DT::renderDataTable({
    return(data.frame(read.table('data/example_NetMHCPAN.xls', 
                                 sep = "\t")))
  }, options = list(scrollX = "300px"))
  
  # NetMHC
  output$d_text_netmhc <- renderUI({
    predictor_url <- htmltools::a("NetMHC", href = "https://services.healthtech.dtu.dk/service.php?NetMHC-4.0")
    list(p(HTML(paste("<ul> 
                      <li>Paper: Gapped sequence alignment using artificial neural networks: application to the MHC class I system</li>
                      <li>Version: 4.0</li>
                      <li>URL: ", predictor_url, "</li>
                      </ul>"))))
  })
  
  output$d_table_netmhc <- DT::renderDataTable({
    return(data.frame(read.table('data/example_NetMHC.xls', 
                                 sep = "\t", row.names = NULL)))
  }, options = list(scrollX = "300px"))
  
  # MHCFlurry
  output$d_text_mhcflurry <- renderUI({
    predictor_url <- htmltools::a("MHCFlurry", href = "https://github.com/openvax/mhcflurry")
    list(p(HTML(paste("<ul> 
                      <li>Paper: MHCflurry 2.0: Improved pan-allele prediction of MHC I-presented peptides by incorporating antigen processing</li>
                      <li>Version: 2.0</li>
                      <li>URL: ", predictor_url, "</li>
                      </ul>"))))
  })
  
  output$d_table_mhcflurry <- DT::renderDataTable({
    return(data.frame(read.table('data/example_MHCFlurry.txt', 
                                 sep = ",")))
  }, options = list(scrollX = "300px"))
  
  # IEDB Consensus
  output$d_text_iedb <- renderUI({
    predictor_url <- htmltools::a("IEDB-Consensus", href = "http://tools.iedb.org/mhci/")
    list(p(HTML(paste("<ul> 
                      <li>Paper: A consensus epitope prediction approach identifies the breadth of murine T(CD8+)-cell responses to vaccinia virus</li>
                      <li>Version: 2.24</li>
                      <li>URL: ", predictor_url, "</li>
                      </ul>"))))
  })
  
  output$d_table_iedb <- DT::renderDataTable({
    return(data.frame(read.table('data/example_IEDB_consensus.txt', 
                                 sep = "\t", row.names = NULL)))
  }, options = list(scrollX = "300px"))
  
  # NetMHCIIpan
  output$d_text_netmhciipan <- renderUI({
    predictor_url <- htmltools::a("NetMHCIIpan", href = "https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.0")
    list(p(HTML(paste("<ul> 
                      <li>Paper: NetMHCpan-4.1 and NetMHCIIpan-4.0: Improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data</li>
                      <li>Version: 4.0</li>
                      <li>URL: ", predictor_url, "</li>
                      </ul>"))))
  })
  
  output$d_table_netmhciipan <- DT::renderDataTable({
    results <- data.frame(read.table('data/example_NetMHCIIPAN.xls',
                                     sep = ",", header = FALSE))
    
    # Obtiene nombres de alelos para Clase II
    alleles_vector = results[1, ] %>% strsplit('\t') %>% unlist()
    alleles_ID = alleles_vector[grepl('D', alleles_vector)]
    
    # Analiza datos
    results = results[c(-1, -2), ] %>% strsplit('\t')
    results = as.data.frame(do.call(rbind, results))
    
    return(results)
  }, options = list(scrollX = "300px"))
}
