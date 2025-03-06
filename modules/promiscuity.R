# promiscuity.R - Módulo para análisis de promiscuidad de epitopos
#
# Este módulo maneja el análisis de epitopos que se unen a múltiples alelos MHC,
# identificando epitopos promiscuos.

#' UI para el módulo de promiscuidad de epitopos
#'
#' @param id ID único para el módulo
#' @return UI para la pestaña de promiscuidad de epitopos
promiscuity_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Promiscuity', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          'Set Promiscuity Analysis Parameters:',
                          br(),
                          numericInput(
                            inputId = ns('min_al'),
                            label = 'Minimum number of alleles',
                            min = 0,
                            max = 10,
                            value = 6),
                          
                          numericInput(
                            inputId = ns('sb_ct'),
                            label = 'Strong Binding cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 0.5),
                          
                          numericInput(
                            inputId = ns('wb_ct'),
                            label = 'Weak Binding cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          
                          actionButton(
                            inputId = ns('go_5'),
                            label = 'Run Analysis')
                          
             ),
             mainPanel(width = 10,
                       h4("Number of alleles and binding rank per epitopes"),
                       bsCollapse(bsCollapsePanel('Description (see more)',
                                                  wellPanel(style = 'overflow-x:auto; height:200px;', 
                                                           uiOutput(ns('text_5'))),
                                                  style = 'primary')),
                       fluidRow(
                         column(6, wellPanel(withSpinner(plotlyOutput(ns("plot_5"), height = 500), 
                                                       type = 5))),
                         column(6, wellPanel(withSpinner(DT::dataTableOutput(ns("table_5")), 
                                                       type = 5)))
                       ),
                       fluidRow(
                         column(2, downloadButton(ns("downloadData5"), "Download Table"), offset = 6)
                       )
             )
           ))
}

#' Servidor para el módulo de promiscuidad de epitopos
#'
#' @param id ID único para el módulo
#' @param parsed_data Reactive con datos analizados
#' @return NULL
promiscuity_server <- function(id, parsed_data) {
  moduleServer(id, function(input, output, session) {
    
    # Texto descriptivo
    output$text_5 <- renderUI({
      list(p(HTML("Epitope Promiscuity shows epitopes predicted to bind a minimum
      number of MHC alleles. The tool depicts a heatmap where alleles are on the x-axis
      and epitopes are on the y-axis. Users should set the strong and weak binding 
      cutoff %rank to filter the strong binder and weak binder epitopes. Strong binder epitopes (SB) are the 
      epitopes with a %rank less or equal than the strong binding cutoff %rank. 
      Weak binder epitopes (WB) are the epitopes with a %rank greater than the strong binding
      cutoff %rank but less or equal than the weak binding cutoff %rank. In the heatmap, 
      SB and WB epitopes are colored in red and orange, respectively. Users also 
      need to select a minimum number of MHC alleles to filter and show only epitopes 
      that are predicted to bind to at least this number of alleles. This tool also shows
      a table containing the amino acid sequence, the position within the protein, the ID 
      of the protein, and the number of alleles to which they are predicted to bind.
      <br>
      <h4> Parameters </h4>
      <ul>
       <li> <b>Minimum number of alleles:</b> This indicates the minimum number of MHC alleles that an epitope needs to bind to be shown in the heatmap.</li>    
       <li> <b>Strong Binding Cutoff % rank:</b> Epitopes with a %rank lower or equal than this cutoff are considered as strong binder epitopes (SB). </li>      
       <li> <b>Weak Binding Cutoff %rank:</b> Epitopes with a %rank lower or equal than this cutoff but greater than Strong Binding Cutoff are considered weak binder epitopes (WB). </li>
       </ul>")))
    })
    
    # Actualiza el número máximo de alelos según los datos cargados
    observe({
      req(parsed_data())
      alleles <- colnames(parsed_data())[c(-1, -2, -3, -4)]
      updateNumericInput(session, "min_al",
                         max = length(alleles),
                         value = length(alleles) - 1)
    })
    
    # Variables reactivas para inputs
    wb_ct_react <- eventReactive(input$go_5, {input$wb_ct})
    sb_ct_react <- eventReactive(input$go_5, {input$sb_ct})
    min_al_react <- eventReactive(input$go_5, {input$min_al})
    
    # Cálculo de datos de promiscuidad
    parsed_data5 <- reactive({
      req(parsed_data(), wb_ct_react(), sb_ct_react(), min_al_react())
      
      ep_data5 <- parsed_data()
      wb_ct = wb_ct_react()
      sb_ct = sb_ct_react()
      min_al = min_al_react()
      
      # Calcula número de alelos que cumplen con el criterio
      ep_data5$prom <- rowSums(ep_data5[, seq(5, ncol(ep_data5))] < wb_ct)
      
      # Filtra epitopos que se unen al número mínimo de alelos
      ep_data5 <- ep_data5 %>% filter(prom >= min_al) 
      
      # Ordena y selecciona columnas relevantes
      ep_data5 <- ep_data5[order(ep_data5[, "prom"], decreasing = TRUE), ]
      ep_data5 <- ep_data5 %>% select("Peptide", "Pos", "Length", "ID", "prom")
      
      return(ep_data5)
    })
    
    # Renderiza tabla
    output$table_5 <- DT::renderDataTable({
      req(parsed_data5())
      parsed_data5()
    }, rownames = FALSE, options = list(scrollX = "300px"))
    
    # Gráfico de heatmap
    output$plot_5 <- renderPlotly({
      req(parsed_data5())
      
      # Comprueba si hay resultados
      if (nrow(parsed_data5()) == 0) return(NULL)
      
      wb_ct = wb_ct_react()
      sb_ct = sb_ct_react()
      min_al = min_al_react()
      
      # Obtiene datos completos y filtra
      ep_data5 <- parsed_data()
      ep_data5$prom <- rowSums(ep_data5[, seq(5, ncol(ep_data5))] < wb_ct)
      ep_data5 <- ep_data5 %>% filter(prom >= min_al) 
      
      # Prepara datos para gráfico
      ids <- ep_data5 %>% select(Peptide) %>% unlist()
      ep_data5 <- ep_data5[order(ep_data5[, "prom"], decreasing = TRUE), ]
      
      # Crea factor para ordenar péptidos
      asd <- factor(levels = ep_data5$Peptide, labels = ep_data5$Peptide)
      
      # Selecciona datos y prepara para heatmap
      mdata <- parsed_data() %>% filter(Peptide %in% ids)
      mdata <- mdata[order(mdata[, 'Peptide']), ]
      mdata <- mdata[!duplicated(mdata$Peptide), ]
      mdata <- melt(mdata, id = c("Peptide", "Pos", "Length", "ID"))
      
      # Clasifica según umbral
      mdata$value2[mdata$value <= wb_ct & mdata$value > sb_ct] = "WB"
      mdata$value2[mdata$value <= sb_ct] = "SB"
      mdata$value2[mdata$value > wb_ct] = "NA"
      mdata$value2 <- as.factor(mdata$value2)
      
      # Crea gráfico
      p <- ggplot(mdata) + 
        geom_tile(
          aes(
            x = variable, 
            y = Peptide, 
            fill = value2, 
            cutoff = value
          ), 
          width = 0.9, 
          height = 0.9, 
          color = "black", 
          size = 0.3
        ) +
        scale_fill_manual(values = c("NA" = "white", "SB" = "red", "WB" = "orange")) +
        ylim(levels(asd)) + 
        xlab("Alleles") + 
        ylab("") +
        theme_bw(base_size = 18) +
        xlab('') +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.title = element_blank()
        ) + 
        xlab("Alleles")
      
      # Convierte a plotly con tooltip
      ggplotly(p, tooltip = c("cutoff", "y"))
    })
    
    # Botón de descarga
    output$downloadData5 <- downloadHandler(
      filename = function() {
        paste("Epitope_promiscuity.txt")
      },
      content = function(file) {
        write.table(parsed_data5(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
  })
}
