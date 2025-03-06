# density.R - Módulo para análisis de densidad de epitopos
#
# Este módulo maneja el análisis de la densidad de epitopos, identificando
# proteínas con alta densidad de epitopos potenciales.

#' UI para el módulo de densidad de epitopos
#'
#' @param id ID único para el módulo
#' @return UI para la pestaña de densidad de epitopos
density_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Density', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          'Density Analysis Parameters:',
                          
                          numericInput(
                            inputId = ns('cutoff'),
                            label = 'Cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          
                          actionButton(
                            inputId = ns('go_3_1'),
                            label = 'Run Analysis'),
                          
                          br(),
                          br(),
                          'Plot parameters:',
                          
                          radioButtons(
                            inputId = ns('plot_type'),
                            label = 'Plot Type',
                            choices = c('Heatmap', 'Bar plot'),
                            selected = 'Heatmap'),
                          
                          radioButtons(
                            inputId = ns('sort_type'),
                            label = 'Sort Type',
                            choices = c('Default', 'Descendent'),
                            selected = 'Descendent'),
                          
                          radioButtons(
                            inputId = ns('fill_type'),
                            label = 'Color by:',
                            choices = c('Number of Epitopes', 'Density of Epitopes'),
                            selected = 'Number of Epitopes'),
                          
                          actionButton(
                            inputId = ns('go_3_2'),
                            label = 'Generate Plot')
                          
             ),
             mainPanel(width = 10,
                       h4("Correlation between number of epitopes and protein length"),
                       bsCollapse(bsCollapsePanel('Description (see more)',
                                                  wellPanel(style = 'overflow-x:auto; height:200px;', 
                                                           uiOutput(ns('text_density'))),
                                                  style = 'primary')),
                       fluidRow(
                         column(6, wellPanel(withSpinner(plotlyOutput(ns("plot_3_1"), height = 500), 
                                                       type = 5))),
                         column(6, wellPanel(withSpinner(DT::dataTableOutput(ns("table_3")), 
                                                       type = 5))),
                         column(2, downloadButton(ns("downloadData3"), "Download Table"), offset = 6)),
                       fluidRow(
                         column(12, wellPanel(withSpinner(plotlyOutput(ns("plot_3_2"), height = 700), 
                                                        type = 5))))
             )
           ))
}

#' Servidor para el módulo de densidad de epitopos
#'
#' @param id ID único para el módulo
#' @param parsed_data Reactive con datos analizados
#' @return NULL
density_server <- function(id, parsed_data) {
  moduleServer(id, function(input, output, session) {
    
    # Texto descriptivo
    output$text_density <- renderUI({
      list(p(HTML("This tool can be used to determine the set of proteins
      containing a high number of predicted epitopes as a first step to
      finding highly immunogenic proteins. The tool displays a scatter 
      plot of protein length versus the number of epitopes predicted to bind 
      an allele or combination of MHC alleles. Hovering over each point shows 
      the name of the protein, number of epitopes, length of the protein, and 
      the epitope density. Also, selecting any protein from the table will
      highlight the respective point in the scatter plot. The tool also shows 
      the absolute number of epitopes within each protein predicted to bind to
      each MHC allele. This visualization can be displayed as a bar plot for a
      small number of proteins, or as a heatmap for several proteins. In both
      cases, users can modify the plots by changing the fill range (i.e., by
      number or by the density of epitopes) and by arranging the set of proteins.
      <br>
       To obtain the scatter plot, users must indicate the maximum cutoff
       %rank to consider a peptide as an epitope and click on the
       'Run analysis' button. This will produce a scatter plot and 
       a table containing the ID of the proteins, their length in amino acids,
       the number of epitopes, and the epitope density. 
      <br>
      To obtain the heatmap, in addition to the maximum cutoff%rank, users must
      select the plot-type (heatmap or bar plot) and the fill-type (by
      number or density of epitopes). When heatmap is selected, this shows
      the alleles on the x-axis and the proteins on the y-axis. The color 
      intensity indicates the log-10 of the number (or the density) of epitopes.
      When 'barplot' is selected, each protein is represented as a bar while
      alleles are on the x-axis and the number (or density) of epitopes are 
      on the y-axis. Hovering over each bar (or cell in the heatmap) will show 
      the protein ID, the allele, the number of epitopes, the length of the 
      protein, and the density of epitopes.
      <br>
      <h4> Parameters </h4>
      <ul>
       <li> <b>Cutoff % rank:</b> Maximum %rank to consider a peptide as a predicted epitope. (Default MHC Class I = 2, MHC Class II = 10) </li>      
       <li> <b>Color by:</b> Fill heatmap by number or density of epitopes. (Default = Epitopes Number)</li>
       <li> <b>Plot type:</b> Heatmap is recommended for several proteins, while a bar plot is recommended for a few proteins. (Default = Heatmap)</li>
       <li> <b>Sort type:</b> A parameter to sort proteins/alleles in a descendent way or as found in the input file. (Default = Descendent)</li>
       </ul>")))
    })
    
    # Variables reactivas para controlar estado
    reactive_3 <- reactiveValues(doTable = FALSE, doPlot = FALSE)
    
    # Maneja eventos de botones
    observeEvent(input$go_3_1, {
      reactive_3$doTable <- input$go_3_1
      reactive_3$doPlot <- FALSE
    })
    
    observeEvent(input$go_3_2, {
      reactive_3$doTable <- input$go_3_1
      reactive_3$doPlot <- input$go_3_2
    })
    
    # Variables reactivas para inputs
    cutoff_3_react <- eventReactive(input$go_3_1, {input$cutoff})
    cutoff_3_2_react <- eventReactive(input$go_3_2, {input$cutoff})
    fill_type_react <- eventReactive(input$go_3_2, {input$fill_type})
    plot_type_react <- eventReactive(input$go_3_2, {input$plot_type})
    sort_type_react <- eventReactive(input$go_3_2, {input$sort_type})
    
    # Calcula tabla de densidad
    d_table <- reactive({
      req(parsed_data(), cutoff_3_react())
      
      cutoff <- cutoff_3_react()
      ep_data5 <- parsed_data()
      
      # Calcula valor mínimo para cada fila y filtra
      ep_data5$min <- apply(ep_data5[seq(5, ncol(ep_data5))], 1, min)
      ep_data5$min <- as.numeric(ep_data5$min)
      
      # Agrupa por ID y Length, cuenta epitopos y calcula densidad
      ep_data5 <- ep_data5 %>% 
        dplyr::filter(min <= cutoff) %>% 
        dplyr::group_by(ID, Length) %>% 
        dplyr::summarize(count = n(), .groups = 'drop')
      
      ep_data5$Density <- round(ep_data5$count / ep_data5$Length, 2)
      
      return(ep_data5)
    })
    
    # Renderiza tabla
    output$table_3 <- DT::renderDataTable({
      req(d_table())
      DT::datatable(d_table(), rownames = FALSE, options = list(scrollX = "300px"))
    })
    
    # Gráfico de correlación
    output$plot_3_1 <- renderPlotly({
      req(d_table())
      
      ep_data5 <- d_table()
      
      # Obtiene filas seleccionadas
      s <- input$table_3_rows_selected
      
      # Define límites del gráfico
      y_min_value <- min(ep_data5$count)
      y_max_value <- max(ep_data5$count)
      x_min_value <- min(ep_data5$Length)
      x_max_value <- max(ep_data5$Length)
      
      # Obtiene IDs de proteínas seleccionadas
      ids <- as.data.frame(ep_data5)[s, ]
      
      # Crea gráfico con puntos destacados
      plot_ly(type = "scatter") %>% 
        # Puntos destacados (seleccionados)
        add_trace(
          data = ep_data5[ep_data5$ID %in% ids$ID, ],
          x = ~as.factor(Length), 
          y = ~as.factor(count),
          marker = list(
            size = 10,
            color = "red",
            line = list(color = "black", width = 2)
          ),
          showlegend = FALSE,
          hoverinfo = "text",
          text = paste(
            "Protein: ", ep_data5[ep_data5$ID %in% ids$ID, ]$ID,
            "<br>",
            "N° epitopes: ", ep_data5[ep_data5$ID %in% ids$ID, ]$count,
            "<br>",
            "Length: ", ep_data5[ep_data5$ID %in% ids$ID, ]$Length, "aminoacids",
            "<br>",
            "Density: ", round(ep_data5[ep_data5$ID %in% ids$ID, ]$Density, 2)
          )
        ) %>% 
        # Puntos normales (no seleccionados)
        add_trace(
          data = ep_data5[!(ep_data5$ID %in% ids$ID), ],
          x = ~Length, 
          y = ~count,
          marker = list(
            size = 10,
            color = "rgba(255, 182, 193, .1)",
            line = list(color = "black", width = 2)
          ),
          showlegend = FALSE,
          hoverinfo = "text",
          text = paste(
            "Protein: ", ep_data5[!(ep_data5$ID %in% ids$ID), ]$ID,
            "<br>",
            "N° epitopes: ", ep_data5[!(ep_data5$ID %in% ids$ID), ]$count,
            "<br>",
            "Length: ", ep_data5[!(ep_data5$ID %in% ids$ID), ]$Length, "aminoacids",
            "<br>",
            "Density: ", round(ep_data5[!(ep_data5$ID %in% ids$ID), ]$Density, 2)
          )
        ) %>%
        # Configuración de diseño
        layout(
          title = "Correlation N° Epitopes and Proteins Length",
          plot_bgcolor = '#a9c8f5',
          yaxis = list(
            title = "N of epitopes", 
            zerolinecolor = 'black',
            range = c(0, y_max_value),
            zerolinewidth = 2,
            gridcolor = 'black',
            tickfont = list(size = 20),  
            titlefont = list(size = 20)
          ),
          xaxis = list(
            title = "Protein Length",
            zerolinecolor = 'black', 
            range = c(0, x_max_value),
            zerolinewidth = 2, 
            gridcolor = 'black',
            tickfont = list(size = 20),  
            titlefont = list(size = 20)
          )
        )
    })
    
    # Gráfico específico (heatmap o barplot)
    output$plot_3_2 <- renderPlotly({
      req(reactive_3$doPlot, parsed_data())
      
      cutoff <- cutoff_3_2_react()
      table <- parsed_data()
      fill_type <- fill_type_react()
      sort_type <- sort_type_react()
      plot_type <- plot_type_react()
      
      # Transforma datos para análisis
      table_melt <- melt(table, id.vars = c("Pos", "Length", "ID", "Peptide"))
      protein_length <- table_melt %>% select(ID, Length) %>% unique()
      
      # Obtiene IDs únicos para ejes
      y_ids <- table$ID %>% unique() %>% unlist() %>% as.character()
      x_ids <- colnames(table)[5:ncol(table)]
      
      # Calcula resumen por alelo y proteína
      table_summarized <- table_melt %>% 
        group_by(ID, variable) %>% 
        filter(value <= cutoff) %>% 
        summarise(count = n(), .groups = 'drop') %>% 
        as.data.frame()
      
      # Asigna longitudes y calcula densidades
      index = match(table_summarized$ID, protein_length$ID)
      
      # Ajusta tipo de relleno (número o densidad)
      if (fill_type == "Number of Epitopes") {
        table_summarized$Length <- protein_length$Length[index] 
        table_summarized$Density <- table_summarized$count / table_summarized$Length
        table_summarized$type_fill <- log10(table_summarized$count)
      } else {
        table_summarized$Length <- protein_length$Length[index] 
        table_summarized$Density <- table_summarized$count / table_summarized$Length
        table_summarized$type_fill <- log10(table_summarized$count / table_summarized$Length)
      }
      
      # Ordena si se seleccionó
      if (sort_type != "Default") {
        y_ids <- table_summarized %>% 
          group_by(ID) %>% 
          summarise(accumulated = sum(10^type_fill)) %>% 
          arrange(desc(accumulated)) %>% 
          select(ID) %>% 
          unlist()
        
        x_ids <- table_summarized %>% 
          group_by(variable) %>% 
          summarise(accumulated = sum(10^type_fill)) %>% 
          arrange(desc(accumulated)) %>% 
          select(variable) %>% 
          unlist() %>% 
          as.character()
      }
      
      # Crea gráfico según tipo seleccionado
      if (plot_type == "Bar plot") {
        if (fill_type == "Number of Epitopes") {
          # Bar plot por número de epitopos
          p <- ggplot(table_summarized) + 
            geom_bar(
              stat = "identity",
              color = "black",
              position = "dodge",
              width = 0.85,
              aes(
                x = variable,
                y = count,
                fill = factor(ID, levels = y_ids),
                text = paste(
                  "Protein:", ID,
                  "<br>",
                  "Allele:", variable,
                  "<br>",
                  "N_epitopes: ", count,
                  "<br>",
                  "Length: ", Length,
                  "<br>",
                  "Density: ", round(Density, 2)
                )
              )
            ) +
            theme_classic(base_size = 20) +
            theme(
              axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15),
              axis.text.y = element_text(size = 15),
              legend.title = element_blank()
            ) +
            xlab("") + 
            ylab("") + 
            xlim(x_ids)
        } else {
          # Bar plot por densidad
          p <- ggplot(table_summarized) + 
            geom_bar(
              stat = "identity",
              color = "black",
              position = "dodge",
              width = 0.85,
              aes(
                x = variable,
                y = Density, 
                text = paste(
                  "Protein:", ID,
                  "<br>",
                  "Allele:", variable,
                  "<br>",
                  "N_epitopes: ", count,
                  "<br>",
                  "Length: ", Length,
                  "<br>",
                  "Density: ", round(Density, 2)
                ),
                fill = factor(ID, levels = y_ids)
              )
            ) +
            theme_classic(base_size = 20) +
            theme(
              axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15),
              axis.text.y = element_text(size = 15),
              legend.title = element_blank()
            ) +
            xlab("") + 
            ylab("") + 
            xlim(x_ids)
        }
      } else {
        # Heatmap
        min_value = min(table_summarized$type_fill) - 0.1
        max_value = round(max(table_summarized$type_fill), 2) + 0.1
        round_dec = ifelse(fill_type == 'Number of Epitopes', 0, 2)
        name_legend = ifelse(fill_type == 'Number of Epitopes', 'Number', 'Density')
        
        p <- ggplot(table_summarized) + 
          geom_tile(
            aes(
              x = variable, 
              y = ID, 
              fill = type_fill, 
              text = paste(
                "Protein:", ID,
                "<br>",
                "Allele:", variable,
                "<br>",
                "N_epitopes: ", count,
                "<br>",
                "Length: ", Length,
                "<br>",
                "Density: ", round(Density, 2)
              )
            )
          ) +
          theme_classic(base_size = 20) +
          scale_fill_gradientn(
            name = name_legend,
            colors = c("white", "#1E74A5"),
            limits = c(min_value, max_value),
            breaks = c(min_value, mean(c(min_value, max_value)), max_value),
            labels = c(
              round(10^min_value, round_dec), 
              round(10^mean(c(min_value, max_value)), round_dec), 
              round(10^max_value, round_dec)
            ),
            guide = guide_colorbar(
              barwidth = 1.5,
              barheight = 8, 
              nbin = 50, 
              ticks = FALSE, 
              frame.colour = "black"
            )
          ) +
          xlab("") + 
          ylab("") +
          ylim(y_ids) + 
          xlim(x_ids) +
          theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15),
            axis.text.y = element_text(size = 15)
          )
      }
      
      # Convierte a plotly con tooltip
      ggplotly(p, tooltip = c("text"))
    })
    
    # Botón de descarga
    output$downloadData3 <- downloadHandler(
      filename = function() {
        paste("Density_Protein.txt")
      },
      content = function(file) {
        write.table(d_table(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
  })
}
