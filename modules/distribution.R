# distribution.R - Módulo para análisis de distribución de epitopos
#
# Este módulo maneja el análisis de la distribución de epitopos en función
# de los rangos percentiles (% rank) o puntajes de unión.

#' UI para el módulo de distribución de epitopos
#'
#' @param id ID único para el módulo
#' @return UI para la pestaña de distribución de epitopos
distribution_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Distribution', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          checkboxGroupInput(
                            inputId = ns('allele'),
                            label = 'MHC alleles',
                            choices = c(),
                            selected = c()),
                          radioButtons(
                            inputId = ns('conditional'),
                            label = 'Shared epitopes',
                            choices = list('Intersection', 'Union'),
                            selected = 'Intersection'),
                          numericInput(
                            inputId = ns('xmin'),
                            label = 'Min %rank',
                            min = 0,
                            max = 10,
                            value = 0),
                          numericInput(
                            inputId = ns('xmax'),
                            label = 'Max %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          numericInput(
                            inputId = ns('step'),
                            label = 'Bin width',
                            min = 0,
                            max = 1,
                            value = 0.1),
                          radioButtons(
                            inputId = ns('conditional_plot'),
                            label = 'Plot type',
                            choices = list('Histogram', 'Cumulative Histogram'),
                            selected = 'Histogram'),
                          actionButton(
                            inputId = ns('go_2_1'),
                            label = 'Run analysis')
             ),
             mainPanel(width = 10,
                       h4("Distribution of epitopes by MHC allele"),
                       bsCollapse(bsCollapsePanel('Description (see more)', 
                                                  wellPanel(style = "overflow-x:auto; height:200px;",
                                                           uiOutput(ns('text_distribution'))),
                                                  style = "primary")),
                       fluidRow(
                         column(5, wellPanel(withSpinner(plotlyOutput(ns('distribution_plot'), 
                                                                    height = 500), type = 5))),
                         column(7, wellPanel(withSpinner(DT::dataTableOutput(ns('distribution_table')), 
                                                       type = 5))),
                         column(2, downloadButton(ns('downloadData2'), 'Download Table'), offset = 6),
                         column(12, wellPanel(withSpinner(plotOutput(ns('heatmap_plot'), height = 500), 
                                                        type = 5)))
                       ))
           ))
}

#' Servidor para el módulo de distribución de epitopos
#'
#' @param id ID único para el módulo
#' @param parsed_data Reactive con datos analizados
#' @return NULL
distribution_server <- function(id, parsed_data) {
  moduleServer(id, function(input, output, session) {
    
    # Texto descriptivo
    output$text_distribution <- renderUI({
      list(p(HTML("Epitope Distribution shows a histogram representing the number of epitopes 
      within a percentile rank or binding score range. The histogram can be shown per MHC
      allele or for the union or intersection of different MHC alleles, which can represent
      the number of epitopes that can be recognized by a heterozygote individual or the 
      number of epitopes that are recognized by the alleles present in a population, 
      respectively. Users can select between plotting a histogram or a cumulative histogram.
      Hovering over any of the histogram bars provides the number of predicted epitopes with
      a percentile rank lower than a selected cutoff. In addition, the tool shows a heatmap
      indicating the number of epitopes predicted to bind to each allele with a percentile
      rank lower than a defined cutoff.

      <br>
      <h4> Parameters </h4>
      <ul>
      <li> <b>MHC alleles:</b> List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = First allele) </li>
      <li> <b>Shared epitopes:</b> When multiple alleles are selected, users must indicate if they are interested in epitopes that bind to all
      the selected MHC alleles (Intersection) or epitopes that bind to at least one of the MHC alleles selected (Union). (Default: Intersection)  </li>
      <li> <b>Min % rank:</b> Minimum percentile rank to filter the epitopes. Peptides with a rank lower
      than the cutoff selected will be ignored. To consider all epitopes, choose 0. If you are interested
      in only weak binder epitopes, you could set 0.5 as cutoff for MHC Class I and 2 for MHC Class II epitopes. (Default = 0)</li>
      <li> <b>Max % rank:</b> Maximum percentile rank to consider a peptide as predicted epitope. Peptides
      with a percentile rank higher than the cutoff selected will be ignored. (Default MHC Class I = 2, MHC Class II = 10).</li>      
      <li> <b>Step:</b> Width of the bins in the histogram/cumulative histogram plot and in the heatmap. (Default = 0.1)</li>
      <li> <b>Plot type:</b> The distribution of epitopes can be shown as a histogram or a cumulative histogram. (Default = histogram)</li>  
      </ul>")))
    })
    
    # Actualiza la lista de alelos cuando se cargan nuevos datos
    observe({
      req(parsed_data())
      alleles <- colnames(parsed_data())[c(-1, -2, -3, -4)]
      updateCheckboxGroupInput(session, inputId = "allele", 
                               choices = alleles, 
                               selected = alleles[1])
    })
    
    # Variables reactivas para inputs
    allele_react <- eventReactive(input$go_2_1, {input$allele})
    xmin_react <- eventReactive(input$go_2_1, {input$xmin})
    xmax_react <- eventReactive(input$go_2_1, {input$xmax})
    step_react <- eventReactive(input$go_2_1, {input$step})
    conditional_plot_react <- eventReactive(input$go_2_1, {input$conditional_plot})
    conditional_react <- eventReactive(input$go_2_1, {input$conditional})
    
    # Cálculo de datos de puntuación
    ep_data.score <- reactive({
      req(parsed_data(), allele_react())
      
      allele <- allele_react()
      conditional <- conditional_react()
      
      ep_data.score <- parsed_data() %>% 
        dplyr::select(Peptide, all_of(allele)) %>% 
        unique() 
      
      # Calcula puntuación según criterio seleccionado
      if (length(allele) != 1) {
        if (conditional == "Intersection") {
          ep_data.score$plot_score <- apply(ep_data.score[, allele, drop = FALSE], 1, max)
        } else if (conditional == "Union") {
          ep_data.score$plot_score <- apply(ep_data.score[, allele, drop = FALSE], 1, min)
        }
      } else {
        ep_data.score$plot_score <- ep_data.score[, allele]
      }
      
      return(ep_data.score)
    })
    
    # Tabla de distribución
    distribution_table <- reactive({
      req(ep_data.score())
      
      allele <- allele_react()
      xmin <- xmin_react()
      xmax <- xmax_react()
      
      ep_data.score <- ep_data.score()
      tmp_data <- parsed_data()
      
      # Filtra péptidos según rango
      keep_peptides = tmp_data$Peptide %in% ep_data.score$Peptide[ep_data.score$plot_score > xmin & 
                                                                    ep_data.score$plot_score <= xmax]
      distribution_table = tmp_data[keep_peptides, c(1, 2, 4)]
      distribution_table$Peptide_length = nchar(distribution_table[, 1])
      rownames(distribution_table) <- NULL
      
      return(distribution_table)
    })
    
    # Renderiza tabla
    output$distribution_table <- DT::renderDataTable({
      req(parsed_data(), distribution_table())
      distribution_table()
    }, rownames = FALSE, options = list(scrollX = "300px"))
    
    # Gráfico de distribución
    output$distribution_plot <- renderPlotly({
      req(ep_data.score())
      
      allele <- allele_react()
      xmin <- xmin_react()
      xmax <- xmax_react()
      step <- step_react()
      conditional_plot <- conditional_plot_react()
      plot_values <- ep_data.score()$plot_score
      
      cumulative = ifelse(conditional_plot == 'Histogram', FALSE, TRUE)
      
      plot_ly(x = plot_values,
              type = "histogram",
              opacity = 0.45,
              hoverinfo = "y",
              marker = list(color = "ligthblue", 
                            line = list(cauto = FALSE, width = 1, color = "black", cmid = 2)),
              xbins = list(start = xmin, end = xmax, size = step),
              cumulative = list(enabled = cumulative)) %>% 
        layout(xaxis = list(title = "% Rank Predictor", 
                            tickfont = list(size = 20), 
                            titlefont = list(size = 20)),
               yaxis = list(title = "Number of epitopes", 
                            tickfont = list(size = 20), 
                            titlefont = list(size = 20)))
    })
    
    # Datos para heatmap
    heatmap_data = reactive({
      req(parsed_data(), allele_react())
      
      allele <- allele_react()
      xmin <- xmin_react()
      xmax <- xmax_react()
      step <- step_react()
      parsed_data <- parsed_data()
      
      cutoffs = seq(xmin, xmax, by = step)
      cutoffs = cutoffs[cutoffs != 0]
      
      heatmap_data = data.frame(alleles = c(),
                                cutoff = c(),
                                number = c())
      
      # Calcula número de epitopos por alelo y cutoff
      for (cutoff in cutoffs) {
        tmp_data = apply(parsed_data %>% dplyr::select(all_of(allele)), 2, 
                         function(x) return(sum(x <= cutoff))) %>% as.data.frame()
        tmp_data$alleles = rownames(tmp_data)
        rownames(tmp_data) = NULL
        tmp_data$cutoff = cutoff
        tmp_data = tmp_data[, c(2, 3, 1)]
        colnames(tmp_data) = c('alleles', 'cutoff', 'number')
        heatmap_data = rbind(heatmap_data, tmp_data)
      }
      
      return(heatmap_data)
    })
    
    # Gráfico de heatmap
    output$heatmap_plot = renderPlot({
      req(heatmap_data())
      
      allele <- allele_react()
      xmin <- xmin_react()
      xmax <- xmax_react()
      step <- step_react()
      parsed_data <- parsed_data()
      
      cutoffs = seq(xmin, xmax, by = step)
      cutoffs = cutoffs[cutoffs != 0]
      
      heatmap_data = heatmap_data()
      
      ggplot(heatmap_data, aes(x = cutoff, y = alleles, fill = number)) + 
        geom_tile(color = 'black') +
        scale_fill_gradient(low = 'white', high = '#0C6394') +
        geom_text(aes(label = number), size = 5) +
        theme_bw(base_size = 22) + 
        scale_x_continuous(breaks = cutoffs) + 
        theme(legend.position = 'none') +
        xlab('') + ylab('')
    })
    
    # Botón de descarga
    output$downloadData2 <- downloadHandler(
      filename = function() {
        paste("Epitope_Distribution.txt")
      },
      content = function(file) {
        write.table(distribution_table(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
  })
}
