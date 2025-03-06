# example_distribution.R - Module for epitope distribution analysis
#
# This module provides distribution analysis for example data

#' UI for the example distribution module
#'
#' @param id Unique ID for the module
#' @return UI for the example distribution tab
example_distribution_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Distribution', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          checkboxGroupInput(
                            inputId = ns('e_allele'),
                            label = 'MHC alleles',
                            choices = colnames(example_data)[c(-1, -2, -3, -4)],
                            selected = colnames(example_data)[c(-1, -2, -3, -4)][c(1,2,3)]),
                          radioButtons(
                            inputId = ns('e_conditional'),
                            label = 'Shared epitopes',
                            choices = list('Intersection', 'Union'),
                            selected = 'Union'),
                          numericInput(
                            inputId = ns('e_xmin'),
                            label = 'Min %rank',
                            min = 0,
                            max = 10,
                            value = 0),
                          numericInput(
                            inputId = ns('e_xmax'),
                            label = 'Max %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          numericInput(
                            inputId = ns('e_step'),
                            label = 'Bin width',
                            min = 0,
                            max = 1,
                            value = 0.1),
                          radioButtons(
                            inputId = ns('e_conditional_plot'),
                            label = 'Plot type',
                            choices = list('Histogram', 'Cumulative Histogram'),
                            selected = 'Histogram'),
                          actionButton(
                            inputId = ns('e_go_2_1'),
                            label = 'Run analysis')
             ),
             mainPanel(width = 10,
                       h4("Distribution of epitopes by MHC allele"),
                       bsCollapse(bsCollapsePanel('Description (see more)', 
                                                  wellPanel(style = "overflow-x:auto; height:200px;",uiOutput(ns('e_text_distribution'))),
                                                  style = "primary")),
                       fluidRow(
                         column(5, wellPanel(withSpinner(plotlyOutput(ns('e_distribution_plot'), height = 500), type = 5))),
                         column(7, wellPanel(withSpinner(DT::dataTableOutput(ns('e_distribution_table')), type = 5))),
                         column(2, downloadButton(ns('e_downloadData2'), 'Download Table'), offset = 6),
                         column(12, wellPanel(withSpinner(plotOutput(ns('e_heatmap_plot'), height = 500), type = 5)))
                       ))
           ))
}

#' Server for the example distribution module
#'
#' @param id Unique ID for the module
#' @param e_data Reactive expression containing example data
#' @return NULL
example_distribution_server <- function(id, e_data) {
  moduleServer(id, function(input, output, session) {
    
    # Text description for epitope distribution
    output$e_text_distribution <- renderUI({
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
    
    # Get parsed data to use in multiple functions
    e_parsed_data <- reactive({
      # Data is already standardized in example_data
      return(e_data())
    })
    
    # Score data calculation
    e_ep_data.score <- reactive({
      req(input$e_allele)
      
      allele <- input$e_allele
      conditional <- input$e_conditional
      
      ep_data.score <- e_parsed_data() %>% 
        dplyr::select(Peptide, all_of(allele)) %>% 
        unique() 
      
      # Calculate score based on selected criteria
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
    
    # Distribution table
    e_distribution_table <- eventReactive(input$e_go_2_1, {
      req(e_ep_data.score())
      
      allele <- input$e_allele
      xmin <- input$e_xmin
      xmax <- input$e_xmax
      
      ep_data.score <- e_ep_data.score()
      tmp_data <- e_parsed_data()
      
      # Filter peptides by range
      keep_peptides = tmp_data$Peptide %in% ep_data.score$Peptide[ep_data.score$plot_score > xmin & 
                                                                  ep_data.score$plot_score <= xmax]
      distribution_table = tmp_data[keep_peptides, c(1, 2, 4)]
      distribution_table$Peptide_length = nchar(distribution_table[, 1])
      rownames(distribution_table) <- NULL
      
      return(distribution_table)
    })
    
    # Render distribution table
    output$e_distribution_table <- DT::renderDataTable({
      req(e_distribution_table())
      e_distribution_table()
    }, rownames = FALSE, options = list(scrollX = "300px"))
    
    # Distribution plot
    output$e_distribution_plot <- renderPlotly({
      req(input$e_go_2_1, e_ep_data.score())
      
      allele <- input$e_allele
      xmin <- input$e_xmin
      xmax <- input$e_xmax
      step <- input$e_step
      conditional_plot <- input$e_conditional_plot
      plot_values <- e_ep_data.score()$plot_score
      
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
    
    # Heatmap data
    e_heatmap_data = eventReactive(input$e_go_2_1, {
      req(input$e_allele)
      
      allele <- input$e_allele
      xmin <- input$e_xmin
      xmax <- input$e_xmax
      step <- input$e_step
      parsed_data <- e_parsed_data()
      
      cutoffs = seq(xmin, xmax, by = step)
      cutoffs = cutoffs[cutoffs != 0]
      
      heatmap_data = data.frame(alleles = c(),
                                cutoff = c(),
                                number = c())
      
      # Calculate epitope count by allele and cutoff
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
    
    # Heatmap plot
    output$e_heatmap_plot = renderPlot({
      req(e_heatmap_data())
      
      xmin <- input$e_xmin
      xmax <- input$e_xmax
      step <- input$e_step
      
      cutoffs = seq(xmin, xmax, by = step)
      cutoffs = cutoffs[cutoffs != 0]
      
      heatmap_data = e_heatmap_data()
      
      ggplot(heatmap_data, aes(x = cutoff, y = alleles, fill = number)) + 
        geom_tile(color = 'black') +
        scale_fill_gradient(low = 'white', high = '#0C6394') +
        geom_text(aes(label = number), size = 5) +
        theme_bw(base_size = 22) + 
        scale_x_continuous(breaks = cutoffs) + 
        theme(legend.position = 'none') +
        xlab('') + ylab('')
    })
    
    # Download button
    output$e_downloadData2 <- downloadHandler(
      filename = function() {
        paste("Epitope_Distribution_example.txt")
      },
      content = function(file) {
        write.table(e_distribution_table(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
    
  })
}