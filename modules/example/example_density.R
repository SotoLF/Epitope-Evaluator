# example_density.R - Module for epitope density analysis
#
# This module provides density analysis for example data

#' UI for the example density module
#'
#' @param id Unique ID for the module
#' @return UI for the example density tab
example_density_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Density', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          'Density Analysis Parameters:',
                          
                          numericInput(
                            inputId = ns('e_cutoff'),
                            label = 'Cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                            
                          actionButton(
                            inputId = ns('e_go_3_1'),
                            label = 'Run Analysis'),
                          
                          br(),
                          br(),
                          'Plot parameters:',
                          
                          radioButtons(
                            inputId = ns('e_plot_type'),
                            label = 'Plot Type',
                            choices = c('Heatmap', 'Bar plot'),
                            selected = 'Heatmap'),
                          
                          radioButtons(
                            inputId = ns('e_sort_type'),
                            label = 'Sort Type',
                            choices = c('Default', 'Descendent'),
                            selected = 'Descendent'),
                          
                          radioButtons(
                            inputId = ns('e_fill_type'),
                            label = 'Color by:',
                            choices = c('Number of Epitopes', 'Density of Epitopes'),
                            selected = 'Number of Epitopes'),
                          
                          actionButton(
                            inputId = ns('e_go_3_2'),
                            label = 'Generate Plot')
                          
             ),
             mainPanel(width = 10,
                       h4("Correlation between number of epitopes and protein length"),
                       bsCollapse(bsCollapsePanel('Description (see more)',
                                                  wellPanel(style= 'overflow-x:auto; height:200px;', 
                                                            uiOutput(ns('e_text_density'))),
                                                  style = 'primary')),
                       fluidRow(column(6, wellPanel(withSpinner(plotlyOutput(ns("e_plot_3_1"), height = 500), 
                                                                type = 5))),
                                column(6, wellPanel(withSpinner(DT::dataTableOutput(ns("e_table_3")), 
                                                                type = 5))),
                                column(2,downloadButton(ns("e_downloadData3"), "Download Table"), offset = 6)),
                       fluidRow(column(12, wellPanel(withSpinner(plotlyOutput(ns("e_plot_3_2"), height = 700), 
                                                                 type = 5))))
             )
           )
  )
}

#' Server for the example density module
#'
#' @param id Unique ID for the module
#' @param e_data Reactive expression containing example data
#' @return NULL
example_density_server <- function(id, e_data) {
  moduleServer(id, function(input, output, session) {
    
    # Descriptive text
    output$e_text_density <- renderUI({
      list(p(HTML("This tool can be used to determine the set of proteins
      containing a high number of predicted epitopes as a first step to
      finding highly immunogenic proteins. The tool displays a scatter 
      plot of protein length versus the number of epitopes predicted to bind 
      an allele or combination of MHC alleles. Hovering over each point shows 
      the name of the protein, number of epitopes, length of the protein, and 
      the epitope density. Also, selecting any protein from the table will
      highlight the respective point in the scatter plot.")))
    })
    
    # Get parsed data to use in multiple functions
    e_parsed_data <- reactive({
      return(e_data())
    })
    
    # Calculate density table
    e_d_table <- eventReactive(input$e_go_3_1, {
      req(e_parsed_data())
      
      cutoff <- input$e_cutoff
      ep_data5 <- e_parsed_data()
      
      # Calculate minimum value per row
      ep_data5$min <- apply(ep_data5[seq(5, ncol(ep_data5))], 1, min)
      ep_data5$min <- as.numeric(ep_data5$min)
      
      # Group by protein and calculate density
      ep_data5 <- ep_data5 %>% 
        dplyr::filter(min <= cutoff) %>% 
        dplyr::group_by(ID, Length) %>% 
        dplyr::summarize(count = n(), .groups = 'drop')
      
      ep_data5$Density <- round(ep_data5$count / ep_data5$Length, 2)
      
      return(ep_data5)
    })
    
    # Render density table
    output$e_table_3 <- DT::renderDataTable({
      req(e_d_table())
      DT::datatable(e_d_table(), rownames = FALSE, options = list(scrollX = "300px"))
    })
    
    # Correlation plot
    output$e_plot_3_1 <- renderPlotly({
      req(e_d_table())
      
      ep_data5 <- e_d_table()
      
      # Get selected rows
      s <- input$e_table_3_rows_selected
      
      # Define plot limits
      y_min_value <- min(ep_data5$count)
      y_max_value <- max(ep_data5$count)
      x_min_value <- min(ep_data5$Length)
      x_max_value <- max(ep_data5$Length)
      
      # Get selected protein IDs
      ids <- as.data.frame(ep_data5)[s, ]
      
      # Create plot
      plot_ly(type = "scatter") %>% 
        # Highlighted points (selected)
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
        # Normal points (not selected)
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
        # Layout configuration
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
    
    # Specific plot (heatmap or barplot)
    output$e_plot_3_2 <- renderPlotly({
      req(input$e_go_3_2, e_parsed_data())
      
      cutoff <- input$e_cutoff
      table <- e_parsed_data()
      fill_type <- input$e_fill_type
      sort_type <- input$e_sort_type
      plot_type <- input$e_plot_type
      
      # Transform data
      table_melt <- melt(table, id.vars = c("Pos", "Length", "ID", "Peptide"))
      protein_length <- table_melt %>% select(ID, Length) %>% unique()
      
      # Get IDs for axes
      y_ids <- table$ID %>% unique() %>% unlist() %>% as.character()
      x_ids <- colnames(table)[5:ncol(table)]
      
      # Calculate summary by allele and protein
      table_summarized <- table_melt %>% 
        group_by(ID, variable) %>% 
        filter(value <= cutoff) %>% 
        summarise(count = n(), .groups = 'drop') %>% 
        as.data.frame()
      
      # Assign length and calculate density
      index = match(table_summarized$ID, protein_length$ID)
      
      # Adjust fill type
      if (fill_type == "Number of Epitopes") {
        table_summarized$Length <- protein_length$Length[index] 
        table_summarized$Density <- table_summarized$count / table_summarized$Length
        table_summarized$type_fill <- log10(table_summarized$count)
      } else {
        table_summarized$Length <- protein_length$Length[index] 
        table_summarized$Density <- table_summarized$count / table_summarized$Length
        table_summarized$type_fill <- log10(table_summarized$count / table_summarized$Length)
      }
      
      # Sort if selected
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
      
      # Create plot based on selected type
      if (plot_type == "Bar plot") {
        if (fill_type == "Number of Epitopes") {
          # Bar plot by number of epitopes
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
          # Bar plot by density
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
      
      # Convert to plotly
      ggplotly(p, tooltip = c("text"))
    })
    
    # Download button
    output$e_downloadData3 <- downloadHandler(
      filename = function() {
        paste("Density_Protein_example.txt")
      },
      content = function(file) {
        req(e_d_table())
        write.table(e_d_table(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
    
  })
}