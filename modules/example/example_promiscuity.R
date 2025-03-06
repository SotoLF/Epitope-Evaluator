# example_promiscuity.R - Module for epitope promiscuity analysis
#
# This module provides promiscuity analysis for example data

#' UI for the example promiscuity module
#'
#' @param id Unique ID for the module
#' @return UI for the example promiscuity tab
example_promiscuity_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Promiscuity', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          'Set Promiscuity Analysis Parameters:',
                          br(),
                          numericInput(
                            inputId = ns('e_min_al'),
                            label = 'Minimum number of alleles',
                            min = 0,
                            max = 10,
                            value = 7),
                          
                          numericInput(
                            inputId = ns('e_sb_ct'),
                            label = 'Strong Binding cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 0.5),
                          
                          numericInput(
                            inputId = ns('e_wb_ct'),
                            label = 'Weak Binding cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          
                          actionButton(
                            inputId = ns('e_go_5'),
                            label = 'Run Analysis')
                          
             ),
             mainPanel(width = 10,
                       h4("Number of alleles and binding rank per epitopes"),
                       bsCollapse(bsCollapsePanel('Description (see more)',
                                                  wellPanel(style= 'overflow-x:auto; height:200px;', 
                                                            uiOutput(ns('e_text_5'))),
                                                  style = 'primary')),
                       fluidRow(column(6, wellPanel(withSpinner(plotlyOutput(ns("e_plot_5"), height = 500), 
                                                                type = 5))),
                                column(6, wellPanel(withSpinner(DT::dataTableOutput(ns("e_table_5")), 
                                                                type = 5)))),
                       fluidRow(column(2, downloadButton(ns("e_downloadData5"), "Download Table"), offset = 6))
             )
           ))
}

#' Server for the example promiscuity module
#'
#' @param id Unique ID for the module
#' @param e_data Reactive expression containing example data
#' @return NULL
example_promiscuity_server <- function(id, e_data) {
  moduleServer(id, function(input, output, session) {
    
    # Descriptive text
    output$e_text_5 <- renderUI({
      list(p(HTML("Epitope Promiscuity shows epitopes predicted to bind a minimum
      number of MHC alleles. The tool depicts a heatmap where alleles are on the x-axis
      and epitopes are on the y-axis. Users should set the strong and weak binding 
      cutoff %rank to filter the strong binder and weak binder epitopes. Strong binder epitopes (SB) are the 
      epitopes with a %rank less or equal than the strong binding cutoff %rank. 
      Weak binder epitopes (WB) are the epitopes with a %rank greater than the strong binding
      cutoff %rank but less or equal than the weak binding cutoff %rank.")))
    })
    
    # Get parsed data to use in multiple functions
    e_parsed_data <- reactive({
      return(e_data())
    })
    
    # Calculate promiscuity data
    e_parsed_data5 <- eventReactive(input$e_go_5, {
      req(e_parsed_data())
      
      ep_data5 <- e_parsed_data()
      wb_ct = input$e_wb_ct
      sb_ct = input$e_sb_ct
      min_al = input$e_min_al
      
      # Calculate number of alleles that meet the criteria
      ep_data5$prom <- rowSums(ep_data5[, seq(5, ncol(ep_data5))] < wb_ct)
      
      # Filter by minimum number of alleles
      ep_data5 <- ep_data5 %>% filter(prom >= min_al)
      
      # Sort and select relevant columns
      ep_data5 <- ep_data5[order(ep_data5[, "prom"], decreasing = TRUE), ]
      ep_data5 <- ep_data5 %>% select("Peptide", "Pos", "Length", "ID", "prom")
      
      return(ep_data5)
    })
    
    # Render promiscuity table
    output$e_table_5 <- DT::renderDataTable({
      req(e_parsed_data5())
      e_parsed_data5()
    }, rownames = FALSE, options = list(scrollX = "300px"))
    
    # Promiscuity plot (heatmap)
    output$e_plot_5 <- renderPlotly({
      req(e_parsed_data5())
      
      # Check if there are results
      if (nrow(e_parsed_data5()) == 0) return(NULL)
      
      wb_ct = input$e_wb_ct
      sb_ct = input$e_sb_ct
      min_al = input$e_min_al
      
      # Prepare data for plot
      ep_data5 <- e_parsed_data()
      ep_data5$prom <- rowSums(ep_data5[, seq(5, ncol(ep_data5))] < wb_ct)
      ep_data5 <- ep_data5 %>% filter(prom >= min_al)
      
      ids <- ep_data5 %>% select(Peptide) %>% unlist()
      ep_data5 <- ep_data5[order(ep_data5[, "prom"], decreasing = TRUE), ]
      
      # Create factor for ordering
      asd <- factor(levels = ep_data5$Peptide, labels = ep_data5$Peptide)
      
      # Prepare data for heatmap
      mdata <- e_parsed_data() %>% filter(Peptide %in% ids)
      mdata <- mdata[order(mdata[, 'Peptide']), ]
      mdata <- mdata[!duplicated(mdata$Peptide), ]
      mdata <- melt(mdata, id = c("Peptide", "Pos", "Length", "ID"))
      
      # Assign categories based on thresholds
      mdata$value2[mdata$value <= wb_ct & mdata$value > sb_ct] = "WB"
      mdata$value2[mdata$value <= sb_ct] = "SB"
      mdata$value2[mdata$value > wb_ct] = "NA"
      mdata$value2 <- as.factor(mdata$value2)
      
      # Create plot
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
      
      # Convert to plotly
      ggplotly(p, tooltip = c("cutoff", "y"))
    })
    
    # Download button
    output$e_downloadData5 <- downloadHandler(
      filename = function() {
        paste("Epitope_promiscuity_example.txt")
      },
      content = function(file) {
        req(e_parsed_data5())
        write.table(e_parsed_data5(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
    
  })
}