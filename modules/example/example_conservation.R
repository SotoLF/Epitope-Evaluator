# example_conservation.R - Module for epitope conservation analysis
#
# This module provides conservation analysis for example data

#' UI for the example conservation module
#'
#' @param id Unique ID for the module
#' @return UI for the example conservation tab
example_conservation_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Conservation', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          checkboxGroupInput(
                            inputId = ns('e_protein_7'),
                            label = 'Proteins',
                            choices =  example_data %>% dplyr::select(ID) %>% unlist() %>% unique(),
                            selected = c('sp|P0DTC2|SPIKE_SARS2', 'sp|P0DTD8|NS7B_SARS2', 'sp|P0DTD1|R1AB_SARS2')),
                          numericInput(
                            inputId = ns('e_cutoff_7'),
                            label = 'Cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          
                          checkboxGroupInput(
                            inputId = ns('e_allele_7'),
                            label = 'MHC alleles',
                            choices = colnames(example_data)[c(-1, -2, -3, -4)],
                            selected = colnames(example_data)[c(-1, -2, -3, -4)][c(1,2,3)]),
                          
                          radioButtons(
                            inputId = ns('e_conditional_7'),
                            label = 'Shared epitopes',
                            choices = c('Intersection', 'Union'),
                            selected = 'Union'),
                          radioButtons(
                            inputId = ns('e_plot_type_7'),
                            label = 'Plot Type',
                            choices = list('Venn Diagram', 'Up-Set'),
                            selected = 'Venn Diagram'),
                          
                          actionButton(
                            inputId = ns('e_go_7'),
                            label = 'Run Analysis')
                          
             ),
             mainPanel(width = 10,
                       h4("Conservation epitopes across variants"),
                       bsCollapse(bsCollapsePanel('Description (see more)',
                                                  wellPanel(style= 'overflow-x:auto; height:200px;', 
                                                            uiOutput(ns('e_text_7'))),
                                                  style = 'primary')),
                       wellPanel(plotlyOutput(ns("e_plot_7"), height = "550px"),
                                 plotOutput(ns("e_plot_7_VD"), height = "550px")),
                       fluidRow(column(12, withSpinner(DT::dataTableOutput(ns("e_table_7")), type = 5)),
                                column(2, downloadButton(ns('e_downloadData7'), 'Download Table'), offset = 6))
             )
           ))
}

#' Server for the example conservation module
#'
#' @param id Unique ID for the module
#' @param e_data Reactive expression containing example data
#' @return NULL
example_conservation_server <- function(id, e_data) {
  moduleServer(id, function(input, output, session) {
    
    # Descriptive text
    output$e_text_7 <- renderUI({
      list(p(HTML("This tool allows users to identify epitopes that are present
      in multiple proteins, which can be useful to identify conserved epitopes
      across different pathogen strains. In addition, this tool allows for
      identifying epitopes gained or lost by diverse mutations. For this
      tool, users need to select the proteins, the MHC alleles of interest,
      and the cutoff percentile rank.")))
    })
    
    # Get parsed data to use in multiple functions
    e_parsed_data <- reactive({
      return(e_data())
    })
    
    # Calculate table of conserved epitopes
    e_table_epitopes <- eventReactive(input$e_go_7, {
      req(e_parsed_data(), input$e_protein_7, input$e_allele_7)
      
      parsed <- e_parsed_data()
      allele_7 <- input$e_allele_7
      cutoff_7 <- input$e_cutoff_7
      protein_7 <- input$e_protein_7
      conditional_7 <- input$e_conditional_7
      
      # Filter by selected proteins
      parsed <- parsed %>% filter(ID %in% protein_7)
      parsed_1 <- parsed %>% select(Peptide, all_of(allele_7)) %>% unique()
      
      # Calculate score based on criteria
      if (length(allele_7) != 1) {
        if (conditional_7 == 'Intersection') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, max)
        } else if (conditional_7 == 'Union') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, min) 
        }
      } else {
        parsed_1$plot_score <- parsed_1[, allele_7]
      }
      
      # Filter candidate peptides
      candidate_peptides <- parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
      parsed_2 <- parsed %>% 
        select(Peptide, ID, all_of(allele_7)) %>% 
        filter(Peptide %in% candidate_peptides)
      
      # Prepare information table
      peptide_info <- dcast(parsed_2[, c(1, 2)], formula = Peptide ~ ID) %>% as.data.frame()
      peptide_info$Proteins <- apply(peptide_info[, -1], 1, function(x) return(paste(x[!is.na(x)], collapse = '-')))
      
      # Group by shared proteins
      peptide_info <- peptide_info %>% 
        group_by(Proteins) %>%
        summarize(sequences = paste(Peptide, collapse = '-')) %>%
        as.data.frame()
      
      # Calculate number of sequences
      peptide_info$N_sequences <- apply(peptide_info %>% select(sequences), 1, 
                                        function(x) return(length(unlist(strsplit(x, split = '-')))))
      
      return(peptide_info)
    })
    
    # Render conservation table
    output$e_table_7 <- DT::renderDataTable({
      req(e_table_epitopes())
      e_table_epitopes()
    }, rownames = FALSE, options = list(scrollX = '300px'))
    
    # Venn Diagram plot for conservation
    output$e_plot_7_VD <- renderPlot({
      req(e_parsed_data(), input$e_protein_7, input$e_allele_7, input$e_go_7)
      
      allele_7 <- input$e_allele_7
      cutoff_7 <- input$e_cutoff_7
      protein_7 <- input$e_protein_7
      conditional_7 <- input$e_conditional_7
      
      # Don't show if too many proteins
      if (length(protein_7) > 6) {
        return(NULL)
      }
      
      # Filter data by proteins and alleles
      parsed <- e_parsed_data() %>% filter(ID %in% protein_7)
      parsed_1 <- parsed %>% select(Peptide, all_of(allele_7)) %>% unique()
      
      # Calculate score based on criteria
      if (length(allele_7) != 1) {
        if (conditional_7 == 'Intersection') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, max)
        } else if (conditional_7 == 'Union') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, min) 
        }
      } else {
        parsed_1$plot_score <- parsed_1[, allele_7]
      }
      
      # Filter candidate peptides
      candidate_peptides <- parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
      parsed_2 <- parsed %>% 
        select(Peptide, ID, all_of(allele_7)) %>% 
        filter(Peptide %in% candidate_peptides)
      
      # Create and process Venn diagram
      myList <- split(parsed_2$Peptide, parsed_2$ID)
      venn <- Venn(myList)
      data <- process_data(venn)
      
      # Generate plot
      ggplot() +
        geom_sf(aes(fill = count), data = venn_region(data)) +
        geom_sf(size = 1.2, color = "black", data = venn_setedge(data), show.legend = TRUE) +
        geom_sf_text(aes(label = name), data = venn_setlabel(data), nudge_y = 0.1, size = 6.5) +
        geom_sf_label(aes(label = count), fontface = "bold", data = venn_region(data), 
                      label.size = NA, fill = NA, size = 7) +
        theme_void(base_size = 20) + 
        scale_fill_gradient(low = 'white', high = 'red') +
        theme(legend.position = "none") +
        scale_x_continuous(expand = expansion(mult = 0.5))
    })
    
    # Up-Set plot for conservation
    output$e_plot_7 <- renderPlotly({
      req(e_parsed_data(), input$e_protein_7, input$e_allele_7, input$e_go_7)
      
      allele_7 <- input$e_allele_7
      cutoff_7 <- input$e_cutoff_7
      protein_7 <- input$e_protein_7
      conditional_7 <- input$e_conditional_7
      
      # Filter data by proteins and alleles
      parsed <- e_parsed_data() %>% filter(ID %in% protein_7)
      parsed_1 <- parsed %>% select(Peptide, all_of(allele_7)) %>% unique()
      
      # Calculate score based on criteria
      if (length(allele_7) != 1) {
        if (conditional_7 == 'Intersection') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, max)
        } else if (conditional_7 == 'Union') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, min) 
        }
      } else {
        parsed_1$plot_score <- parsed_1[, allele_7]
      }
      
      # Filter candidate peptides
      candidate_peptides <- parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
      parsed_2 <- parsed %>% 
        select(Peptide, ID, all_of(allele_7)) %>% 
        filter(Peptide %in% candidate_peptides)
      
      # Prepare data for Up-Set plot
      new_data2 <- parsed_2 %>% dcast(Peptide ~ ID)
      new_data2[, -1] <- apply(new_data2[, -1], 2, function(x) return(as.integer(ifelse(is.na(x), 0, 1)))) %>% 
        as.data.frame()
      
      # Normalize column names
      colnames(new_data2) <- gsub("\\|", "_", colnames(new_data2))
      
      # List of proteins
      protein <- colnames(new_data2[, -1])
      
      # Calculate intersections
      i <- 1
      x <- vector("list", 2^length(protein) - 1)
      y <- vector("list", 2^length(protein) - 1)
      
      for (j in 1:length(protein)) {
        combinations <- combn(protein, j)
        for (k in 1:ncol(combinations)) {
          x[[i]] <- sum(colSums(get_intersect_members(new_data2, combinations[, k]))) / j
          y[[i]] <- paste(combinations[, k], collapse = "-")
          i <- i + 1
        }
      }
      
      z <- do.call(rbind, Map(data.frame, names = y, values = x))
      z <- z[z$values != 0, ]
      
      # Prepare data for plot
      nintersections <- nrow(z)
      nalleles <- length(protein)
      
      set_size <- data.frame(alleles = protein, size = colSums(new_data2[, 2:ncol(new_data2)]))
      set_size <- set_size[order(set_size[, 'size']), ]
      
      ordered_alleles <- as.character(set_size$alleles)
      dict <- vector("list")
      i <- length(ordered_alleles)
      
      for (allele in ordered_alleles) {
        dict[[allele]] <- i
        i <- i - 1
      }
      
      # Set size chart
      set_size_chart <- plot_ly(
        x = set_size$size,
        y = set_size$alleles,
        type = "bar",
        orientation = "h",
        marker = list(color = "black")) %>% 
        layout(bargap = 0.4,
               xaxis = list(tickfont = list(size = 15)),
               yaxis = list(tickfont = list(size = 15),
                            categoryarray = rev(set_size$alleles),
                            categoryorder = "array"))
      
      zz <- z[order(z[, "values"], decreasing = TRUE), ]
      
      # Intersection size chart
      intersect_size_chart <- plot_ly(showlegend = FALSE) %>% 
        add_trace(
          x = 1:nintersections,
          y = zz$values,
          type = "bar",
          marker = list(color = "black", hoverinfo = "none")
        ) %>% 
        layout(xaxis = list(tickfont = list(size = 15)),
               yaxis = list(tickfont = list(size = 15))) %>% 
        add_trace(
          type = "scatter",
          mode = "text",
          x = 1:nintersections,
          y = zz$values + max(zz$values) * 0.05,
          text = zz$values,
          textfont = list(color = "black", size = 15)
        )
      
      # Prepare data for grid
      names_intersection <- as.character(zz$names)
      n_grids_rows <- 0
      
      for (name in names_intersection) {
        name_splitted <- unlist(strsplit(name, ","))
        n_grids_rows <- n_grids_rows + length(name_splitted)
      }
      
      y_axis <- vector("numeric", n_grids_rows)
      x_axis <- vector("numeric", n_grids_rows)
      combo_axis <- vector("numeric", n_grids_rows)
      name_axis <- vector("character", n_grids_rows)
      j <- 1
      x <- 1
      
      for (name in names_intersection) {
        name_splitted <- unlist(strsplit(name, "-"))
        for (i in 1:length(name_splitted)) {
          just_allele <- name_splitted[i]
          y_value <- dict[[just_allele]]
          y_axis[j] <- y_value
          x_axis[j] <- x
          combo_axis[j] <- x
          name_axis[j] <- just_allele
          j <- j + 1
        }
        x <- x + 1
      }
      
      lines2 <- data.frame(
        combo = combo_axis,
        x = x_axis,
        y = y_axis,
        name = name_axis
      )
      
      # Grid chart
      grid <- plot_ly(
        type = "scatter",
        mode = "markers",
        marker = list(color = "lightgrey", size = 8)
      ) %>% add_trace(
        type = "scatter",
        x = rep(1:nintersections, length(protein)), 
        y = unlist(lapply(1:length(protein), function(x) rep(x - 0.5, nintersections))),
        hoverinfo = "none"
      ) %>% add_trace(
        data = group_by(lines2, combo),
        type = "scatter",
        mode = "lines+markers",
        x = lines2$x,
        y = lines2$y - 0.5,
        hoverinfo = "text",
        text = ~name,
        line = list(color = "black", width = 3),
        marker = list(color = "black", size = 10)
      ) %>% 
        layout(
          xaxis = list(
            showticklabels = FALSE,
            showgrid = FALSE,
            zeroline = FALSE
          ),
          yaxis = list(
            showticklabels = FALSE,
            showgrid = TRUE,
            range = c(0, length(protein)),
            zeroline = FALSE,
            range = 1:10
          ),
          margin = list(t = 0, b = 40)
        )
      
      # Adjust y-axis title
      intersect_size_chart <- intersect_size_chart %>% 
        layout(yaxis = list(title = "Intersections size"))
      
      # Combine subplots
      s1 <- subplot(
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        set_size_chart,
        nrows = 2,
        widths = c(0.6, 0.4)
      )
      
      s2 <- subplot(
        intersect_size_chart,
        grid,
        nrows = 2,
        shareX = TRUE
      ) %>% layout(showlegend = FALSE)
      
      # Combine everything
      subplot(s1, s2, widths = c(0.3, 0.7))
    })
    
    # Control visibility for plots
    observe({
      req(input$e_plot_type_7)
      
      if (input$e_plot_type_7 == "Venn Diagram") {
        if (length(input$e_protein_7) <= 6) {
          shinyjs::hide(session$ns("e_plot_7"))
          shinyjs::show(session$ns("e_plot_7_VD"))
        } else {
          shinyjs::hide(session$ns("e_plot_7_VD"))
          shinyjs::show(session$ns("e_plot_7"))
        }
      } else {
        shinyjs::hide(session$ns("e_plot_7_VD"))
        shinyjs::show(session$ns("e_plot_7"))
      }
    })
    
    # Download button
    output$e_downloadData7 <- downloadHandler(
      filename = function() {
        paste("Conservation_epitopes_example.txt")
      },
      content = function(file) {
        req(e_table_epitopes())
        write.table(e_table_epitopes(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
  })
}