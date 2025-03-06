# example_intersection.R - Module for epitope intersection analysis
#
# This module provides intersection analysis for example data

#' UI for the example intersection module
#'
#' @param id Unique ID for the module
#' @return UI for the example intersection tab
example_intersection_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Intersection', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          numericInput(
                            inputId = ns('e_cutoff_6'),
                            label = 'Cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          
                          checkboxGroupInput(
                            inputId = ns('e_allele_6'),
                            label = 'MHC alleles',
                            choices = colnames(example_data)[c(-1, -2, -3, -4)],
                            selected = colnames(example_data)[c(-1, -2, -3, -4)][c(1,2,3,4)]),
                          
                          radioButtons(
                            inputId = ns('e_plot_type_6'),
                            label = 'Plot Type',
                            choices = c('Venn Diagram', 'Up-Set'),
                            selected = 'Venn Diagram'),
                          
                          actionButton(
                            inputId = ns('e_go_6'),
                            label = 'Run Analysis')
             ),
             mainPanel(width = 10,
                      h4("Set of epitopes shared between different MHC allele combinations"),
                      bsCollapse(bsCollapsePanel('Description (see more)',
                                                 wellPanel(style= 'overflow-x:auto; height:200px;', 
                                                           uiOutput(ns('e_text_Uset'))),
                                                 style = 'primary')),
                      wellPanel(plotlyOutput(ns("e_plot_6"), height = "550px"),
                                plotOutput(ns("e_plot_6_VD"), height = "550px")),
                      wellPanel(fluidRow(column(12, withSpinner(DT::dataTableOutput(ns("e_table_6")), 
                                                                type = 5)))),
                      fluidRow(column(2, downloadButton(ns("e_downloadData6"), "Download Table"), 
                                      offset = 6))
             )
           )
  )
}

#' Server for the example intersection module
#'
#' @param id Unique ID for the module
#' @param e_data Reactive expression containing example data
#' @return NULL
example_intersection_server <- function(id, e_data) {
  moduleServer(id, function(input, output, session) {
    
    # Descriptive text
    output$e_text_Uset <- renderUI({
      list(p(HTML("Different populations are represented by distinct combinations
      of MHC alleles. This tool enables us to identify the set of epitopes that
      could be used in epitope vaccines that are potentially recognized by all
      MHC alleles in the population, as well as to identify epitopes restricted
      to a particular set of alleles. This tool shows the number of epitopes
      predicted to bind to different MHC allele combinations represented as a 
      Venn Diagram if 6 or fewer alleles are selected, or an Up-Set plot if 
      more than 6 alleles are selected. In addition to the downloadable plots,
      the tool provides a table containing the epitope sequences and the number
      of epitopes within each combination of alleles.")))
    })
    
    # Get parsed data to use in multiple functions
    e_parsed_data <- reactive({
      return(e_data())
    })
    
    # Calculate data for intersection
    e_ep_data6 <- eventReactive(input$e_go_6, {
      req(e_parsed_data(), input$e_allele_6)
      
      allele_6 <- input$e_allele_6
      cutoff_6 <- input$e_cutoff_6
      ep_data5 <- e_parsed_data()
      
      # Create matrix for analysis
      new_data <- data.frame(matrix(ncol = length(allele_6) + 1, nrow = nrow(ep_data5)))
      colnames(new_data) <- c("Peptide", allele_6)
      new_data$Peptide <- ep_data5$Peptide
      
      # Mark alleles that meet the cutoff criteria
      for (allele in allele_6) {
        new_data[allele] <- ifelse(ep_data5[allele] <= cutoff_6, allele, NA)
      }
      
      # Remove duplicates and process
      new_data <- new_data[!duplicated(new_data$Peptide), ]
      new_data$Alleles = apply(new_data[, -1], 1, function(x) return(paste(x[!is.na(x)], collapse = '-')))
      new_data = new_data %>% filter(Alleles != '')
      
      # Group by alleles
      new_data = new_data %>% 
        group_by(Alleles) %>% 
        summarize(Sequences = paste(Peptide, collapse = '-')) %>% 
        as.data.frame()
      
      # Calculate additional metrics
      new_data$N_epitopes = apply(new_data %>% select(Sequences), 1, 
                                  function(x) return(length(unlist(strsplit(x, split = '-')))))
      new_data$N_alleles = apply(new_data %>% select(Alleles), 1, 
                                 function(x) return(length(unlist(strsplit(x, split = '-')))))
      
      # Reorder columns
      new_data = new_data[, c(4, 1, 2, 3)]
      
      return(new_data)
    })
    
    # Render intersection table
    output$e_table_6 <- DT::renderDataTable({
      req(e_ep_data6())
      e_ep_data6() %>% arrange(desc(N_alleles))
    }, rownames = FALSE, options = list(scrollX = "300px"))
    
    # Venn Diagram plot
    output$e_plot_6_VD <- renderPlot({
      req(e_parsed_data(), input$e_allele_6, input$e_go_6)
      
      allele_6 <- input$e_allele_6
      cutoff_6 <- input$e_cutoff_6
      
      # Don't show if too many alleles
      if (length(allele_6) > 6) {
        return(NULL)
      }
      
      ep_data5 <- e_parsed_data()
      
      # Prepare data for Venn Diagram
      ep_data6 = ep_data5 %>% select(1, all_of(allele_6))
      ep_data6 = melt(ep_data6, id.vars = 'Peptide') %>% as.data.frame()
      ep_data6 = ep_data6 %>% filter(value <= cutoff_6)
      
      # Create Venn Diagram
      myList = split(ep_data6$Peptide, ep_data6$variable)
      venn <- Venn(myList)
      data <- process_data(venn)
      
      # Generate plot
      ggplot() +
        geom_sf(aes(fill = count), data = venn_region(data)) +
        geom_sf(size = 1.2, color = "black", data = venn_setedge(data), show.legend = TRUE) +
        geom_sf_text(aes(label = name), data = venn_setlabel(data), size = 6.5) +
        geom_sf_label(aes(label = count), fontface = "bold", data = venn_region(data), 
                      label.size = NA, fill = NA, size = 7) +
        theme_void(base_size = 20) + 
        scale_fill_gradient(low = 'white', high = 'red') +
        theme(legend.position = "none") +
        scale_x_continuous(expand = expansion(mult = 0.7))
    })
    
    # Up-Set plot
    output$e_plot_6 <- renderPlotly({
      req(e_parsed_data(), input$e_allele_6, input$e_go_6)
      
      allele_6 <- input$e_allele_6
      cutoff_6 <- input$e_cutoff_6
      
      ep_data5 <- e_parsed_data()
      
      # Prepare data for Up-Set plot
      new_data <- data.frame(matrix(ncol = length(allele_6) + 1, nrow = nrow(ep_data5)))
      colnames(new_data) <- c("Peptide", allele_6)
      new_data$Peptide <- ep_data5$Peptide
      
      for (allele in allele_6) {
        new_data[allele] <- ep_data5[allele] <= cutoff_6
        new_data[, allele] <- as.integer(new_data[, allele])
      }
      
      new_data <- new_data[order(new_data[, 'Peptide']), ]
      new_data <- new_data[!duplicated(new_data$Peptide), ]
      
      # Calculate intersections
      i <- 1
      x = vector("list", 2^length(allele_6) - 1)
      y = vector("list", 2^length(allele_6) - 1)
      
      for (j in 1:length(allele_6)) {
        combinations <- combn(allele_6, j)
        for (k in 1:ncol(combinations)) {
          x[[i]] <- sum(colSums(get_intersect_members(new_data, combinations[, k]))) / j
          y[[i]] <- paste(combinations[, k], collapse = ",")
          i = i + 1
        }
      }
      
      z <- do.call(rbind, Map(data.frame, names = y, values = x))
      z <- z[z$values != 0, ]
      
      # Prepare data for the plot
      nintersections = nrow(z)
      nalleles = length(allele_6)
      
      set_size <- data.frame(alleles = allele_6, size = colSums(new_data[, 2:ncol(new_data)]))
      set_size <- set_size[order(set_size[, 'size']), ]
      
      ordered_alleles <- as.character(set_size$alleles)
      dict <- vector("list")
      i = length(ordered_alleles)
      for (allele in ordered_alleles) {
        dict[[allele]] <- i
        i = i - 1
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
          marker = list(color = "black", hoverinfo = "none")) %>% 
        layout(xaxis = list(tickfont = list(size = 15)),
               yaxis = list(tickfont = list(size = 15))) %>% 
        add_trace(
          type = "scatter",
          mode = "text",
          x = 1:nintersections,
          y = zz$values + max(zz$values) * 0.05,
          text = zz$values,
          textfont = list(color = "black"))
      
      # Prepare data for grid
      names_intersection <- as.character(zz$names)
      n_grids_rows = 0
      for (name in names_intersection) {
        name_splitted <- unlist(strsplit(name, ","))
        n_grids_rows = n_grids_rows + length(name_splitted)
      }
      
      y_axis <- vector("numeric", n_grids_rows)
      x_axis <- vector("numeric", n_grids_rows)
      combo_axis <- vector("numeric", n_grids_rows)
      name_axis <- vector("character", n_grids_rows)
      j = 1
      x = 1
      
      for (name in names_intersection) {
        name_splitted <- unlist(strsplit(name, ","))
        for (i in 1:length(name_splitted)) {
          just_allele <- name_splitted[i]
          y_value <- dict[[just_allele]]
          y_axis[j] <- y_value
          x_axis[j] <- x
          combo_axis[j] <- x
          name_axis[j] <- just_allele
          j = j + 1
        }
        x = x + 1
      }
      
      lines2 <- data.frame(combo = combo_axis,
                           x = x_axis,
                           y = y_axis,
                           name = name_axis)
      
      # Grid chart
      grid <- plot_ly(
        type = "scatter",
        mode = "markers",
        marker = list(color = "lightgrey", size = 8)
      ) %>% add_trace(
        type = "scatter",
        x = rep(1:nintersections, length(allele_6)), 
        y = unlist(lapply(1:length(allele_6), function(x) rep(x - 0.5, nintersections))),
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
        marker = list(color = "black", size = 10)) %>% 
        layout(xaxis = list(showticklabels = FALSE,
                            showgrid = FALSE,
                            zeroline = FALSE),
               yaxis = list(showticklabels = FALSE,
                            showgrid = TRUE,
                            range = c(0, length(allele_6)),
                            zeroline = FALSE,
                            range = 1:10),
               margin = list(t = 0, b = 40))
      
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
      req(input$e_plot_type_6)
      
      if (input$e_plot_type_6 == "Venn Diagram") {
        if (length(input$e_allele_6) <= 6) {
          shinyjs::hide(session$ns("e_plot_6"))
          shinyjs::show(session$ns("e_plot_6_VD"))
        } else {
          shinyjs::hide(session$ns("e_plot_6_VD"))
          shinyjs::show(session$ns("e_plot_6"))
        }
      } else {
        shinyjs::hide(session$ns("e_plot_6_VD"))
        shinyjs::show(session$ns("e_plot_6"))
      }
    })
    
    # Download button
    output$e_downloadData6 <- downloadHandler(
      filename = function() {
        paste("Up_Set_epitopes_example.txt")
      },
      content = function(file) {
        req(e_ep_data6())
        write.table(e_ep_data6(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
    
  })
}