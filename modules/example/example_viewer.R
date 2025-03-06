# example_viewer.R - Module for epitope viewer analysis
#
# This module provides viewer functionality for example data

#' UI for the example viewer module
#'
#' @param id Unique ID for the module
#' @return UI for the example viewer tab
example_viewer_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Viewer', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          selectInput(
                            inputId = ns('e_protein_4'),
                            label = 'Protein',
                            choices =  example_data %>% dplyr::select(ID) %>% unlist() %>% unique(),
                            selected = 'sp|P0DTC6|NS6_SARS2'),
                          
                          numericInput(
                            inputId = ns('e_cutoff_4'),
                            label = 'Cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          
                          checkboxGroupInput(
                            inputId = ns('e_allele_4'),
                            label = 'MHC alleles',
                            choices = colnames(example_data)[c(-1, -2, -3, -4)],
                            selected = colnames(example_data)[c(-1, -2, -3, -4)][c(1,2,3,4)]),
                          
                          radioButtons(
                            inputId = ns('e_conditional_4'),
                            label = 'Shared epitopes',
                            choices = c('Intersection', 'Union'),
                            selected = 'Union'),
                          
                          actionButton(
                            inputId = ns('e_go_4'),
                            label = 'Run Analysis')
                          
             ),
             mainPanel(width = 10,
                       h4("Amino acid position of the epitopes in proteins"),
                       bsCollapse(bsCollapsePanel('Description (see more)',
                                                  wellPanel(style= 'overflow-x:auto; height:200px;', 
                                                            uiOutput(ns('e_text_4'))),
                                                  style = 'primary')),
                       fluidRow(column(12, wellPanel(withSpinner(plotlyOutput(ns("e_plot_4")), type = 5)))),
                       fluidRow(column(12, wellPanel(withSpinner(DT::dataTableOutput(ns("e_table_4")), 
                                                                 type = 5))),
                                column(2,downloadButton(ns("e_downloadData4"), "Download Table"), offset = 6)))
           )
  )
}

#' Server for the example viewer module
#'
#' @param id Unique ID for the module
#' @param e_data Reactive expression containing example data
#' @return NULL
example_viewer_server <- function(id, e_data) {
  moduleServer(id, function(input, output, session) {
    
    # Descriptive text
    output$e_text_4 <- renderUI({
      list(p(HTML("This tool can be used to identify the location of epitopes
      predicted to bind single and multiple MHC alleles, facilitating the visual
      identification of regions enriched with promiscuous epitopes. This tool 
      graphically represents the position of each predicted epitope within 
      each protein. Epitope Viewer shows the protein as a light blue bar 
      while epitopes are represented on a color gradient from yellow to red 
      reflecting the number of MHC alleles they are predicted to bind to.")))
    })
    
    # Get parsed data to use in multiple functions
    e_parsed_data <- reactive({
      return(e_data())
    })
    
    # Calculate location table
    e_location_table <- eventReactive(input$e_go_4, {
      req(e_parsed_data(), input$e_allele_4, input$e_protein_4)
      
      ep_data5 <- e_parsed_data()
      allele_4 <- input$e_allele_4
      cutoff_4 <- input$e_cutoff_4
      protein_4 <- input$e_protein_4
      conditional_4 <- input$e_conditional_4
      
      # Filter by protein and alleles
      peptide_counts = e_parsed_data() %>% 
        filter(ID == protein_4) %>% 
        select(Peptide, Pos, Length, all_of(allele_4))
      
      # Calculate score based on criteria
      if (length(allele_4) != 1) {
        if (conditional_4 == "Intersection") {
          peptide_counts$score <- apply(peptide_counts[, allele_4, drop = FALSE], 1, max)
        } else if (conditional_4 == "Union") {
          peptide_counts$score <- apply(peptide_counts[, allele_4, drop = FALSE], 1, min)
        }
      } else {
        peptide_counts$score <- peptide_counts[, allele_4]
      }
      
      # Process data
      peptide_counts <- peptide_counts %>% 
        filter(score <= cutoff_4) %>% 
        dplyr::select(Peptide, Pos, Length, all_of(allele_4)) %>% 
        melt(id.vars = c("Peptide", "Pos", "Length")) %>% 
        group_by(Peptide) %>% 
        filter(value <= cutoff_4) %>% 
        mutate(
          count = sum(value <= cutoff_4),
          alleles = paste(variable, collapse = " - ")
        ) %>% 
        select(Peptide, count, alleles) %>% 
        as.data.frame() %>%
        distinct()
      
      # Check if there are results
      if (nrow(peptide_counts) == 0) return(NULL)
      
      # Get peptide details
      ep_data2 = e_parsed_data() %>%
        filter(ID == protein_4) %>% 
        select(Peptide, Pos, Length)
      
      # Combine information
      index = match(ep_data2$Peptide, peptide_counts$Peptide)
      ep_data2$Value = peptide_counts$count[index]
      ep_data2$alleles = peptide_counts$alleles[index]
      ep_data2 <- ep_data2 %>% filter(!is.na(Value) & Value > 0)
      
      # Prepare table for visualization
      length_protein <- ep_data2 %>% dplyr::select(Length) %>% unlist() %>% unique()
      local_table <- data.frame(
        type = c("Protein", rep("Epitope", nrow(ep_data2))),
        start = c(1, ep_data2$Pos + 1),
        end = c(length_protein + 1, ep_data2$Pos + 1 + nchar(ep_data2$Peptide)),
        sequence = c(protein_4, ep_data2$Peptide),
        nivel = c(1, rep(0, nrow(ep_data2))),
        id = seq(0, nrow(ep_data2)),
        alleles = c("-", ep_data2$alleles),
        count = c("NA", ep_data2$Value)
      )
      
      # Algorithm to organize visualization levels
      nivel = 0
      while (nivel %in% local_table$nivel) {
        subdata <- local_table[local_table$type == "Epitope" & local_table$nivel == nivel, ]
        changed = TRUE
        
        while (changed) {
          number_rows = nrow(subdata)
          
          # Special case: only one epitope
          if (number_rows == 1) {
            changed = FALSE
            nivel = nivel - 1
            break
          }
          
          conteo = 0
          for (i in 1:number_rows) {
            # Static row
            static_row <- subdata[i, ]
            # Static vector (positions)
            v_s <- seq(static_row$start, static_row$end)
            new_data <- subdata[-i, ]
            
            # Compare with all other epitopes
            conteo = conteo + 1
            for (row in 1:nrow(new_data)) {
              row_data <- new_data[row, ]
              v_x <- seq(row_data$start, row_data$end)
              id <- row_data$id
              
              # If there's overlap, adjust level
              if (length(intersect(v_s, v_x)) != 0) {
                local_table[local_table$id == id, ]$nivel <- row_data$nivel - 1
                subdata[subdata$id == id, ]$nivel <- row_data$nivel - 1
              }
            }
            
            # Update subdata with new level
            subdata <- subdata[subdata$nivel == nivel, ]
            final_number_rows = nrow(subdata)
            
            if (number_rows != final_number_rows) {
              break
            }
            
            if (number_rows == final_number_rows & conteo == number_rows) {
              changed = FALSE
            }
          }
        }
        nivel = nivel - 1
      }
      
      return(local_table)
    })
    
    # Render position table
    output$e_table_4 <- DT::renderDataTable({
      req(e_location_table())
      
      local_table <- e_location_table()
      if (is.null(local_table)) {
        return(data.frame(
          sequence = character(),
          start = integer(),
          end = integer(),
          alleles = character(),
          count = integer()
        ))
      }
      
      # Exclude protein row and select columns
      local_table <- local_table[-1, c(4, 2, 3, 7, 8)]
      
      DT::datatable(local_table, rownames = FALSE, options = list(scrollX = "300px"))
    })
    
    # Position plot
    output$e_plot_4 <- renderPlotly({
      req(e_location_table())
      
      protein_4 <- input$e_protein_4
      local_table <- e_location_table()
      
      if (is.null(local_table)) return(NULL)
      
      # Convert count to numeric
      local_table$count <- as.numeric(local_table$count)
      
      # Special case: only protein without epitopes
      if (nrow(local_table) == 1) {
        p <- ggplot() +
          geom_rect(
            data = local_table, 
            aes(
              xmin = start,
              xmax = end,
              ymin = nivel - 0.3,
              ymax = nivel + 0.3,
              fill = count,
              sequence = "-",
              text = paste("Class:", type, "<br>", "Seq:", sequence)
            ),
            alpha = 0.9,
            size = 0.3,
            colour = "black"
          ) + 
          geom_text(
            data = local_table[local_table$type == 'Protein', ],
            aes(x = start/2 + end/2, y = (nivel - 0.3)/2 + (nivel + 0.3)/2), 
            label = protein_4
          ) +
          scale_fill_gradientn(
            colors = c("yellow", "red"),
            limits = c(0, max(local_table$count)),
            na.value = "lightblue"
          ) +
          theme_bw(base_size = 12) +
          theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            panel.border = element_blank(),
            legend.position = "none"
          ) +
          scale_y_continuous(expand = c(0, 0), limits = c(-5, 1.5))
        
        return(ggplotly(p, tooltip = c("text")) %>% layout(dragmode = "zoom", showlegend = FALSE))
      } else {
        # Normal case: protein with epitopes
        steps = c(5, 10, 25, 50, 100, 150, 2000, 250, 500, 1000, 1500, 2500, 5000)
        n_splits = 15
        max_x = local_table[local_table$type == "Protein", ]$end
        step_x_1 = max_x / n_splits
        step_x = steps[which(abs(steps - step_x_1) == min(abs(steps - step_x_1)))]
        
        p <- ggplot() +
          geom_rect(
            data = local_table, 
            aes(
              xmin = start,
              xmax = end,
              ymin = nivel - 0.3,
              ymax = nivel + 0.3,
              fill = count,
              text = paste(
                "Seq: ", sequence, 
                "<br>",
                "Position: ", start, "-", end - 1,
                "<br>",
                "Alleles: ", alleles
              )
            ),
            alpha = 0.9,
            size = 0.3,
            colour = "black"
          ) +
          geom_text(
            data = local_table[local_table$type == 'Protein', ],
            aes(x = start/2 + end/2, y = (nivel - 0.3)/2 + (nivel + 0.3)/2), 
            label = protein_4
          ) +
          scale_fill_gradientn(
            colors = c("yellow", "red"),
            limits = c(0, max(local_table$count)),
            na.value = "lightblue"
          ) +
          scale_y_continuous(expand = c(0, 0), limits = c(min(local_table$nivel) - 1, 1.5)) +
          scale_x_continuous(limits = c(-1, max_x), breaks = seq(0, max_x, by = step_x)) +
          theme_bw(base_size = 12) +
          theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            panel.border = element_blank(),
            legend.position = "none"
          )
        
        return(ggplotly(p, tooltip = c("text")) %>% layout(dragmode = "zoom", showlegend = FALSE))
      }
    })
    
    # Download button
    output$e_downloadData4 <- downloadHandler(
      filename = function() {
        paste("Epitopes_Plot_Location_example.txt")
      },
      content = function(file) {
        req(e_location_table())
        write.table(e_location_table(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
    
  })
}