# viewer.R - Módulo para visualización de posición de epitopos
#
# Este módulo maneja la visualización de la posición de epitopos dentro de proteínas,
# permitiendo identificar regiones enriquecidas con epitopos.

#' UI para el módulo de visualización de epitopos
#'
#' @param id ID único para el módulo
#' @return UI para la pestaña de visualización de epitopos
viewer_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Viewer', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          selectInput(
                            inputId = ns('protein_4'),
                            label = 'Protein',
                            choices = c(),
                            selected = c()),
                          
                          numericInput(
                            inputId = ns('cutoff_4'),
                            label = 'Cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          
                          checkboxGroupInput(
                            inputId = ns('allele_4'),
                            label = 'MHC alleles',
                            choices = c(),
                            selected = c()),
                          
                          radioButtons(
                            inputId = ns('conditional_4'),
                            label = 'Shared epitopes',
                            choices = c('Intersection', 'Union'),
                            selected = 'Union'),
                          
                          actionButton(
                            inputId = ns('go_4'),
                            label = 'Run Analysis')
                          
             ),
             mainPanel(width = 10,
                       h4("Amino acid position of the epitopes in proteins"),
                       bsCollapse(bsCollapsePanel('Description (see more)',
                                                  wellPanel(style = 'overflow-x:auto; height:200px;', 
                                                           uiOutput(ns('text_4'))),
                                                  style = 'primary')),
                       fluidRow(column(12, wellPanel(withSpinner(plotlyOutput(ns("plot_4")), type = 5)))),
                       fluidRow(column(12, wellPanel(withSpinner(DT::dataTableOutput(ns("table_4")), type = 5))),
                                column(2, downloadButton(ns('downloadData4'), 'Download Table'), offset = 6)))
           ))
}

#' Servidor para el módulo de visualización de epitopos
#'
#' @param id ID único para el módulo
#' @param parsed_data Reactive con datos analizados
#' @return NULL
viewer_server <- function(id, parsed_data) {
  moduleServer(id, function(input, output, session) {
    
    # Texto descriptivo
    output$text_4 <- renderUI({
      list(p(HTML("This tool can be used to identify the location of epitopes
      predicted to bind single and multiple MHC alleles, facilitating the visual
      identification of regions enriched with promiscuous epitopes. This tool 
      graphically represents the position of each predicted epitope within 
      each protein. Epitope Viewer shows the protein as a light blue bar 
      while epitopes are represented on a color gradient from yellow to red 
      reflecting the number of MHC alleles they are predicted to bind to.  
      Hovering on each epitope will show the amino acid sequence, the position 
      within the protein, and the MHC alleles that are predicted to recognize it.
      Users must select a cutoff to identify epitopes and the MHC alleles of
      interest. Choosing the 'Intersection' option will show only epitopes
      predicted to bind to all of the selected MHC alleles. Selecting the 
      'Union' option displays all the epitopes which are colored ranging from
      yellow to red indicating the number of MHC alleles predicted to bind. 
      This tool also shows a table containing the sequence of the epitope, 
      the start and end amino acid position within the protein, the alleles,
      and the number of alleles to which they were predicted to bind. This
      analysis may take a few minutes depending on the number of epitopes and
      the length of the protein.

      <br>
      <h4> Parameters </h4>
      <ul>
       <li> <b>Protein:</b> List of proteins obtained from the input file. Users should select the protein of interest. (Default = First protein)</li>    
       <li> <b>Cutoff % rank:</b> Maximum percentile rank to consider a peptide as a predicted epitope. (Default MHC Class I = 2, MHC Class II = 10) </li>      
       <li> <b>Shared epitopes:</b> When multiple alleles are selected, users must indicate if they are interested in epitopes
       that bind to all the selected MHC alleles (intersection) or epitopes that bind to at least one of the MHC alleles selected (union).</li>
       <li> <b>MHC alleles:</b>  List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = first allele) </li>
       </ul>")))
    })
    
    # Actualiza UI cuando cambian los datos
    observe({
      req(parsed_data())
      # Actualiza lista de proteínas
      updateSelectInput(session, "protein_4", 
                        choices = parsed_data() %>% dplyr::select(ID) %>% unlist() %>% unique(), 
                        selected = parsed_data()["ID"][1, ])
      
      # Actualiza lista de alelos
      alleles <- colnames(parsed_data())[c(-1, -2, -3, -4)]
      updateCheckboxGroupInput(session, "allele_4", 
                               choices = alleles, 
                               selected = alleles[1])
    })
    
    # Variables reactivas para inputs
    allele_4_react <- eventReactive(input$go_4, {input$allele_4})
    cutoff_4_react <- eventReactive(input$go_4, {input$cutoff_4})
    protein_4_react <- eventReactive(input$go_4, {input$protein_4})
    conditional_4_react <- eventReactive(input$go_4, {input$conditional_4})
    
    # Cálculo de tabla de ubicación
    location_table <- reactive({
      req(parsed_data(), allele_4_react(), protein_4_react())
      
      ep_data5 <- parsed_data()
      allele_4 <- allele_4_react()
      cutoff_4 <- cutoff_4_react()
      protein_4 <- protein_4_react()
      conditional_4 <- conditional_4_react()
      
      # Filtra por proteína y selecciona datos relevantes
      peptide_counts = parsed_data() %>% 
        filter(ID == protein_4) %>% 
        select(Peptide, Pos, Length, all_of(allele_4))
      
      # Calcula puntuación según criterio seleccionado
      if (length(allele_4) != 1) {
        if (conditional_4 == "Intersection") {
          peptide_counts$score <- apply(peptide_counts[, allele_4, drop = FALSE], 1, max)
        } else if (conditional_4 == "Union") {
          peptide_counts$score <- apply(peptide_counts[, allele_4, drop = FALSE], 1, min)
        }
      } else {
        peptide_counts$score <- peptide_counts[, allele_4]
      }
      
      # Filtra y procesa datos
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
      
      # Comprueba si hay resultados
      if (nrow(peptide_counts) == 0) return(NULL)
      
      # Obtiene detalles de los péptidos
      ep_data2 = parsed_data() %>%
        filter(ID == protein_4) %>% 
        select(Peptide, Pos, Length)
      
      # Combina información
      index = match(ep_data2$Peptide, peptide_counts$Peptide)
      ep_data2$Value = peptide_counts$count[index]
      ep_data2$alleles = peptide_counts$alleles[index]
      ep_data2 <- ep_data2 %>% filter(!is.na(Value) & Value > 0)
      
      # Prepara tabla para visualización
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
      
      # Algoritmo para organizar visualización en niveles
      nivel = 0
      while (nivel %in% local_table$nivel) {
        subdata <- local_table[local_table$type == "Epitope" & local_table$nivel == nivel, ]
        changed = TRUE
        
        while (changed) {
          number_rows = nrow(subdata)
          
          # Caso especial: solo un epitopo
          if (number_rows == 1) {
            changed = FALSE
            nivel = nivel - 1
            break
          }
          
          conteo = 0
          for (i in 1:number_rows) {
            # Fila estática
            static_row <- subdata[i, ]
            # Vector estático (posiciones)
            v_s <- seq(static_row$start, static_row$end)
            new_data <- subdata[-i, ]
            
            # Compara con todos los demás epitopos
            conteo = conteo + 1
            for (row in 1:nrow(new_data)) {
              row_data <- new_data[row, ]
              v_x <- seq(row_data$start, row_data$end)
              id <- row_data$id
              
              # Si hay solapamiento, ajusta nivel
              if (length(intersect(v_s, v_x)) != 0) {
                local_table[local_table$id == id, ]$nivel <- row_data$nivel - 1
                subdata[subdata$id == id, ]$nivel <- row_data$nivel - 1
              }
            }
            
            # Actualiza subdatos con nuevo nivel
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
    
    # Renderiza tabla
    output$table_4 <- DT::renderDataTable({
      req(location_table())
      
      local_table <- location_table()
      if (is.null(local_table)) {
        return(data.frame(
          sequence = character(),
          start = integer(),
          end = integer(),
          alleles = character(),
          count = integer()
        ))
      }
      
      # Excluye fila de proteína y selecciona columnas relevantes
      local_table <- local_table[-1, c(4, 2, 3, 7, 8)]
      
      DT::datatable(local_table, rownames = FALSE, options = list(scrollX = "300px"))
    })
    
    # Gráfico de visualización
    output$plot_4 <- renderPlotly({
      req(location_table())
      
      ep_data5 <- parsed_data()
      protein_4 <- protein_4_react()
      local_table <- location_table()
      
      if (is.null(local_table)) return(NULL)
      
      # Convierte count a numérico para color
      local_table$count <- as.numeric(local_table$count) 
      
      # Caso especial: solo una proteína (sin epitopos)
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
            label = 'protein_4'
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
        # Caso normal: proteína con epitopos
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
    
    # Botón de descarga
    output$downloadData4 <- downloadHandler(
      filename = function() {
        paste("Epitopes_Plot_Location.txt")
      },
      content = function(file) {
        write.table(location_table(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
  })
}
