# conservation.R - Módulo para análisis de conservación de epitopos
#
# Este módulo maneja el análisis de epitopos conservados entre múltiples proteínas,
# permitiendo identificar epitopos compartidos o específicos.

#' UI para el módulo de conservación de epitopos
#'
#' @param id ID único para el módulo
#' @return UI para la pestaña de conservación de epitopos
conservation_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Conservation', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          checkboxGroupInput(
                            inputId = ns('protein_7'),
                            label = 'Proteins',
                            choices = c(),
                            selected = c()),
                          numericInput(
                            inputId = ns('cutoff_7'),
                            label = 'Cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          
                          checkboxGroupInput(
                            inputId = ns('allele_7'),
                            label = 'MHC alleles',
                            choices = c(),
                            selected = c()),
                          
                          radioButtons(
                            inputId = ns('conditional_7'),
                            label = 'Shared epitopes',
                            choices = c('Intersection', 'Union'),
                            selected = 'Union'),
                          radioButtons(
                            inputId = ns('plot_type_7'),
                            label = 'Plot Type',
                            choices = list('Venn Diagram', 'Up-Set'),
                            selected = 'Venn Diagram'),
                          
                          actionButton(
                            inputId = ns('go_7'),
                            label = 'Run Analysis')
                          
             ),
             mainPanel(width = 10,
                       h4("Conservation epitopes across variants"),
                       bsCollapse(bsCollapsePanel('Description (see more)',
                                                  wellPanel(style = 'overflow-x:auto; height:200px;', 
                                                           uiOutput(ns('text_7'))),
                                                  style = 'primary')),
                       wellPanel(
                         plotlyOutput(ns("plot_7"), height = "550px"),
                         plotOutput(ns("plot_7_VD"), height = "550px")
                       ),
                       fluidRow(
                         column(12, withSpinner(DT::dataTableOutput(ns("table_7")), type = 5)),
                         column(2, downloadButton(ns('downloadData7'), 'Download Table'), offset = 6)
                       )
             )
           ))
}

#' Servidor para el módulo de conservación de epitopos
#'
#' @param id ID único para el módulo
#' @param parsed_data Reactive con datos analizados
#' @return NULL
conservation_server <- function(id, parsed_data) {
  moduleServer(id, function(input, output, session) {
    
    # Texto descriptivo
    output$text_7 <- renderUI({
      list(p(HTML("This tool allows users to identify epitopes that are present
      in multiple proteins, which can be useful to identify conserved epitopes
      across different pathogen strains. In addition, this tool allows for
      identifying epitopes gained or lost by diverse mutations. For this
      tool, users need to select the proteins, the MHC alleles of interest,
      and the cutoff percentile rank. The number of shared epitopes is
      represented as Venn Diagrams (<= 6 proteins) or Up-Set plot (> 6 proteins).
      <br>
      <h4> Parameters </h4>
      <ul>
       <li> <b>Protein:</b>  List of proteins obtained from the input file. Users should select the set of proteins
       of interest. (Default = First two proteins)</li>    
       <li> <b>Cutoff:</b> Maximum percentile rank to consider a peptide as a predicted epitope. (Default MHC Class I = 2, MHC Class II = 10) </li>      
       <li> <b>MHC alleles:</b> List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = First three alleles) </li>
       <li> <b>Shared epitopes:</b>  When multiple alleles are selected, users must indicate if they are
       interested in epitopes that bind to all the selected MHC alleles (intersection) or epitopes that
       bind to at least one of the MHC alleles selected (union).</li>
       </ul>")))
    })
    
    # Actualiza UI cuando cambian los datos
    observe({
      req(parsed_data())
      
      # Actualiza lista de proteínas
      proteins <- parsed_data() %>% dplyr::select(ID) %>% unlist() %>% unique()
      updateCheckboxGroupInput(session, "protein_7", 
                               choices = proteins, 
                               selected = proteins[1:min(3, length(proteins))])
      
      # Actualiza lista de alelos
      alleles <- colnames(parsed_data())[c(-1, -2, -3, -4)]
      updateCheckboxGroupInput(session, "allele_7", 
                               choices = alleles, 
                               selected = alleles[1:min(3, length(alleles))])
    })
    
    # Variables reactivas para inputs
    allele_7_react <- eventReactive(input$go_7, {input$allele_7})
    cutoff_7_react <- eventReactive(input$go_7, {input$cutoff_7})
    protein_7_react <- eventReactive(input$go_7, {input$protein_7})
    plot_type_7_react <- eventReactive(input$go_7, {input$plot_type_7})
    conditional_7_react <- eventReactive(input$go_7, {input$conditional_7})
    
    # Cálculo de tabla de epitopos
    table_epitopes <- reactive({
      req(parsed_data(), allele_7_react(), protein_7_react())
      
      parsed <- parsed_data()
      allele_7 <- allele_7_react()
      cutoff_7 <- cutoff_7_react()
      protein_7 <- protein_7_react()
      conditional_7 <- conditional_7_react()
      
      # Filtra por proteínas seleccionadas
      parsed <- parsed %>% filter(ID %in% protein_7)
      
      # Obtiene puntuación para alelos seleccionados
      parsed_1 <- parsed %>% select(Peptide, all_of(allele_7)) %>% unique()
      
      # Calcula puntuación según criterio seleccionado
      if (length(allele_7) != 1) {
        if (conditional_7 == 'Intersection') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, max)
        } else if (conditional_7 == 'Union') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, min) 
        }
      } else {
        parsed_1$plot_score <- parsed_1[, allele_7]
      }
      
      # Filtra péptidos candidatos
      candidate_peptides <- parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
      parsed_2 <- parsed %>% 
        select(Peptide, ID, all_of(allele_7)) %>% 
        filter(Peptide %in% candidate_peptides)
      
      # Prepara tabla de información
      peptide_info <- dcast(parsed_2[, c(1, 2)], formula = Peptide ~ ID) %>% as.data.frame()
      peptide_info$Proteins <- apply(peptide_info[, -1], 1, function(x) return(paste(x[!is.na(x)], collapse = '-')))
      
      # Agrupa por proteínas
      peptide_info <- peptide_info %>% 
        group_by(Proteins) %>%
        summarize(sequences = paste(Peptide, collapse = '-')) %>%
        as.data.frame()
      
      # Calcula número de secuencias
      peptide_info$N_sequences <- apply(peptide_info %>% select(sequences), 1, 
                                       function(x) return(length(unlist(strsplit(x, split = '-')))))
      
      return(peptide_info)
    })
    
    # Renderiza tabla
    output$table_7 <- DT::renderDataTable({
      req(table_epitopes())
      table_epitopes()
    }, rownames = FALSE, options = list(scrollX = '300px'))
    
    # Gráfico Venn Diagram
    output$plot_7_VD <- renderPlot({
      req(parsed_data(), allele_7_react(), protein_7_react())
      
      allele_7 <- allele_7_react()
      cutoff_7 <- cutoff_7_react()
      protein_7 <- protein_7_react()
      conditional_7 <- conditional_7_react()
      
      # No mostrar si hay demasiadas proteínas
      if (length(protein_7) > 6) {
        return(NULL)
      }
      
      # Filtra datos por proteínas seleccionadas
      parsed <- parsed_data() %>% filter(ID %in% protein_7)
      parsed_1 <- parsed %>% select(Peptide, all_of(allele_7)) %>% unique()
      
      # Calcula puntuación según criterio seleccionado
      if (length(allele_7) != 1) {
        if (conditional_7 == 'Intersection') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, max)
        } else if (conditional_7 == 'Union') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, min) 
        }
      } else {
        parsed_1$plot_score <- parsed_1[, allele_7]
      }
      
      # Filtra péptidos candidatos
      candidate_peptides <- parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
      parsed_2 <- parsed %>% 
        select(Peptide, ID, all_of(allele_7)) %>% 
        filter(Peptide %in% candidate_peptides)
      
      # Crea y procesa diagrama de Venn
      myList <- split(parsed_2$Peptide, parsed_2$ID)
      venn <- Venn(myList)
      data <- process_data(venn)
      
      # Genera gráfico
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
    
    # Gráfico Up-Set
    output$plot_7 <- renderPlotly({
      req(parsed_data(), allele_7_react(), protein_7_react())
      
      allele_7 <- allele_7_react()
      cutoff_7 <- cutoff_7_react()
      protein_7 <- protein_7_react()
      conditional_7 <- conditional_7_react()
      
      # Filtra datos por proteínas seleccionadas
      parsed <- parsed_data() %>% filter(ID %in% protein_7)
      parsed_1 <- parsed %>% select(Peptide, all_of(allele_7)) %>% unique()
      
      # Calcula puntuación según criterio seleccionado
      if (length(allele_7) != 1) {
        if (conditional_7 == 'Intersection') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, max)
        } else if (conditional_7 == 'Union') {
          parsed_1$plot_score <- apply(parsed_1[, allele_7, drop = FALSE], 1, min) 
        }
      } else {
        parsed_1$plot_score <- parsed_1[, allele_7]
      }
      
      # Filtra péptidos candidatos
      candidate_peptides <- parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
      parsed_2 <- parsed %>% 
        select(Peptide, ID, all_of(allele_7)) %>% 
        filter(Peptide %in% candidate_peptides)
      
      # Prepara datos para Up-Set plot
      new_data2 <- parsed_2 %>% dcast(Peptide ~ ID)
      new_data2[, -1] <- apply(new_data2[, -1], 2, function(x) return(as.integer(ifelse(is.na(x), 0, 1)))) %>% 
                         as.data.frame()
      
      # Normaliza nombres de columnas
      colnames(new_data2) <- gsub("\\|", "_", colnames(new_data2))
      
      # Lista de proteínas
      protein <- colnames(new_data2[, -1])
      
      # Calcula intersecciones
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
      
      # Prepara datos para gráfico
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
      
      # Gráfico de tamaños de conjuntos
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
      
      # Gráfico de tamaños de intersección
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
      
      # Prepara datos para la cuadrícula
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
      
      # Gráfico de cuadrícula
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
      
      # Ajusta título de eje Y
      intersect_size_chart <- intersect_size_chart %>% 
        layout(yaxis = list(title = "Intersections size"))
      
      # Combina subgráficos
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
      
      # Combina todo
      subplot(s1, s2, widths = c(0.3, 0.7))
    })
    
    # Controla visibilidad de gráficos según tipo y número de proteínas
    observe({
      protein_7 <- input$protein_7
      plot_type_7 <- input$plot_type_7
      
      if (length(protein_7) <= 6 & plot_type_7 == 'Venn Diagram') {
        hideElement(id = "plot_7")
        showElement(id = "plot_7_VD")
      } else {
        hideElement(id = "plot_7_VD")
        showElement(id = "plot_7")
      }
    })
    
    # Botón de descarga
    output$downloadData7 <- downloadHandler(
      filename = function() {
        paste("Conservation_epitopes.txt")
      },
      content = function(file) {
        write.table(table_epitopes(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
  })
}
