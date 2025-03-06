# intersection.R - Módulo para análisis de intersección de epitopos
#
# Este módulo maneja el análisis de la intersección de epitopos entre múltiples alelos MHC,
# permitiendo identificar epitopos compartidos o específicos para combinaciones particulares.

#' UI para el módulo de intersección de epitopos
#'
#' @param id ID único para el módulo
#' @return UI para la pestaña de intersección de epitopos
intersection_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span('Epitope Intersection', style = 'font-size:120%'), br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          numericInput(
                            inputId = ns('cutoff_6'),
                            label = 'Cutoff %rank',
                            min = 0,
                            max = 10,
                            value = 2),
                          checkboxGroupInput(
                            inputId = ns('allele_6'),
                            label = 'MHC alleles',
                            choices = c(),
                            selected = c()),
                          radioButtons(
                            inputId = ns('plot_type_6'),
                            label = 'Plot Type',
                            choices = c('Venn Diagram', 'Up-Set'),
                            selected = 'Venn Diagram'),
                          actionButton(
                            inputId = ns('go_6'),
                            label = 'Run Analysis')
             ),
             mainPanel(width = 10,
                       h4("Set of epitopes shared between different MHC allele combinations"),
                       bsCollapse(bsCollapsePanel('Description (see more)',
                                                  wellPanel(style = 'overflow-x:auto; height:200px;', 
                                                           uiOutput(ns('text_Uset'))),
                                                  style = 'primary')),
                       wellPanel(
                         plotlyOutput(ns("plot_6"), height = "550px"),
                         plotOutput(ns("plot_6_VD"), height = "550px")
                       ),
                       wellPanel(fluidRow(column(12, withSpinner(DT::dataTableOutput(ns("table_6")), type = 5)))),
                       fluidRow(column(2, downloadButton(ns("downloadData6"), "Download Table"), offset = 6))
             )
           ))
}

#' Servidor para el módulo de intersección de epitopos
#'
#' @param id ID único para el módulo
#' @param parsed_data Reactive con datos analizados
#' @return NULL
intersection_server <- function(id, parsed_data) {
  moduleServer(id, function(input, output, session) {
    
    # Texto descriptivo
    output$text_Uset <- renderUI({
      list(p(HTML("Different populations are represented by distinct combinations
      of MHC alleles. This tool enables us to identify the set of epitopes that
      could be used in epitope vaccines that are potentially recognized by all
      MHC alleles in the population, as well as to identify epitopes restricted
      to a particular set of alleles. This tool shows the number of epitopes
      predicted to bind to different MHC allele combinations represented as a 
      Venn Diagram if 6 or fewer alleles are selected, or an Up-Set plot if 
      more than 6 alleles are selected. In addition to the downloadable plots,
      the tool provides a table containing the epitope sequences and the number
      of epitopes within each combination of alleles. In the Up-Set plot, the 
      selected MHC alleles are represented with bars on the left side indicating
      the number of predicted epitopes for each MHC allele. The number of epitopes
      for each MHC allele or intersection of MHC alleles is represented with bars
      at the top. Individual points in the grid indicate epitopes binding to a 
      specific MHC allele while connected points indicate epitopes that can bind
      to multiple MHC alleles. 
      <br>
      <h4> Parameters </h4>
      <ul>
       <li> <b>Cutoff %rank:</b> Maximum percentile rank to consider a peptide
       as a predicted epitope. Peptides with a percentile rank higher
       than the cutoff selected will be ignored. (Default MHC Class I = 2, MHC Class II = 10)</li>    
       <li> <b>MHC alleles:</b>  List of MHC alleles obtained from
       the input file. Users can choose more than one allele. (Default = First 3 alleles) </li>      
       </ul>")))
    })
    
    # Actualiza la lista de alelos cuando se cargan nuevos datos
    observe({
      req(parsed_data())
      alleles <- colnames(parsed_data())[c(-1, -2, -3, -4)]
      updateCheckboxGroupInput(session, inputId = "allele_6", 
                               choices = alleles, 
                               selected = alleles[1:min(length(alleles), 4)])
    })
    
    # Variables reactivas para inputs
    allele_6_react <- eventReactive(input$go_6, {input$allele_6})
    cutoff_6_react <- eventReactive(input$go_6, {input$cutoff_6})
    plot_type_6_react <- eventReactive(input$go_6, {input$plot_type_6})
    
    # Cálculo de datos de epitopos compartidos
    ep_data6 <- reactive({
      req(parsed_data(), allele_6_react())
      
      allele_6 <- allele_6_react()
      cutoff_6 <- cutoff_6_react()
      ep_data5 <- parsed_data()
      
      # Crea nuevo DataFrame para análisis
      new_data <- data.frame(matrix(ncol = length(allele_6) + 1, nrow = nrow(ep_data5)))
      colnames(new_data) <- c("Peptide", allele_6)
      new_data$Peptide <- ep_data5$Peptide
      
      # Marca alelos que cumplen con el criterio de cutoff
      for (allele in allele_6) {
        new_data[allele] <- ifelse(ep_data5[allele] <= cutoff_6, allele, NA)
      }
      
      # Elimina duplicados y procesa para obtener información de alelos
      new_data <- new_data[!duplicated(new_data$Peptide), ]
      new_data$Alleles = apply(new_data[, -1], 1, function(x) return(paste(x[!is.na(x)], collapse = '-')))
      new_data = new_data %>% filter(Alleles != '')
      
      # Agrupa por combinación de alelos
      new_data = new_data %>% 
        group_by(Alleles) %>% 
        summarize(Sequences = paste(Peptide, collapse = '-')) %>% 
        as.data.frame()
      
      # Calcula métricas adicionales
      new_data$N_epitopes = apply(new_data %>% select(Sequences), 1, 
                                  function(x) return(length(unlist(strsplit(x, split = '-')))))
      new_data$N_alleles = apply(new_data %>% select(Alleles), 1, 
                                function(x) return(length(unlist(strsplit(x, split = '-')))))
      
      # Reordena columnas
      new_data = new_data[, c(4, 1, 2, 3)]
      
      return(new_data)
    })
    
    # Renderiza tabla
    output$table_6 <- DT::renderDataTable({
      req(ep_data6())
      ep_data6() %>% arrange(desc(N_alleles))
    }, rownames = FALSE, options = list(scrollX = "300px"))
    
    # Gráfico Venn Diagram
    output$plot_6_VD <- renderPlot({
      req(parsed_data(), allele_6_react())
      
      allele_6 <- allele_6_react()
      cutoff_6 <- cutoff_6_react()
      
      # No mostrar si hay demasiados alelos
      if (length(allele_6) > 6) {
        return(NULL)
      }
      
      ep_data5 <- parsed_data()
      
      # Filtra datos para el diagrama
      ep_data6 = ep_data5 %>% select(1, all_of(allele_6))
      ep_data6 = melt(ep_data6, id.vars = 'Peptide') %>% as.data.frame()
      ep_data6 = ep_data6 %>% filter(value <= cutoff_6)
      
      # Crea y procesa diagrama de Venn
      myList = split(ep_data6$Peptide, ep_data6$variable)
      venn <- Venn(myList)
      data <- process_data(venn)
      
      # Genera gráfico
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
    
    # Gráfico Up-Set
    output$plot_6 <- renderPlotly({
      req(parsed_data(), allele_6_react())
      
      allele_6 <- allele_6_react()
      cutoff_6 <- cutoff_6_react()
      
      ep_data5 <- parsed_data()
      
      # Prepara datos para Up-Set plot
      new_data <- data.frame(matrix(ncol = length(allele_6) + 1, nrow = nrow(ep_data5)))
      colnames(new_data) <- c("Peptide", allele_6)
      new_data$Peptide <- ep_data5$Peptide
      
      for (allele in allele_6) {
        new_data[allele] <- ep_data5[allele] <= cutoff_6
        new_data[, allele] <- as.integer(new_data[, allele])
      }
      
      new_data <- new_data[order(new_data[, 'Peptide']), ]
      new_data <- new_data[!duplicated(new_data$Peptide), ]
      
      # Calcula intersecciones
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
      
      # Prepara datos para gráfico
      nintersections = nrow(z)
      nalleles = length(allele_6)
      
      set_size <- data.frame(alleles = allele_6, 
                             size = colSums(new_data[, 2:ncol(new_data)]))
      set_size <- set_size[order(set_size[, 'size']), ]
      
      ordered_alleles <- as.character(set_size$alleles)
      dict <- vector("list")
      i = length(ordered_alleles)
      for (allele in ordered_alleles) {
        dict[[allele]] <- i
        i = i - 1
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
          marker = list(color = "black", hoverinfo = "none")) %>% 
        layout(xaxis = list(tickfont = list(size = 15)),
               yaxis = list(tickfont = list(size = 15))) %>% 
        add_trace(type = "scatter",
                  mode = "text",
                  x = 1:nintersections,
                  y = zz$values + max(zz$values) * 0.05,
                  text = zz$values,
                  textfont = list(color = "black"))
      
      # Prepara datos para la cuadrícula
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
      
      # Gráfico de cuadrícula
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
    
    # Controla visibilidad de gráficos según tipo y número de alelos
    observe({
      allele_6 <- input$allele_6
      plot_type_6 <- input$plot_type_6
      
      if (length(allele_6) <= 6 & plot_type_6 == 'Venn Diagram') {
        hideElement(id = "plot_6")
        showElement(id = "plot_6_VD")
      } else {
        hideElement(id = "plot_6_VD")
        showElement(id = "plot_6")
      }
    })
    
    # Botón de descarga
    output$downloadData6 <- downloadHandler(
      filename = function() {
        paste("Up_Set_epitopes.txt")
      },
      content = function(file) {
        write.table(ep_data6(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      })
  })
}
