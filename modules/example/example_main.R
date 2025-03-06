# example_main.R - Main module for running example analysis
#
# This module coordinates the various example submodules and
# provides a demonstration of the application's capabilities
# using preconfigured example data.

#' UI for the example module
#'
#' @param id Unique ID for the module
#' @return UI for the example tab
example_ui <- function(id) {
  ns <- NS(id)
  
  tabPanel(title = span("Run Example", style = 'font-size:130%'),
           fluidPage(
             tabsetPanel(
               # Include all example tab panels
               example_data_ui(ns("data")),
               example_distribution_ui(ns("distribution")),
               example_intersection_ui(ns("intersection")),
               example_density_ui(ns("density")),
               example_viewer_ui(ns("viewer")),
               example_promiscuity_ui(ns("promiscuity")),
               example_conservation_ui(ns("conservation"))
             )
           ))
}

#' Server for the example module
#'
#' @param id Unique ID for the module
#' @return NULL
example_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Get example data (already loaded from global.R)
    e_data <- reactive({
      example_data
    })
    
    # Initialize all example submodules
    example_data_server("data", e_data)
    example_distribution_server("distribution", e_data)
    example_intersection_server("intersection", e_data) 
    example_density_server("density", e_data)
    example_viewer_server("viewer", e_data)
    example_promiscuity_server("promiscuity", e_data)
    example_conservation_server("conservation", e_data)
  })
}