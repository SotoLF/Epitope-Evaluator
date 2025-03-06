# app.R - Archivo principal de la aplicación Epitope-Evaluator
#
# Este archivo integra todos los componentes de la aplicación y define
# la estructura general de la interfaz de usuario.

# Carga configuración global y librerías
source("global.R")

# Carga módulos UI
source("modules/data_input.R")
source("modules/distribution.R")
source("modules/intersection.R")
source("modules/density.R")
source("modules/viewer.R")
source("modules/promiscuity.R")
source("modules/conservation.R")
source("modules/example.R")

# Carga componentes UI
source("ui/ui_about.R")
source("ui/ui_documentation.R")
source("ui/ui_tutorial.R")

# UI principal
ui <- navbarPage(
  title = span(
    "Epitope-Evaluator",
    style = 'color: white; font-size:140%',
    tags$head(HTML("<title>Epitope-Evaluator</title>"))
  ), 
  collapsible = TRUE, 
  inverse = TRUE, 
  
  # Panel Home con análisis principal
  tabPanel(
    useShinyjs(), 
    title = span("Home", style = 'font-size:130%'),
    fluidPage(
      setBackgroundColor("#ecf0f5"),
      tabsetPanel(
        # Carga módulos
        data_input_ui("data_input"),
        distribution_ui("distribution"),
        intersection_ui("intersection"),
        density_ui("density"),
        viewer_ui("viewer"),
        promiscuity_ui("promiscuity"),
        conservation_ui("conservation")
      )
    )
  ),
  
  # Carga componentes UI adicionales
  about_ui(),
  example_ui("example"),
  documentation_ui(),
  tutorial_ui()
)

# Server principal
server <- function(session, input, output) {
  
  # Inicializa servidores de los componentes UI
  about_server(output)
  documentation_server(output)
  
  # Inicializa servidores de módulos
  parsed_data <- data_input_server("data_input", session)
  distribution_server("distribution", parsed_data)
  intersection_server("intersection", parsed_data)
  density_server("density", parsed_data)
  viewer_server("viewer", parsed_data)
  promiscuity_server("promiscuity", parsed_data)
  conservation_server("conservation", parsed_data)
  
  # Inicializa servidor de ejemplo
  example_server("example")
}

# Inicia la aplicación
shinyApp(ui = ui, server = server)
