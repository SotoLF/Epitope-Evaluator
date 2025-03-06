# global.R - Configuración global y carga de librerías para Epitope-Evaluator
#
# Este archivo contiene todas las librerías necesarias, opciones de configuración
# y carga los archivos de utilidades necesarios para la aplicación.

# Carga de librerías ----
library(shiny)
library(shinythemes)
library(plotly)
library(shinyjs)
library(shinycssloaders)
library(shinyBS)
library(rsconnect)

# Librerías para visualización y manipulación de datos
library(ggplot2)
library(dplyr)
library(readr)
library(grid)
library(gridExtra)
library(reshape)
library(shinydashboard)
library(tidyselect)
library(rlist)
library(tibble)
library(seqinr)
library(phylotools)
library(reshape2)
library(ggVennDiagram)
library(stringr)
library(shinyWidgets)

# Configuración global de la aplicación ----
options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=1)
options(shiny.maxRequestSize = 30*1024^2)  # Permite subir archivos de hasta 30MB

# Carga de archivos de utilidades
source("utils/parsing_functions.R")
source("utils/helper_functions.R")

# Carga del conjunto de datos de ejemplo
example_data = Parse_NetMHCPAN('data/example.xls', 'data/example.fasta', 'Rank')
colnames(example_data)[1:4] = c("Peptide", "Pos", "Length", "ID")  # Estandariza nombres de columnas
