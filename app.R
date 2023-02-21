#setwd('D:/BostonLab/ShinyApp/scripts')
#source("Functions.R")

library(shiny)
library(shinythemes)
library(plotly)
library(shinyjs)
library(shinycssloaders)
library(shinyBS)
library(rsconnect)

library(ggplot2)
library(dplyr)
library(readr)
library(grid)
library(gridExtra)
library(reshape)
library(shinydashboard)
library(tidyselect)
library(rlist)
library(dplyr)
library(tibble)
library(seqinr)
library(phylotools)
library(reshape2)
library(ggVennDiagram)
library(stringr)
library(shinyWidgets)


options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=1)
options(shiny.maxRequestSize = 30*1024^2)

get_intersect_members <- function (x, ...){
  
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}

Parse_NetMHCPAN <- function(prediction_file, fasta_file, score_type){
  
  ## Read prediction file
  
  results <- data.frame(read.table(prediction_file, 
                                   sep = "\t"))
  
  ## Get length from fasta file
  faas<- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = T))-8
  
  ## Get names from fasta file
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"),sep = " ")
  
  ## Create metadata DF
  metadata_df$length <- lengths
  metadata_df$prev <- results[c(-1,-2),3] %>% unique() %>% unlist()
  
  ## Get alleles names for Class I
  raw_header <- results[1,]
  alleles_ID <- raw_header[grepl("HLA", raw_header)]
  
  ## Add names to metadata DF
  results$ID_2 <- metadata_df$ID[match(results[,3], metadata_df$prev)]
  results$length <- metadata_df$length[match(results[,3], metadata_df$prev)]
  
  ## Prediction adding Binding affinity Score (BA)
  if((ncol(results) - 7) / length(alleles_ID) == 6){
    if (score_type == "Rank"){
      # Get index for % Ranks Class I
      col_IDs <- seq(from = 7, by = 6, length.out = length(alleles_ID))
    } else{
      # Index for nM concentration Class I
      col_IDs <- seq(from = 9, by = 6, length.out = length(alleles_ID))
    }
  } else{ ## Prediction without Binding affinity Score (BA)
    if (score_type == "Rank"){
      # Get index for % Ranks Class I
      col_IDs <- seq(from = 7, by = 4, length.out = length(alleles_ID))
    } else{
      # Get index for % Ranks Class I (There is no BA column in this file type)
      col_IDs <- seq(from = 7, by = 3, length.out = length(alleles_ID))
    }
  }
  
  ## Create parsed DF
  df_parsed <- results[c(-1,-2), c(2, 1, ncol(results), ncol(results)-1, col_IDs)]
  colnames(df_parsed) <- c("Peptide", "Position", "Length", "ID", alleles_ID)
  df_parsed$Position = as.numeric(df_parsed$Position) + 1
  df_parsed[c(5:ncol(df_parsed))] <-sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)
  colnames(df_parsed) = str_replace_all(colnames(df_parsed), ':', '.') %>% str_replace_all('-','.')
    
  return(df_parsed)
}
Parse_NetMHC <- function(prediction_file, fasta_file, score_type){
  ## Read prediction file
  results <- data.frame(read.table(prediction_file,
                                   sep = "\t", row.names = NULL))
  
  ## Get length from fasta file
  faas<- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = T))-8
  
  ## Get names from fasta file
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"),sep = " ")
  
  ## Create metadata DF
  metadata_df$length <- lengths
  metadata_df$prev <- results[-1,3] %>% unique() %>% unlist()
  
  ## Get alleles names for Class I
  raw_header <- colnames(results)
  alleles_ID <- raw_header[grepl("HLA", raw_header)]
  
  ## Add names to metadata DF
  results$ID_2 <- metadata_df$ID[match(results[,3], metadata_df$prev)]
  results$length <- metadata_df$length[match(results[,3], metadata_df$prev)]
  
  if (score_type == "Rank"){
    # Get index for % Ranks Class I
    col_IDs <- seq(from = 5, by = 3, length.out = length(alleles_ID))
  } else{
    # Index for nM concentration Class I
    col_IDs <- seq(from = 4, by = 3, length.out = length(alleles_ID))
  }
  
  ## Create parsed DF
  df_parsed <- results[-1, c(2, 1, ncol(results), ncol(results)-1, col_IDs)]
  colnames(df_parsed) <- c("Peptide", "Position", "Length", "ID", alleles_ID)
  df_parsed$Position = as.numeric(df_parsed$Position) + 1
  
  df_parsed[c(5:ncol(df_parsed))] <-sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)  
  return(df_parsed)
}
Parse_NetMHCIIPAN <- function(prediction_file, fasta_file, score_type){
  
  ## Read prediction file
  
  results <- data.frame(read.table(prediction_file,
                                   sep = ",", header = FALSE))
  
  ## Get alleles names for Class II
  alleles_vector = results[1,] %>% strsplit('\t') %>% unlist()
  alleles_ID = alleles_vector[grepl('D', alleles_vector)]
  
  # Parse Data
  results = results[c(-1,-2),] %>% strsplit('\t')
  results = as.data.frame(do.call(rbind, results))
  
  ## Get length from fasta file
  faas<- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = T))-8
  
  ## Get names from fasta file
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"),sep = " ")
  
  ## Create metadata DF
  metadata_df$length <- lengths
  metadata_df$prev <- results[,3] %>% unique() %>% unlist()
  
  ## Add names to metadata DF
  results$ID_2 <- metadata_df$ID[match(results[,3], metadata_df$prev)]
  results$length <- metadata_df$length[match(results[,3], metadata_df$prev)]
  
  ## Prediction adding Binding affinity Score (BA)
  if((ncol(results) - 8) / length(alleles_ID) == 6){
    if (score_type == "Rank"){
      # Get index for % Ranks Class II
      col_IDs <- seq(from = 7, by = 6, length.out = length(alleles_ID))
    } else{
      # Index for nM concentration Class II
      col_IDs <- seq(from = 10, by = 6, length.out = length(alleles_ID))
    }
  } else{ ## Prediction without Binding affinity Score (BA)
    if (score_type == "Rank"){
      # Get index for % Ranks Class II
      col_IDs <- seq(from = 7, by = 3, length.out = length(alleles_ID))
    } else{
      # Get index for % Ranks Class II (There is no BA column in this file type)
      col_IDs <- seq(from = 7, by = 3, length.out = length(alleles_ID))
    }
  }
  
  ## Create parsed DF
  df_parsed <- results[, c(2, 1, ncol(results), ncol(results)-1, col_IDs)]
  colnames(df_parsed) <- c("Peptide", "Position", "Length", "ID", alleles_ID)
  df_parsed$Position = as.numeric(df_parsed$Position)
  df_parsed[c(5:ncol(df_parsed))] <-sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)  
  
  return(df_parsed)
}
Parse_MHCFlurry <- function(prediction_file, fasta_file, score_type){
  ## Read prediction file
  results <- data.frame(read_delim(prediction_file,
                                   delim = ",", escape_double = FALSE, 
                                   trim_ws = TRUE))
  
  if (score_type == "Rank"){
    # Get index for % Ranks Class I
    results <- results[, c(3,2,1,8,9)]
  } else{
    # Index for nM concentration Class I
    results <- results[, c(3,2,1,8,7)]
  }
  
  results<-reshape2::dcast(data = results, formula = peptide +pos + sequence_name ~best_allele)
  
  ## Get length from fasta file
  faas<- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = T))-8
  
  ## Get names from fasta file
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"),sep = " ")
  
  ## Create metadata DF
  metadata_df$length <- lengths
  
  ## Add names to metadata DF
  results$length <- metadata_df$length[match(results[,3], metadata_df$ID)]
  
  ## Create parsed DF
  df_parsed <- results[c(1, 2, ncol(results), 3, 4:(ncol(results)-1))]
  colnames(df_parsed)[1:4] <- c("Peptide", "Position", "Length", "ID")
  df_parsed$Position = as.numeric(df_parsed$Position) + 1
  df_parsed[c(5:ncol(df_parsed))] <-sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)  
  
  return(df_parsed)
}
Parse_IEDB <- function(prediction_file, fasta_file, score_type){
  
  results <- read.delim(prediction_file, sep = '\t')
  
  if (score_type == "Rank"){
    # Get index for % Ranks Class I
    results <- results[, c(6, 3, 2, 1, 7)]
  } else{
    # Index for nM concentration Class I
    results <- results[, c(6, 3, 2, 1, 8)]
  }
  
  results <- reshape2::dcast(data = results, formula = peptide + start + seq_num ~ allele)
  
  ## Get length from fasta file
  faas<- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = T)) - 8
  
  ## Get names from fasta file
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"),sep = " ")
  
  ## Create metadata DF
  metadata_df$length <- lengths
  
  results$seq_num <- metadata_df[results$seq_num,1]  
  results$length <- metadata_df$length[match(results$seq_num, metadata_df$ID)]
  
  ## Create parsed DF
  df_parsed <- results[c(1, 2, ncol(results), 3, 4:(ncol(results)-1))]
  colnames(df_parsed)[1:4] <- c("Peptide", "Position", "Length", "ID")
  df_parsed[c(5:ncol(df_parsed))] <-sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)  
  
  return(df_parsed)
}



example_data = Parse_NetMHCPAN('example.xls','example.fasta', 'Rank')
colnames(example_data)[1:4] = c("Peptide", "Pos", "Length", "ID")


ui <- navbarPage(title = span("Epitope-Evaluator",
                              style ='color: white; font-size:140%',
                              tags$head(HTML("<title>Epitope-Evaluator</title>"))), collapsible = TRUE, inverse = TRUE, 
             tabPanel(useShinyjs(), title = span("Home", style ='font-size:130%'),
                      fluidPage( setBackgroundColor("#ecf0f5"),
                        tabsetPanel(
                          tabPanel(title = span('Input Data', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                       fileInput(
                                         inputId = 'input_file',
                                         label = 'Upload Prediction file'
                                       ),
                                       fileInput(
                                         inputId = 'input_fasta',
                                         label = 'Upload Fasta file',
                                       ),
                                       radioButtons(
                                         inputId = 'cond_pred',
                                         label = 'Predictor',
                                         choices = list('NetMHC', 'NetMHCpan', 'NetMHCIIpan', 'MHCFlurry', 'IEDB Consensus', 'Other'),
                                         selected = 'NetMHCpan'
                                       ),
                                       radioButtons(
                                         inputId = 'method',
                                         label = 'Score Type',
                                         choices = list('Score', 'Rank'),
                                         selected = 'Rank'
                                       ),
                                       actionButton(
                                         inputId = 'go_1',
                                         label = 'Run Analysis')
                                     ),
                                     mainPanel(width = 10, wellPanel(DT::dataTableOutput('contents_1')))
                                   )),
                          
                          tabPanel(title = span('Epitope Distribution', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  checkboxGroupInput(
                                                    inputId = 'allele',
                                                    label = 'MHC alleles',
                                                    choices = c(),
                                                    selected = c()),
                                                  radioButtons(
                                                    inputId = 'conditional',
                                                    label = 'Shared epitopes',
                                                    choices = list('Intersection', 'Union')),
                                                  numericInput(
                                                    inputId = 'xmin',
                                                    label = 'Min %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 0),
                                                  numericInput(
                                                    inputId = 'xmax',
                                                    label = 'Max %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  numericInput(
                                                    inputId = 'step',
                                                    label = 'Bin width',
                                                    min = 0,
                                                    max = 1,
                                                    value = 0.1),
                                                  radioButtons(
                                                    inputId = 'conditional_plot',
                                                    label = 'Plot type',
                                                    choices = list('Histogram', 'Cumulative Histogram')),
                                                  actionButton(
                                                    inputId = 'go_2_1',
                                                    label = 'Run analysis')
                                                  ),
                                     mainPanel(width = 10,
                                               h4("Distribution of epitopes by MHC allele"),
                                               bsCollapse(bsCollapsePanel('Description (see more)', 
                                                                          wellPanel(style = "overflow-x:auto; height:200px;",uiOutput('text_distribution')),
                                                                          style = "primary")),
                                               fluidRow(
                                                 column(5, wellPanel(withSpinner(plotlyOutput('distribution_plot', height = 500), type = 5))),
                                                 column(7, wellPanel(withSpinner(DT::dataTableOutput('distribution_table'), type = 5))),
                                                 column(2, downloadButton('downloadData2', 'Download Table'), offset = 6),
                                                 column(12, wellPanel(withSpinner(plotOutput('heatmap_plot', height = 500), type = 5)))
                                     )))),
                          
                          tabPanel(title = span('Epitope Intersection', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  
                                                  numericInput(
                                                    inputId = 'cutoff_6',
                                                    label = 'Cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  checkboxGroupInput(
                                                    inputId = 'allele_6',
                                                    label = 'MHC alleles',
                                                    choices = c(),
                                                    selected = c()),
                                                  
                                                  radioButtons(
                                                    inputId = 'plot_type_6',
                                                    label = 'Plot Type',
                                                    choices = c('Venn Diagram', 'Up-Set'),
                                                    selected = 'Venn Diagram'),
                                                  
                                                  actionButton(
                                                    inputId = 'go_6',
                                                    label = 'Run Analysis')
                                                  ),
                                     mainPanel(width = 10,
                                               h4("Set of epitopes shared between different MHC allele combinations"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('text_Uset')),
                                                                          style = 'primary')),
                                               wellPanel(plotlyOutput("plot_6", height = "550px"),
                                                         plotOutput("plot_6_VD", height = "550px")),
                                               wellPanel(fluidRow(column(12, withSpinner(DT::dataTableOutput("table_6"), type = 5)))),
                                               fluidRow(column(2, downloadButton("downloadData6", "Download Table"), offset = 6))
                                               
                                     )
                                   )
                          ),
                          
                          tabPanel(title = span('Epitope Density', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  'Density Analysis Parameters:',
                                                  
                                                  numericInput(
                                                    inputId = 'cutoff',
                                                    label = 'Cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  actionButton(
                                                    inputId = 'go_3_1',
                                                    label = 'Run Analysis'),
                                                  
                                                  br(),
                                                  br(),
                                                  'Plot parameters:',
                                                  
                                                  radioButtons(
                                                    inputId = 'plot_type',
                                                    label = 'Plot Type',
                                                    choices = c('Heatmap', 'Bar plot'),
                                                    selected = 'Heatmap'),
                                                  
                                                  radioButtons(
                                                    inputId = 'sort_type',
                                                    label = 'Sort Type',
                                                    choices = c('Default', 'Descendent'),
                                                    selected = 'Descendent'),
                                                  
                                                  radioButtons(
                                                    inputId = 'fill_type',
                                                    label = 'Color by:',
                                                    choices = c('Number of Epitopes', 'Density of Epitopes'),
                                                    selected = 'Number of Epitopes'),
                                                  
                                                  actionButton(
                                                    inputId = 'go_3_2',
                                                    label = 'Generate Plot')
                                                  
                                                  ),
                                     mainPanel(width = 10,
                                               h4("Correlation between number of epitopes and protein length"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('text_density')),
                                                                          style = 'primary')),
                                               fluidRow(column(6, wellPanel(withSpinner(plotlyOutput("plot_3_1", height = 500), type = 5))),
                                                        column(6, wellPanel(withSpinner(DT::dataTableOutput("table_3"), type = 5))),
                                                        column(2,downloadButton("downloadData3", "Download Table"), offset = 6)),
                                               fluidRow(column(12, wellPanel(withSpinner(plotlyOutput("plot_3_2", height = 700), type = 5))))
                                               )
                                   )
                                   ),
                          tabPanel(title = span('Epitope Viewer', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  
                                                  selectInput(
                                                    inputId = 'protein_4',
                                                    label = 'Protein',
                                                    choices = c(),
                                                    selected = c()),
                                                  
                                                  numericInput(
                                                    inputId = 'cutoff_4',
                                                    label = 'Cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  checkboxGroupInput(
                                                    inputId = 'allele_4',
                                                    label = 'MHC alleles',
                                                    choices = c(),
                                                    selected = c()),
                                                  
                                                  radioButtons(
                                                    inputId = 'conditional_4',
                                                    label = 'Shared epitopes',
                                                    choices = c('Intersection', 'Union'),
                                                    selected = 'Union'),
                                                  
                                                  actionButton(
                                                    inputId = 'go_4',
                                                    label = 'Run Analysis')
                                                  
                                                  ),
                                     mainPanel(width = 10,
                                               h4("Amino acid position of the epitopes in proteins"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('text_4')),
                                                                          style = 'primary')),
                                               fluidRow(column(12, wellPanel(withSpinner(plotlyOutput("plot_4"), type = 5)))),
                                               fluidRow(column(12, wellPanel(withSpinner(DT::dataTableOutput("table_4"), type = 5))),
                                                        column(2,downloadButton("downloadData4", "Download Table"), offset = 6)))
                                     )
                                   ),
                          
                          tabPanel(title = span('Epitope Promiscuity', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  
                                                  'Set Promiscuity Analysis Parameters:',
                                                  br(),
                                                  numericInput(
                                                    inputId = 'min_al',
                                                    label = 'Minimum number of alleles',
                                                    min = 0,
                                                    max = 10,
                                                    value = 6),
                                                  
                                                  numericInput(
                                                    inputId = 'sb_ct',
                                                    label = 'Strong Binding cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 0.5),
                                                  
                                                  numericInput(
                                                    inputId = 'wb_ct',
                                                    label = 'Weak Binding cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  actionButton(
                                                    inputId = 'go_5',
                                                    label = 'Run Analysis')
                                                  
                                                  ),
                                     mainPanel(width = 10,
                                               h4("Number of alleles and binding rank per epitopes"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('text_5')),
                                                                          style = 'primary')),
                                               fluidRow(column(6, wellPanel(withSpinner(plotlyOutput("plot_5", height = 500), type = 5))),
                                                        column(6, wellPanel(withSpinner(DT::dataTableOutput("table_5"), type = 5)))),
                                               fluidRow(column(2, downloadButton("downloadData5", "Download Table"), offset = 6))
                                               )
                                   )),
                          tabPanel(title = span('Epitope Conservation', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  
                                                  checkboxGroupInput(
                                                    inputId = 'protein_7',
                                                    label = 'Proteins',
                                                    choices = c(),
                                                    selected = c()),
                                                  numericInput(
                                                    inputId = 'cutoff_7',
                                                    label = 'Cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  checkboxGroupInput(
                                                    inputId = 'allele_7',
                                                    label = 'MHC alleles',
                                                    choices = c(),
                                                    selected = c()),
                                                  
                                                  radioButtons(
                                                    inputId = 'conditional_7',
                                                    label = 'Shared epitopes',
                                                    choices = c('Intersection', 'Union'),
                                                    selected = 'Union'),
                                                  radioButtons(
                                                    inputId = 'plot_type_7',
                                                    label = 'Plot Type',
                                                    choices = list('Venn Diagram', 'Up-Set'),
                                                    selected = 'Venn Diagram'),
                                                  
                                                  actionButton(
                                                    inputId = 'go_7',
                                                    label = 'Run Analysis')
                                                  
                                                  ),
                                     mainPanel(width = 10,
                                               h4("Conservation epitopes across variants"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('text_7')),
                                                                          style = 'primary')),
                                               wellPanel(plotlyOutput("plot_7", height = "550px"),
                                                         plotOutput("plot_7_VD", height = "550px")),
                                               fluidRow(column(12, withSpinner(DT::dataTableOutput("table_7"), type = 5)),
                                                        column(2, downloadButton('downloadData7', 'Download Table'), offset = 6))
                                               )
                                   )))
                      )),
             tabPanel(title = span("About", style ='font-size:130%'),
                      sidebarLayout(sidebarPanel(width = 4,
                                                 h3('Epitope-Evaluator'),
                                                 wellPanel(style= 'overflow-x:auto; height:70;', uiOutput('text_1_left'))), 
                                    mainPanel(width = 8,
                                              h1('What is Epitope-Evaluator'),
                                              wellPanel(uiOutput('text_1'),imageOutput('image_1',width = '100%', height = '100%'))
                                              )
                      )
                      ),
             tabPanel(title = span("Run Example", style ='font-size:130%'),
                      fluidPage(
                        tabsetPanel(
                          tabPanel(title = span('Example Data', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 4, h3('Input Data'),
                                                  wellPanel(uiOutput('e_text_1'))),
                                     mainPanel(width = 8, downloadButton('e_downloadData_example', 'Download Table'),
                                               wellPanel(DT::dataTableOutput('e_contents_1')))
                                   )),
                          
                          tabPanel(title = span('Epitope Distribution', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  checkboxGroupInput(
                                                    inputId = 'e_allele',
                                                    label = 'MHC alleles',
                                                    choices = colnames(example_data)[c(-1, -2, -3, -4)],
                                                    selected = colnames(example_data)[c(-1, -2, -3, -4)][c(1,2,3)]),
                                                  radioButtons(
                                                    inputId = 'e_conditional',
                                                    label = 'Shared epitopes',
                                                    choices = list('Intersection', 'Union'),
                                                    selected = 'Union'),
                                                  numericInput(
                                                    inputId = 'e_xmin',
                                                    label = 'Min %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 0),
                                                  numericInput(
                                                    inputId = 'e_xmax',
                                                    label = 'Max %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  numericInput(
                                                    inputId = 'e_step',
                                                    label = 'Bin width',
                                                    min = 0,
                                                    max = 1,
                                                    value = 0.1),
                                                  radioButtons(
                                                    inputId = 'e_conditional_plot',
                                                    label = 'Plot type',
                                                    choices = list('Histogram', 'Cumulative Histogram'),
                                                    selected = 'Histogram'),
                                                  actionButton(
                                                    inputId = 'e_go_2_1',
                                                    label = 'Run analysis')
                                     ),
                                     mainPanel(width = 10,
                                               h4("Distribution of epitopes by MHC allele"),
                                               bsCollapse(bsCollapsePanel('Description (see more)', 
                                                                          wellPanel(style = "overflow-x:auto; height:200px;",uiOutput('e_text_distribution')),
                                                                          style = "primary")),
                                               fluidRow(
                                                 column(5, wellPanel(withSpinner(plotlyOutput('e_distribution_plot', height = 500), type = 5))),
                                                 column(7, wellPanel(withSpinner(DT::dataTableOutput('e_distribution_table'), type = 5))),
                                                 column(2, downloadButton('e_downloadData2', 'Download Table'), offset = 6),
                                                 column(12, wellPanel(withSpinner(plotOutput('e_heatmap_plot', height = 500), type = 5)))
                                               )))),
                          
                          tabPanel(title = span('Epitope Intersection', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  
                                                  numericInput(
                                                    inputId = 'e_cutoff_6',
                                                    label = 'Cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  checkboxGroupInput(
                                                    inputId = 'e_allele_6',
                                                    label = 'MHC alleles',
                                                    choices = colnames(example_data)[c(-1, -2, -3, -4)],
                                                    selected = colnames(example_data)[c(-1, -2, -3, -4)][c(1,2,3,4)]),
                                                  
                                                  radioButtons(
                                                    inputId = 'e_plot_type_6',
                                                    label = 'Plot Type',
                                                    choices = c('Venn Diagram', 'Up-Set'),
                                                    selected = 'Venn Diagram'),
                                                  
                                                  actionButton(
                                                    inputId = 'e_go_6',
                                                    label = 'Run Analysis')
                                     ),
                                     mainPanel( width = 10,
                                               h4("Set of epitopes shared between different MHC allele combinations"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('e_text_Uset')),
                                                                          style = 'primary')),
                                               wellPanel(plotlyOutput("e_plot_6", height = "550px"),
                                                         plotOutput("e_plot_6_VD", height = "550px")),
                                               wellPanel(fluidRow(column(12, withSpinner(DT::dataTableOutput("e_table_6"), type = 5)))),
                                               fluidRow(column(2, downloadButton("e_downloadData6", "Download Table"), offset = 6))
                                               
                                     )
                                   )
                          ),
                          tabPanel(title = span('Epitope Density', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  'Density Analysis Parameters:',
                                                  
                                                  numericInput(
                                                    inputId = 'e_cutoff',
                                                    label = 'Cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  actionButton(
                                                    inputId = 'e_go_3_1',
                                                    label = 'Run Analysis'),
                                                  
                                                  br(),
                                                  br(),
                                                  'Plot parameters:',
                                                  
                                                  radioButtons(
                                                    inputId = 'e_plot_type',
                                                    label = 'Plot Type',
                                                    choices = c('Heatmap', 'Bar plot'),
                                                    selected = 'Heatmap'),
                                                  
                                                  radioButtons(
                                                    inputId = 'e_sort_type',
                                                    label = 'Sort Type',
                                                    choices = c('Default', 'Descendent'),
                                                    selected = 'Descendent'),
                                                  
                                                  radioButtons(
                                                    inputId = 'e_fill_type',
                                                    label = 'Color by:',
                                                    choices = c('Number of Epitopes', 'Density of Epitopes'),
                                                    selected = 'Number of Epitopes'),
                                                  
                                                  actionButton(
                                                    inputId = 'e_go_3_2',
                                                    label = 'Generate Plot')
                                                  
                                     ),
                                     mainPanel(width = 10,
                                               h4("Correlation between number of epitopes and protein length"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('e_text_density')),
                                                                          style = 'primary')),
                                               fluidRow(column(6, wellPanel(withSpinner(plotlyOutput("e_plot_3_1", height = 500), type = 5))),
                                                        column(6, wellPanel(withSpinner(DT::dataTableOutput("e_table_3"), type = 5))),
                                                        column(2,downloadButton("e_downloadData3", "Download Table"), offset = 6)),
                                               fluidRow(column(12, wellPanel(withSpinner(plotlyOutput("e_plot_3_2", height = 700), type = 5))))
                                     )
                                   )
                          ),
                          tabPanel(title = span('Epitope Viewer', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  
                                                  selectInput(
                                                    inputId = 'e_protein_4',
                                                    label = 'Protein',
                                                    choices =  example_data %>% dplyr::select(ID) %>% unlist() %>% unique(),
                                                    selected = 'sp|P0DTC6|NS6_SARS2'),
                                                  
                                                  numericInput(
                                                    inputId = 'e_cutoff_4',
                                                    label = 'Cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  checkboxGroupInput(
                                                    inputId = 'e_allele_4',
                                                    label = 'MHC alleles',
                                                    choices = colnames(example_data)[c(-1, -2, -3, -4)],
                                                    selected = colnames(example_data)[c(-1, -2, -3, -4)][c(1,2,3,4)]),
                                                  
                                                  radioButtons(
                                                    inputId = 'e_conditional_4',
                                                    label = 'Shared epitopes',
                                                    choices = c('Intersection', 'Union'),
                                                    selected = 'Union'),
                                                  
                                                  actionButton(
                                                    inputId = 'e_go_4',
                                                    label = 'Run Analysis')
                                                  
                                     ),
                                     mainPanel(width = 10,
                                               h4("Amino acid position of the epitopes in proteins"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('e_text_4')),
                                                                          style = 'primary')),
                                               fluidRow(column(12, wellPanel(withSpinner(plotlyOutput("e_plot_4"), type = 5)))),
                                               fluidRow(column(12, wellPanel(withSpinner(DT::dataTableOutput("e_table_4"), type = 5))),
                                                        column(2,downloadButton("e_downloadData4", "Download Table"), offset = 6)))
                                   )
                          ),
                          tabPanel(title = span('Epitope Promiscuity', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  
                                                  'Set Promiscuity Analysis Parameters:',
                                                  br(),
                                                  numericInput(
                                                    inputId = 'e_min_al',
                                                    label = 'Minimum number of alleles',
                                                    min = 0,
                                                    max = 10,
                                                    value = 7),
                                                  
                                                  numericInput(
                                                    inputId = 'e_sb_ct',
                                                    label = 'Strong Binding cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 0.5),
                                                  
                                                  numericInput(
                                                    inputId = 'e_wb_ct',
                                                    label = 'Weak Binding cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  actionButton(
                                                    inputId = 'e_go_5',
                                                    label = 'Run Analysis')
                                                  
                                     ),
                                     mainPanel(width = 10,
                                               h4("Number of alleles and binding rank per epitopes"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('e_text_5')),
                                                                          style = 'primary')),
                                               fluidRow(column(6, wellPanel(withSpinner(plotlyOutput("e_plot_5", height = 500), type = 5))),
                                                        column(6, wellPanel(withSpinner(DT::dataTableOutput("e_table_5"), type = 5)))),
                                               fluidRow(column(2, downloadButton("e_downloadData5", "Download Table"), offset = 6))
                                     )
                                   )),
                          tabPanel(title = span('Epitope Conservation', style ='font-size:120%'), br(),
                                   sidebarLayout(
                                     sidebarPanel(width = 2,
                                                  
                                                  checkboxGroupInput(
                                                    inputId = 'e_protein_7',
                                                    label = 'Proteins',
                                                    choices =  example_data %>% dplyr::select(ID) %>% unlist() %>% unique(),
                                                    selected = c('sp|P0DTC2|SPIKE_SARS2', 'sp|P0DTD8|NS7B_SARS2', 'sp|P0DTD1|R1AB_SARS2')),
                                                  numericInput(
                                                    inputId = 'e_cutoff_7',
                                                    label = 'Cutoff %rank',
                                                    min = 0,
                                                    max = 10,
                                                    value = 2),
                                                  
                                                  checkboxGroupInput(
                                                    inputId = 'e_allele_7',
                                                    label = 'MHC alleles',
                                                    choices = colnames(example_data)[c(-1, -2, -3, -4)],
                                                    selected = colnames(example_data)[c(-1, -2, -3, -4)][c(1,2,3)]),
                                                  
                                                  radioButtons(
                                                    inputId = 'e_conditional_7',
                                                    label = 'Shared epitopes',
                                                    choices = c('Intersection', 'Union'),
                                                    selected = 'Union'),
                                                  radioButtons(
                                                    inputId = 'e_plot_type_7',
                                                    label = 'Plot Type',
                                                    choices = list('Venn Diagram', 'Up-Set'),
                                                    selected = 'Venn Diagram'),
                                                  
                                                  actionButton(
                                                    inputId = 'e_go_7',
                                                    label = 'Run Analysis')
                                                  
                                     ),
                                     mainPanel(width = 10,
                                               h4("Conservation epitopes across variants"),
                                               bsCollapse(bsCollapsePanel('Description (see more)',
                                                                          wellPanel(style= 'overflow-x:auto; height:200px;', uiOutput('e_text_7')),
                                                                          style = 'primary')),
                                               wellPanel(plotlyOutput("e_plot_7", height = "550px"),
                                                         plotOutput("e_plot_7_VD", height = "550px")),
                                               fluidRow(column(12, withSpinner(DT::dataTableOutput("e_table_7"), type = 5)),
                                                        column(2, downloadButton('e_downloadData7', 'Download Table'), offset = 6))
                                     )
                                   )))
                      )),
             tabPanel(title = span("Documentation", style ='font-size:130%'), 
                      sidebarLayout(sidebarPanel(width = 4, h2('Github Documentation'),
                                                 wellPanel(uiOutput('d_text_1'))
                                                 ),
                                    mainPanel(width = 8, 
                                              h2('Supported predictors of T-cell epitopes'),
                                              h3('MHC Class-I epitopes'),
                                              bsCollapse(bsCollapsePanel('NetMHCpan',
                                                                         wellPanel(style= 'overflow-x:auto;',
                                                                                   uiOutput('d_text_netmhcpan')),
                                                                         wellPanel(
                                                                           DT::dataTableOutput('d_table_netmhcpan')),
                                                                         style = 'primary')),
                                              
                                              bsCollapse(bsCollapsePanel('NetMHC',
                                                                         wellPanel(style= 'overflow-x:auto;', uiOutput('d_text_netmhc')),
                                                                         wellPanel(DT::dataTableOutput('d_table_netmhc')),
                                                                         style = 'primary')),
                                              
                                              bsCollapse(bsCollapsePanel('MHCFlurry',
                                                                         wellPanel(style= 'overflow-x:auto', uiOutput('d_text_mhcflurry')),
                                                                         wellPanel(DT::dataTableOutput('d_table_mhcflurry')),
                                                                         style = 'primary')),
                                              
                                              bsCollapse(bsCollapsePanel('IEDB Consensus',
                                                                         wellPanel(style= 'overflow-x:auto;', uiOutput('d_text_iedb')),
                                                                         wellPanel(DT::dataTableOutput('d_table_iedb')),
                                                                         style = 'primary')),
                                              
                                              h3('MHC Class-II epitopes'),
                                              
                                              bsCollapse(bsCollapsePanel('NetMHCIIpan',
                                                                         wellPanel(style= 'overflow-x:auto;', uiOutput('d_text_netmhciipan')),
                                                                         wellPanel(DT::dataTableOutput('d_table_netmhciipan')),
                                                                         style = 'primary')),
                                              
                                              
                                              ))), 
             tabPanel(title = span("Tutorial", style ='font-size:130%'), 
                      sidebarLayout(sidebarPanel(width = 4),
                                    mainPanel(width = 8 ,
                                              bsCollapse(bsCollapsePanel('Epitope-Evaluator Sections',
                                                                         wellPanel(tags$video(id="video1", type = "video/mp4",src = "Sections.mp4", controls = "controls", width = "800px", height = "450px")),
                                                                         style = 'primary')),
                                              bsCollapse(bsCollapsePanel('Upload Files',
                                                                         wellPanel(tags$video(id="video2", type = "video/mp4",src = "UploadFiles.mp4", controls = "controls", width = "800px", height = "450px")),
                                                                         style = 'primary')),
                                              bsCollapse(bsCollapsePanel('Epitope-Distribution',
                                                                         wellPanel(tags$video(id="video3", type = "video/mp4",src = "EpitopeDistribution.mp4", controls = "controls", width = "800px", height = "450px")),
                                                                         style = 'primary')),
                                              bsCollapse(bsCollapsePanel('Epitope-Intersection',
                                                                         wellPanel(tags$video(id="video4", type = "video/mp4",src = "EpitopeIntersection.mp4", controls = "controls", width = "800px", height = "450px")),
                                                                         style = 'primary')),
                                              bsCollapse(bsCollapsePanel('Epitope-Density',
                                                                         wellPanel(tags$video(id="video5", type = "video/mp4",src = "EpitopeDensity.mp4", controls = "controls", width = "800px", height = "450px")),
                                                                         style = 'primary')),
                                              bsCollapse(bsCollapsePanel('Epitope-Viewer',
                                                                         wellPanel(tags$video(id="video6", type = "video/mp4",src = "EpitopeLocation.mp4", controls = "controls", width = "800px", height = "450px")),
                                                                         style = 'primary')),
                                              bsCollapse(bsCollapsePanel('Epitope-Promiscuity',
                                                                         wellPanel(tags$video(id="video7", type = "video/mp4",src = "EpitopePromiscuity.mp4", controls = "controls", width = "800px", height = "450px")),
                                                                         style = 'primary')),
                                              bsCollapse(bsCollapsePanel('Epitope-Conservation',
                                                                         wellPanel(tags$video(id="video8", type = "video/mp4",src = "EpitopeConservation.mp4", controls = "controls", width = "800px", height = "450px")),
                                                                         style = 'primary')),
                                              )))
  )


server <- function(session, input, output){
  
  
  output$text_1_left <-renderUI({
    "Epitope-Evaluator: an interactive web application to study predicted T-cell epitopes"})
  output$d_text_1 <- renderUI({
    github_code <- htmltools::a("GitHub", href="https://github.com/SotoLF/Epitope-Evaluator")
    list(p(HTML(paste("<b>The documentation can be found on</b>", github_code))))
  })
  
  output$d_text_netmhcpan <- renderUI({
    predictor_url <- htmltools::a("NetMHCpan", href="https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1")
    list(p(HTML(paste("<ul> 
                      <li>Paper: NetMHCpan-4.1 and NetMHCIIpan-4.0: Improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data</li>
                      <li>Version: 4.1</li>
                      <li>URL: ", predictor_url,"</li>
                      </ul>"))))
  })
  output$d_table_netmhcpan <- DT::renderDataTable({
    return(data.frame(read.table('example_NetMHCPAN.xls', 
                                 sep = "\t")))
  }, options = list(scrollX = "300px"))
  
  output$d_text_netmhc <- renderUI({
    predictor_url <- htmltools::a("NetMHC", href="https://services.healthtech.dtu.dk/service.php?NetMHC-4.0")
    list(p(HTML(paste("<ul> 
                      <li>Paper: Gapped sequence alignment using artificial neural networks: application to the MHC class I system</li>
                      <li>Version: 4.0</li>
                      <li>URL: ", predictor_url,"</li>
                      </ul>"))))
  })
  output$d_table_netmhc <- DT::renderDataTable({
    return(data.frame(read.table('example_NetMHC.xls', 
                                 sep = "\t",row.names = NULL)))
  }, options = list(scrollX = "300px"))
  
  output$d_text_mhcflurry <- renderUI({
    predictor_url <- htmltools::a("MHCFlurry", href="https://github.com/openvax/mhcflurry")
    list(p(HTML(paste("<ul> 
                      <li>Paper: MHCflurry 2.0: Improved pan-allele prediction of MHC I-presented peptides by incorporating antigen processing</li>
                      <li>Version: 2.0</li>
                      <li>URL: ", predictor_url,"</li>
                      </ul>"))))
  })
  output$d_table_mhcflurry <- DT::renderDataTable({
    return(data.frame(read.table('example_MHCFlurry.txt', 
                                 sep = ",")))
  }, options = list(scrollX = "300px"))
  
  output$d_text_iedb <- renderUI({
    predictor_url <- htmltools::a("IEDB-Consensus", href="http://tools.iedb.org/mhci/")
    list(p(HTML(paste("<ul> 
                      <li>Paper: A consensus epitope prediction approach identifies the breadth of murine T(CD8+)-cell responses to vaccinia virus</li>
                      <li>Version: 2.24</li>
                      <li>URL: ", predictor_url,"</li>
                      </ul>"))))
  })
  output$d_table_iedb <- DT::renderDataTable({
    return(data.frame(read.table('example_IEDB_consensus.txt', 
                                 sep = "\t", row.names = NULL)))
  }, options = list(scrollX = "300px"))
  
  output$d_text_netmhciipan <- renderUI({
    predictor_url <- htmltools::a("NetMHCIIpan", href="https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.0")
    list(p(HTML(paste("<ul> 
                      <li>Paper: NetMHCpan-4.1 and NetMHCIIpan-4.0: Improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data</li>
                      <li>Version: 4.0</li>
                      <li>URL: ", predictor_url,"</li>
                      </ul>"))))
  })
  output$d_table_netmhciipan <- DT::renderDataTable({
    results <- data.frame(read.table('example_NetMHCIIPAN.xls',
                                     sep = ",", header = FALSE))
    
    ## Get alleles names for Class II
    alleles_vector = results[1,] %>% strsplit('\t') %>% unlist()
    alleles_ID = alleles_vector[grepl('D', alleles_vector)]
    
    # Parse Data
    results = results[c(-1,-2),] %>% strsplit('\t')
    results = as.data.frame(do.call(rbind, results))
    
    return(results)
  }, options = list(scrollX = "300px"))
  

  output$text_1 <-renderUI({
    github_code <- htmltools::a("GitHub", href="https://github.com/SotoLF/Epitope-Evaluator")
    
    list(p(HTML(paste("<b> Epitope-Evaluator </b> is a Shiny app aimed to filter and analyze predicted T-cell epitopes. 
    Epitope-Evaluator includes 6 tools: <ol> <li>Epitope Distribution </li> 
    <li>Epitope Intersection </li>
    <li>Epitope Density </li>
    <li>Epitope Viewer </li>
    <li>Epitope Promiscuity </li>
    <li>Epitope Conservation </li>
    </ol> <br>
    Each of the tools shows
    (1) plots that can be downloaded as png files by clicking on the camera icon ('Download plot as png'),
    and (2) tables that can be downloaded by clicking on the 'Download table' button. <br>
    The R scripts of Epitope-Evaluator can be freely downloaded from" , github_code, "and launched locally.
    <br>
    <br>
    <h4><u> <b>How to run Epitope-Evaluator</b> </u></h4>
    It needs two input files: 
    <ol>
    <li>The output prediction file obtained by any T-cell predictor </li>
    <li>the fasta file used during the prediction</li>
    </ol>
    Then, users should select the predictor name. If they are using any predictor not listed, they must upload a tab-delimited output prediction file
    with the following format: The first 4 columns are (i) the peptide sequence, (ii) the position of the peptide within the protein, (iii) the length 
    of the protein, and (iv) the name of the protein. The following columns should be named with the allele containing the %rank score assigned to each
    peptide for each MHC allele. Finally, users should select the score type to filter the epitopes in all the analyses.
    They can choose between the score value and the %rank score. The %rank score is a number between 0 and 100 defined as the rank of the predicted binding score 
    compared to a set of random natural peptides by several predictors of T-cell epitopes. By default, Epitope-Evaluator considers 2 and 10 as the %rank 
    cutoff for MHC Class I and MHC Class II epitopes, respectively.  Once all the input files are uploaded and the predictor and the score type are 
    selected, the users should click on the 'Run analysis' button. We also included a Run example button that loads a file containing the prediction
    of MHC Class I Epitopes from SARS-CoV-2 proteome using NetMHCPan and MHCFlurry.
                       <p>
                       <h4> <b>How to cite Epitope-Evaluator</b></h4>
                       Epitope-Evaluator: an interactive web application to study predicted T-cell epitopes
                       <p>
                       <p>
                       <h4> <b>Contact </h4></b>
                       <ul>
                       <li>Juan Fuxman Bass - fuxman@bu.edu</li>
                       <li>Luis F. Soto - lfs4003@med.cornell.edu</li>
                       </ul>"))))})
  
  output$image_1 <- renderImage({
    
    list(src = 'Github1.png', contentType = 'image/png', res = 300, width = '100%',
         alt = "This is alternate text")    
    
  }, deleteFile = FALSE)
  
  
  cond_pred <- eventReactive(input$go_1, {input$cond_pred})
  method <- eventReactive(input$go_1, {input$method})
  

  
  data <- eventReactive(input$go_1, {
    
    # Check if Output is empty
    validate(need(!is.null(input$input_file), message =  "Need Output Prediction File."))
    
    # Check if fasta is empty
    validate(need(!is.null(input$input_fasta), message = "Need Fasta File."))
    
    # Check if fasta is OK
    validate(need(!is.na(tryCatch(read.fasta(input$input_fasta$datapath), error = function(e) return(NA))),
                  message = 'File submitted is not in Fasta format (See Documentation)'))
    
    # Check if # Proteins in fasta == # Proteins in output prediction
    
    # Check if # Predicion Table is Ok:
    
    
    if (cond_pred() == "NetMHC"){
      validate(need(!is.na(tryCatch(Parse_NetMHC(input$input_file$datapath, input$input_fasta$datapath, method()), error = function(e) return(NA))), 
                    message = 'Prediction File does not correspond to a *NetMHC* Output format (See Documentation)'))
      data <- Parse_NetMHC(input$input_file$datapath, input$input_fasta$datapath, method())
      
    } else if(cond_pred() == "NetMHCpan"){
      validate(need(!is.na(tryCatch(Parse_NetMHCPAN(input$input_file$datapath, input$input_fasta$datapath, method()), error = function(e) return(NA))), 
                    message = 'Prediction File does not correspond to a *NetMHCPan* Output format (See Documentation)'))
      data <- Parse_NetMHCPAN(input$input_file$datapath,input$input_fasta$datapath, method())
      
    } else if(cond_pred() == "NetMHCIIpan"){
      validate(need(!is.na(tryCatch(Parse_NetMHCIIPAN(input$input_file$datapath, input$input_fasta$datapath, method()), error = function(e) return(NA))), 
                    message = 'Prediction File does not correspond to a *NetMHCIIPan* Output format (See Documentation)'))
      data <- Parse_NetMHCIIPAN(input$input_file$datapath,input$input_fasta$datapath, method())
      
      
    } else if(cond_pred() == "MHCFlurry"){
      validate(need(!is.na(tryCatch(Parse_MHCFlurry(input$input_file$datapath, input$input_fasta$datapath, method()), error = function(e) return(NA))), 
                    message = 'Prediction File does not correspond to a *MHCFlurry* Output format (See Documentation)'))
      data <- Parse_MHCFlurry(input$input_file$datapath, input$input_fasta$datapath, method())
      
      
    } else if(cond_pred() == "IEDB Consensus"){
      validate(need(!is.na(tryCatch(Parse_IEDB(input$input_file$datapath, input$input_fasta$datapath, method()), error = function(e) return(NA))), 
                    message = 'Prediction File does not correspond to a *IEDB-Consensus* Output format (See Documentation)'))
      data <- Parse_IEDB(input$input_file$datapath, input$input_fasta$datapath, method())
      
    } else if (cond_pred() == "Other"){
      
      validate(need(!is.na(tryCatch(read.delim(input$input_file$datapath,sep = '\t') %>% as.data.frame(), error = function(e) return(NA))), 
                    message = 'Prediction File does not correspond to the *specified* Output format (See Documentation)'))
      
      data <- read.delim(input$input_file$datapath,
                         sep = "\t") %>% as.data.frame()
      
    }
    data[c(5:ncol(data))] <- sapply(data[c(5:ncol(data))], as.numeric)  
    data[c(5:ncol(data))] <- round(data[c(5:ncol(data))], 2)
    
    return(data)
    
  })
  
  output$contents_1 <- DT::renderDataTable({
    data()
  }, options = list(scrollX = "300px"))
  
  
  #################
  # Parsing data #
  #################
  parsed_data <- reactive({
    
    if (is.null(data())){
      return()}
    
    
    parsed_data <-data()
    
    ## Rename columns
    colnames(parsed_data)[1:4] = c("Peptide", "Pos", "Length", "ID")
    
    ## Updating buttons
    
    updateCheckboxGroupInput(session = session, inputId = "allele", choices = colnames(parsed_data)[c(-1, -2, -3, -4)], selected = colnames(parsed_data)[c(-1, -2, -3, -4)][1])
    updateCheckboxGroupInput(session = session, inputId = "allele_4", choices = colnames(parsed_data)[c(-1, -2, -3, -4)], selected = colnames(parsed_data)[c(-1, -2, -3, -4)][1])
    updateCheckboxGroupInput(session = session, inputId = "allele_6", choices = colnames(parsed_data)[c(-1, -2, -3, -4)], selected = colnames(parsed_data)[c(-1, -2, -3, -4)][1])
    updateCheckboxGroupInput(session = session, inputId = "allele_7", choices = colnames(parsed_data)[c(-1, -2, -3, -4)], selected = colnames(parsed_data)[c(-1, -2, -3, -4)][1])
    
    allele_example <- colnames(parsed_data)[5]
    if(startsWith(toupper(allele_example), "HLA")){
      updateNumericInput(session = session, inputId = "xmin", min = 0, max = 2, value = 0)
      updateNumericInput(session = session, inputId = "xmax", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "cutoff", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "cutoff_4", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "sb_ct", min = 0, max = 2, value = 0.5)
      updateNumericInput(session = session, inputId = "wb_ct", min = 0, max = 10, value = 2)
      updateNumericInput(session = session, inputId = "cutoff_6", min = 0, max = 10, value = 2)
      updateNumericInput(session = session, inputId = "cutoff_7", min = 0, max = 10, value = 2)
    }else{
      updateNumericInput(session = session, inputId = "xmin", min = 0, max = 10, value = 0)
      updateNumericInput(session = session, inputId = "xmax", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutoff", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutoff_4", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "sb_ct", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "wb_ct", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutoff_6", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutoff_7", min = 0, max = 10, value = 10)
    }
    updateNumericInput(session = session, inputId = "n_prots", min = 1, max = parsed_data %>% 
                         dplyr::select(ID) %>% 
                         unlist() %>% 
                         unique() %>% 
                         length(), value = 2)
    updateSelectInput(session = session, inputId = "protein_4", choices = parsed_data %>% dplyr::select(ID) %>% unlist() %>% unique(), selected = parsed_data["ID"][1,])
    updateCheckboxGroupInput(session = session, inputId = "protein_7", choices = parsed_data %>% dplyr::select(ID) %>% unlist() %>% unique(), selected = (parsed_data %>% dplyr::select(ID) %>% unlist() %>% unique())[1:3] )
    
      
    updateNumericInput(session = session,
                       inputId = "min_al",
                       max = length(colnames(parsed_data)[c(-1, -2, -3, -4)]),
                       value = length(colnames(parsed_data)[c(-1, -2, -3, -4)])-1)

    return(parsed_data)
  })
  
  ####################
  # Distribution tab #
  ####################
  
  output$text_distribution <- renderUI({
    list(p(HTML("Epitope Distribution shows a histogram representing the number of epitopes 
    within a percentile rank or binding score range. The histogram can be shown per MHC
    allele or for the union or intersection of different MHC alleles, which can represent
    the number of epitopes that can be recognized by a heterozygote individual or the 
    number of epitopes that are recognized by the alleles present in a population, 
    respectively.  Users can select between plotting a histogram or a cumulative histogram.
    Hovering over any of the histogram bars provides the number of predicted epitopes with
    a percentile rank lower than a selected cutoff. In addition, the tool shows a heatmap
    indicating the number of epitopes predicted to bind to each allele with a percentile
    rank lower than a defined cutoff.

    <br>
    <h4> Parameters </h4>
    <ul>
    <li> <b>MHC alleles:</b> List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = First allele) </li>
    <li> <b>Shared epitopes:</b> When multiple alleles are selected, users must indicate if they are interested in epitopes that bind to all
    the selected MHC alleles (Intersection) or epitopes that bind to at least one of the MHC alleles selected (Union). (Default: Intersection)  </li>
    <li> <b>Min % rank:</b> Minimum percentile rank to filter the epitopes. Peptides with a rank lower
    than the cutoff selected will be ignored. To consider all epitopes, choose 0. If you are interested
    in only weak binder epitopes, you could set 0.5 as cutoff for MHC Class I and 2 for MHC Class II epitopes. (Default = 0)</li>
    <li> <b>Max % rank:</b> Maximum percentile rank to consider a peptide as predicted epitope. Peptides
    with a percentile rank higher than the cutoff selected will be ignored. (Default MHC Class I = 2, MHC Class II = 10).</li>      
    <li> <b>Step:</b> Width of the bins in the histogram/cumulative histogram plot and in the heatmap. (Default = 0.1)</li>
    <li> <b>Plot type:</b> The distribution of epitopes can be shown as a histogram or a cumulative histogram. (Default = histogram)</li>  
    </ul>")))})
  
  reactive_2 <- reactiveValues(doTable = FALSE, doPlot = FALSE)
  
  observeEvent(input$go_2_1, {
    reactive_2$doTable <- input$go_2_1
    reactive_2$doPlot <- input$go_2_1})
  
  allele_react <- eventReactive(input$go_2_1, {input$allele})
  xmin_react <- eventReactive(input$go_2_1, {input$xmin})
  xmax_react <- eventReactive(input$go_2_1, {input$xmax})
  step_react <- eventReactive(input$go_2_1, {input$step})
  conditional_plot_react <- eventReactive(input$go_2_1, {input$conditional_plot})
  conditional_react <- eventReactive(input$go_2_1, {input$conditional})

  # Parsing
  ep_data.score <- reactive({
    
    # Updating click button
    allele <- allele_react()
    conditional <- conditional_react()
    
    ep_data.score <- parsed_data() %>% 
      dplyr::select(Peptide, allele) %>% 
      unique() 
    
    if (length(allele) != 1){
      if ( conditional == "Intersection"){
        ep_data.score$plot_score <- apply(ep_data.score[, allele], 1, max)
      } else if(conditional == "Union"){
        ep_data.score$plot_score <- apply(ep_data.score[, allele], 1, min)
      }
    } else{
      ep_data.score$plot_score <- ep_data.score[,allele]}
    
    return(ep_data.score)
  })
  
  # Parsing
  distribution_table <- reactive({
    
    # Updating click button
    allele <- allele_react()
    xmin <- xmin_react()
    xmax <- xmax_react()
    
    ep_data.score <- ep_data.score()
    tmp_data <- parsed_data()
    
    keep_peptides = tmp_data$Peptide %in% ep_data.score$Peptide[ep_data.score$plot_score > xmin & ep_data.score$plot_score <= xmax]
    distribution_table = tmp_data[keep_peptides, c(1,2,4)]
    distribution_table$Peptide_length = nchar(distribution_table[,1])
    rownames(distribution_table) <- NULL
    return(distribution_table)
  })
  
  # Table
  output$distribution_table <- DT::renderDataTable({
    
    # Updating click button
    
    # Required to update alleles information
    parsed_data()
    
    
    distribution_table()}, rownames = FALSE, options = list(scrollX = "300px"))
  
  # Plot
  output$distribution_plot <- renderPlotly({
    
    # Updating click button
    allele <- allele_react()
    xmin <- xmin_react()
    xmax <- xmax_react()
    step <- step_react()
    conditional_plot <- conditional_plot_react()
    plot_values <- ep_data.score()$plot_score
    
    cumulative = ifelse(conditional_plot == 'Histogram', FALSE, TRUE)
    
    plot_ly(x = plot_values,
            type = "histogram",
            opacity = 0.45,
            hoverinfo = "y",
            marker = list(color = "ligthblue", line = list(cauto = FALSE, width= 1, color = "black", cmid = 2)),
            xbins = list(start = xmin, end = xmax, size = step),
            cumulative = list(enabled = cumulative)) %>% 
      layout(xaxis = list(title = "% Rank Predictor", tickfont = list(size = 20), titlefont = list(size = 20)),
             yaxis = list(title = "Number of epitopes", tickfont = list(size = 20), titlefont = list(size = 20)))
    
    
  })
  
  
  heatmap_data = reactive({
    #if (reactive_2$doTable == FALSE) return()
    allele <- allele_react()
    xmin <- xmin_react()
    xmax <- xmax_react()
    step <- step_react()
    parsed_data <- parsed_data()
    
    cutoffs = seq(xmin, xmax, by = step)
    cutoffs = cutoffs[cutoffs!=0]
    
    heatmap_data = data.frame(alleles = c(),
                              cutoff = c(),
                              number = c())
    
    for(cutoff  in cutoffs){
      tmp_data = apply(parsed_data %>% select(allele), 2, function(x) return(sum(x <= cutoff))) %>% as.data.frame()
      tmp_data$alleles = rownames(tmp_data)
      rownames(tmp_data) = NULL
      tmp_data$cutoff = cutoff
      tmp_data =  tmp_data[, c(2,3,1)]
      colnames(tmp_data) = c('alleles', 'cutoff', 'number')
      heatmap_data = rbind(heatmap_data, tmp_data)
    }
    
    return(heatmap_data)
  })
  
  output$heatmap_plot = renderPlot({
    #if (reactive_2$doPlot == FALSE) return()
    
    allele <- allele_react()
    xmin <- xmin_react()
    xmax <- xmax_react()
    step <- step_react()
    parsed_data <- parsed_data()
    
    cutoffs = seq(xmin, xmax, by = step)
    cutoffs = cutoffs[cutoffs!=0]
    
    heatmap_data = heatmap_data()
    
    p = ggplot(heatmap_data, aes( x= cutoff, y = alleles, fill = number)) + 
      geom_tile(color = 'black')+
      scale_fill_gradient(low = 'white', high = '#0C6394')+
      geom_text(aes(label = number), size = 5)+
      theme_bw(base_size = 22) + 
      scale_x_continuous(breaks = cutoffs) + 
      theme(legend.position = 'none') +
      xlab('') + ylab('')
    return(p)
  })
  
  
  
  # Download button
  output$downloadData2 <- downloadHandler(
    filename = function(){
      paste("Epitope_Distribution.txt")
    },
    content = function(file){
      write.table(distribution_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  
  ####################
  # Intersection tab #
  ####################
  
  # Text
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
     <li> <b>Cutoff %rank:</b> aximum percentile rank to consider a peptide
     as a predicted epitope. Peptides with a percentile rank higher
     than the cutoff selected will be ignored. (Default MHC Class I = 2, MHC Class II = 10)</li>    
     <li> <b>MHC alleles:</b>  List of MHC alleles obtained from
     the input file. Users can choose more than one allele. (Default = First 3 alleles) </li>      
     </ul>")))})
  
  # Reactive
  reactive_6 <- reactiveValues(doPlot = FALSE)
  observeEvent(input$go_6, {reactive_6$doPlot <- input$go_6})
  
  allele_6_react <- eventReactive(input$go_6, {input$allele_6})
  cutoff_6_react <- eventReactive(input$go_6, {input$cutoff_6})
  plot_type_6_react <- eventReactive(input$go_6, {input$plot_type_6})
  
  # Parse
  ep_data6 <- shiny::reactive({
    
    
    # Updating click button
    #if (reactive_6$doPlot == FALSE) return()
    allele_6 <- allele_6_react()
    cutoff_6 <- cutoff_6_react()
    
    ep_data5 <- parsed_data()
    
    new_data <- data.frame(matrix(ncol = length(allele_6)+1, nrow = nrow(ep_data5)))
    colnames(new_data) <- c("Peptide", allele_6)
    new_data$Peptide<- ep_data5$Peptide
    for (allele in allele_6){
      new_data[allele] <- ifelse(ep_data5[allele] <= cutoff_6, allele, NA)
    }
    
    new_data <- new_data[!duplicated(new_data$Peptide),]
    new_data$Alleles = apply(new_data[,-1], 1, function(x) return(paste(x[!is.na(x)], collapse = '-')))
    new_data = new_data %>% filter(Alleles != '')
    new_data = new_data %>% group_by(Alleles) %>% 
      summarize(Sequences = paste(Peptide, collapse = '-')) %>% 
      as.data.frame()
    colnames(new_data)
    new_data %>% select(Alleles) %>% unlist() %>% unique()
    new_data$N_epitopes = apply(new_data %>% select(Sequences), 1, function(x) return(length(unlist(strsplit(x, split = '-')))))
    new_data$N_alleles = apply(new_data %>% select(Alleles), 1, function(x) return(length(unlist(strsplit(x, split = '-')))))
    new_data = new_data[,c(4,1,2,3)]
    return(new_data)
  })
  
  # Table
  output$table_6 <- DT::renderDataTable({
    
    ep_data6() %>% arrange(desc(N_alleles))
    
  }, rownames = FALSE, options = list(scrollX = "300px"))
  
  # Plot (Venn Diagram)
  output$plot_6_VD <-renderPlot({
    
    
    
    ep_data5 <- parsed_data()
    
    
    allele_6 <- allele_6_react()
    cutoff_6 <- cutoff_6_react()
    
    if(length(allele_6)> 6){
      return(NULL)
    }
    
    #print(head(ep_data5))
    #print(allele_6)
    ep_data6 = ep_data5 %>% select(1, allele_6) ## Selecting alleles 
    #print(head(ep_data6))
    #print('middle')
    ep_data6 = melt(ep_data6, id.vars = 'Peptide') %>% as.data.frame()
    #print(head(ep_data6))
    ep_data6 = ep_data6 %>% filter(value <= cutoff_6) ## Filtering based on cutoff
    
    myList = split(ep_data6$Peptide, ep_data6$variable)
    venn <- Venn(myList)
    data <- process_data(venn)
    
    ggplot() +
      geom_sf(aes(fill=count), data = venn_region(data)) +
      geom_sf(size = 1.2, color = "black", data = venn_setedge(data), show.legend = T) +
      geom_sf_text(aes(label = name), data = venn_setlabel(data), size = 6.5) +
      geom_sf_label(aes(label=count), fontface = "bold", data = venn_region(data), label.size = NA, fill = NA, size = 7) +
      theme_void(base_size = 20) + 
      scale_fill_gradient(low = 'white', high = 'red') +
      theme(legend.position = "none")+
      scale_x_continuous(expand = expansion(mult = 0.7))
  })
  
  # Plot (Up-set)
  output$plot_6 <- renderPlotly({
    # Updating click button
    #if (reactive_6$doPlot == FALSE) return()
    
    ep_data5 <- parsed_data()
    
    
    allele_6 <- allele_6_react()
    cutoff_6 <- cutoff_6_react()
    
    # Parsing
    new_data <- data.frame(matrix(ncol = length(allele_6)+1, nrow = nrow(ep_data5)))
    colnames(new_data) <- c("Peptide", allele_6)
    new_data$Peptide<- ep_data5$Peptide
    for (allele in allele_6){
      new_data[allele] <- ep_data5[allele] <= cutoff_6
      new_data[,allele] <- as.integer(new_data[,allele])}
    
    new_data <- new_data[order(new_data[,'Peptide']),]
    new_data <- new_data[!duplicated(new_data$Peptide),]
    i <- 1
    x = vector("list", 2**length(allele_6)-1)
    y = vector("list", 2**length(allele_6)-1)
    
    
    for(j in 1:length(allele_6)){
      combinations <- combn(allele_6,j)
      for (k in 1: ncol(combinations)){
        x[[i]] <- sum(colSums(get_intersect_members(new_data, combinations[,k])))/j
        y[[i]] <- paste(combinations[,k], collapse = ",")
        i = i + 1
      }
    }
    z <- do.call(rbind, Map(data.frame, names = y, values= x))
    z <-z[z$values!=0,]
    nintersections = nrow(z)
    nalleles = length(allele_6)
    set_size <- data.frame(alleles = allele_6, size = colSums(new_data[,2:ncol(new_data)]))
    set_size <- set_size[order(set_size[, 'size']),]
    
    ordered_alleles <- as.character(set_size$alleles)
    dict <- vector("list")
    i = length(ordered_alleles)
    for (allele in ordered_alleles){
      dict[[allele]] <- i
      i = i -1
    }
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
    
    zz <- z[order(z[,"values"], decreasing = TRUE),]
    
    #### Union
    intersect_size_chart<-plot_ly(showlegend = FALSE) %>% add_trace(
      x = 1:nintersections, # Number of bars
      y = zz$values, # Height bar
      type = "bar",
      marker = list(color = "black",
                    hoverinfo = "none")) %>% 
      layout(xaxis = list(tickfont = list(size = 15)),
             yaxis = list(tickfont = list(size = 15))) %>% 
      add_trace(type = "scatter",
                mode = "text",
                x = 1:nintersections,
                y = zz$values+ max(zz$values) * 0.05,
                text = zz$values,
                textfont = list(color = "black"))
    
    ### Grid
    names_intersection <- as.character(zz$names)
    n_grids_rows = 0
    for (name in names_intersection){
      name_splitted <- unlist(strsplit(name, ","))
      n_grids_rows = n_grids_rows + length(name_splitted)
    }
    
    y_axis <- vector("numeric", n_grids_rows)
    x_axis <- vector("numeric", n_grids_rows)
    combo_axis <- vector("numeric", n_grids_rows)
    name_axis <- vector("numeric", n_grids_rows)
    j = 1
    x = 1
    for(name in names_intersection){
      name_splitted <- unlist(strsplit(name, ","))
      for (i in 1:length(name_splitted)){
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
    
    grid <- plot_ly(
      type = "scatter",
      mode = "markers",
      marker = list(color = "lightgrey", size = 8)
    ) %>% add_trace(
      type = "scatter",
      x = rep(1:nintersections,
              length(allele_6)), 
      y = unlist(lapply(1:length(allele_6), function(x)
        rep(x - 0.5, nintersections))),
      hoverinfo = "none"
    ) %>% add_trace(
      data = group_by(lines2, combo),
      type = "scatter",
      mode = "lines+markers",
      x = lines2$x,
      y = lines2$y-0.5,
      hoverinfo = "text",
      text = ~name,
      line = list(color = "black", width = 3),
      marker = list(color = "black",
                    size = 10))%>% 
      layout(xaxis = list(showticklabels = FALSE,
                          showgrid = FALSE,
                          zeroline = FALSE),
             yaxis = list(showticklabels = FALSE,
                          showgrid = TRUE,
                          range = c(0, length(allele_6)),
                          zeroline = FALSE,
                          range = 1:10),
             margin = list(t = 0, b = 40))
    
    intersect_size_chart <-
      intersect_size_chart %>% layout(yaxis = list(title = "Intersections size"))
    
    s1 <-
      subplot(
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        set_size_chart,
        nrows = 2,
        widths = c(0.6, 0.4)
      )
    s2 <-
      subplot(intersect_size_chart,
              grid,
              nrows = 2,
              shareX = TRUE) %>% layout(showlegend = FALSE)
    subplot(s1, s2, widths = c(0.3, 0.7))
  })
  
  # Showing plots
  observeEvent(input$go_6, {
    allele_6_react = input$allele_6
    plot_type_6 = input$plot_type_6
    if (length(allele_6_react)<=6 & plot_type_6 == 'Venn Diagram') {
      hide("plot_6")
      show("plot_6_VD")
    } else {
      hide("plot_6_VD")
      show("plot_6")
    }
  })
  
  # Download button
  output$downloadData6 <- downloadHandler(
    filename = function(){
      paste("Up_Set_epitopes.txt")
    },
    content = function(file){
      write.table(ep_data6(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  
  ###############
  # Density tab #
  ###############
  
  # Text input
  output$text_density <- renderUI({
    list(p(HTML("This tool can be used to determine the set of proteins
    containing a high number of predicted epitopes as a first step to
    finding highly immunogenic proteins. The tool displays a scatter 
    plot of protein length versus the number of epitopes predicted to bind 
    an allele or combination of MHC alleles. Hovering over each point shows 
    the name of the protein, number of epitopes, length of the protein, and 
    the epitope density. Also, selecting any protein from the table will
    highlight the respective point in the scatter plot. The tool also shows 
    the absolute number of epitopes within each protein predicted to bind to
    each MHC allele. This visualization can be displayed as a bar plot for a
    small number of proteins, or as a heatmap for several proteins. In both
    cases, users can modify the plots by changing the fill range (i.e., by
    number or by the density of epitopes) and by arranging the set of proteins.
    <br>
     To obtain the scatter plot, users must indicate the maximum cutoff
     %rank to consider a peptide as an epitope and click on the
     'Run analysis' button. This will produce a scatter plot and 
     a table containing the ID of the proteins, their length in amino acids,
     the number of epitopes, and the epitope density. 
    <br>
    To obtain the heatmap, in addition to the maximum cutoff%rank, users must
    select the plot-type (heatmap or bar plot) and the fill-type (by
    number or density of epitopes). When heatmap is selected, this shows
    the alleles on the x-axis and the proteins on the y-axis. The color 
    intensity indicates the log-10 of the number (or the density) of epitopes.
    When 'barplot' is selected, each protein is represented as a bar while
    alleles are on the x-axis and the number (or density) of epitopes are 
    on the y-axis. Hovering over each bar (or cell in the heatmap) will show 
    the protein ID, the allele, the number of epitopes, the length of the 
    protein, and the density of epitopes.
    <br>
    <h4> Parameters </h4>
    <ul>
     <li> <b>Cutoff % rank:</b> Maximum %rank to consider a peptide as a predicted epitope. (Default MHC Class I = 2, MHC Class II = 10) </li>      
     <li> <b>Color by:</b> Fill heatmap by number or density of epitopes. (Default = Epitopes Number)</li>
     <li> <b>Plot type:</b> Heatmap is recommended for several proteins, while a bar plot is recommended for a few proteins. (Default = Heatmap)</li>
     <li> <b>Sort type:</b> A parameter to sort proteins/alleles in a descendent way or as found in the input file. (Default = Descendent)</li>
     </ul>")))})
  
  # Reactive
  reactive_3 <- reactiveValues(doTable = FALSE, doPlot = FALSE)
  
  observeEvent(input$go_3_1, {
    reactive_3$doTable <- input$go_3_1
    reactive_3$doPlot <- FALSE})
  
  observeEvent(input$go_3_2, {
    reactive_3$doTable <- input$go_3_1
    reactive_3$doPlot <- input$go_3_2})
  
  cutoff_3_react <- eventReactive(input$go_3_1, {input$cutoff})
  cutoff_3_2_react <- eventReactive(input$go_3_2, {input$cutoff})
  
  fill_type_react <- eventReactive(input$go_3_2, {input$fill_type})
  plot_type_react <- eventReactive(input$go_3_2, {input$plot_type})
  
  sort_type_react <- eventReactive(input$go_3_2, {input$sort_type})
  
  # Parsing data
  d_table <- reactive({
    
    # Updating click button
    #if (reactive_3$doTable == FALSE) return()
    cutoff = cutoff_3_react()
    ep_data5 <- parsed_data()
    
    ep_data5$min <- apply(ep_data5[seq(5, ncol(ep_data5))], 1, min)
    ep_data5$min <- as.numeric(ep_data5$min)
    ep_data5 <- ep_data5 %>% dplyr::filter(min <= cutoff) %>% dplyr::group_by(ID,Length) %>% dplyr::summarize(count=n())
    ep_data5$Density <- round(ep_data5$count/ep_data5$Length,2)
    return(ep_data5)
  })
  
  # Table
  output$table_3 <- DT::renderDataTable({
    
    # Updating click button
    #if (reactive_3$doTable == FALSE) return()
    ep_data3 <- d_table()
    DT::datatable(ep_data3, rownames = FALSE, options= list(scrollX = "300px"))})
  
  # Plot
  output$plot_3_1 <- renderPlotly({
    
    # Updating click button
    #if (reactive_3$doTable == FALSE) return()
    ep_data5 <- d_table()
    
    s = input$table_3_rows_selected
    y_min_value <- min(ep_data5$count)
    y_max_value <- max(ep_data5$count)
    x_min_value <- min(ep_data5$Length)
    x_max_value <- max(ep_data5$Length)
    ids = as.data.frame(ep_data5)[s, ]
    
    plot_ly(type = "scatter") %>% 
      add_trace(data = ep_data5[ep_data5$ID  %in% ids$ID,],
                x = ~as.factor(Length), y = ~as.factor(count),
                marker = list(size = 10,
                              color = "red",
                              line = list(color = "black",
                                          width = 2)),
                showlegend = FALSE,
                
                hoverinfo = "text",
                
                text = paste("Protein: ", ep_data5[ep_data5$ID  %in%  ids$ID,]$ID,
                             "<br>",
                             "N° epitopes: ", ep_data5[ep_data5$ID  %in%  ids$ID,]$count,
                             "<br>",
                             "Length: ", ep_data5[ep_data5$ID  %in%  ids$ID,]$Length, "aminoacids",
                             "<br>",
                             "Density: ", round(ep_data5[ep_data5$ID %in% ids$ID,]$Density, 2))) %>% 
      
      add_trace(data = ep_data5[!(ep_data5$ID  %in%  ids$ID),],
                x = ~Length, y = ~count,
                marker = list(size = 10,
                              color = "rgba(255, 182, 193, .1)",
                              line = list(color = "black",
                                          width = 2)),
                showlegend = FALSE,
                
                hoverinfo = "text",
                
                text = paste("Protein: ", ep_data5[!(ep_data5$ID  %in%  ids$ID),]$ID,
                             "<br>",
                             "N° epitopes: ", ep_data5[!(ep_data5$ID  %in%  ids$ID),]$count,
                             "<br>",
                             "Length: ", ep_data5[!(ep_data5$ID  %in%  ids$ID),]$Length, "aminoacids",
                             "<br>",
                             "Density: ", round(ep_data5[!(ep_data5$ID %in% ids$ID),]$Density, 2))) %>%
      
      layout(title = "Correlation N° Epitopes and Proteins Length",
             plot_bgcolor='#a9c8f5',
             yaxis = list(title = "N of epitopes", zerolinecolor = 'black',range = c(0, y_max_value),
                          zerolinewidth = 2,gridcolor = 'black',
                          tickfont = list(size = 20),  titlefont = list(size = 20)),
             xaxis = list(title = "Protein Length",zerolinecolor = 'black', range = c(0, x_max_value),
                          zerolinewidth = 2, gridcolor = 'black',tickfont = list(size = 20),  titlefont = list(size = 20)))
    
  })
  
  # Plot
  output$plot_3_2 <- renderPlotly({
    
    # Updating click button
    #if (reactive_3$doPlot == FALSE) return()
    cutoff = cutoff_3_2_react()
    table <- parsed_data()
    fill_type <- fill_type_react()
    sort_type <- sort_type_react()
    plot_type <- plot_type_react()
    table_melt <- melt(table, id.vars = c("Pos", "Length", "ID", "Peptide"))
    protein_length <- table_melt %>% select(ID, Length) %>% unique()
    
    y_ids <- table$ID %>% unique() %>% unlist() %>% as.character()
    x_ids <- colnames(table)[5:ncol(table)]
    
    table_summarized <- table_melt %>% group_by(ID, variable) %>% filter(value <= cutoff) %>% 
      summarise(count = n()) %>% as.data.frame()
    index = match(table_summarized$ID, protein_length$ID)
    
    
    if (fill_type == "Number of Epitopes"){
      table_summarized$Length <- protein_length$Length[index] 
      table_summarized$Density <- table_summarized$count/table_summarized$Length
      table_summarized$type_fill <- log10(table_summarized$count)
    }
    else{
      table_summarized$Length <- protein_length$Length[index] 
      table_summarized$Density <- table_summarized$count/table_summarized$Length
      table_summarized$type_fill <- log10(table_summarized$count/table_summarized$Length)
    }
    
    
    if (sort_type != "Default"){
      y_ids <- table_summarized %>% group_by(ID) %>% summarise(accumulated = sum(10**type_fill)) %>% 
        arrange(desc(accumulated)) %>% select(ID) %>% unlist()
      x_ids <- table_summarized %>% group_by(variable) %>% summarise(accumulated = sum(10**type_fill)) %>% 
        arrange(desc(accumulated)) %>% select(variable) %>% unlist() %>% as.character()
    }
    
    if (plot_type == "Bar plot"){
      if(fill_type == "Number of Epitopes"){
        p <-ggplot(table_summarized)+ 
          geom_bar(stat = "identity" ,
                   color ="black",
                   position = "dodge",
                   width = 0.85,
                   aes(x = variable,
                       y = count,
                       fill = factor(ID, levels = y_ids),
                       text = paste("Protein:", ID,
                                    "<br>",
                                    "Allele:", variable,
                                    "<br>",
                                    "N_epitopes: ", count,
                                    "<br>",
                                    "Length: ", Length,
                                    "<br>",
                                    "Density: ", round(Density, 2))))+
          theme_classic(base_size = 20)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15),
                axis.text.y = element_text(size = 15),
                legend.title = element_blank())+
          xlab("") + ylab("") + 
          xlim(x_ids)
      }else{
        p <-ggplot(table_summarized)+ 
          geom_bar(stat = "identity" ,
                   color ="black",
                   position = "dodge",
                   width = 0.85,
                   aes(x = variable,
                       y = Density, 
                       text = paste("Protein:", ID,
                                    "<br>",
                                    "Allele:", variable,
                                    "<br>",
                                    "N_epitopes: ", count,
                                    "<br>",
                                    "Length: ", Length,
                                    "<br>",
                                    "Density: ", round(Density, 2)),
                       fill = factor(ID, levels = y_ids)))+
          theme_classic(base_size = 20)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15),
                axis.text.y = element_text(size = 15),
                legend.title = element_blank())+
          xlab("") + ylab("") + 
          xlim(x_ids)
      }
      
    }else{
      min_value = min(table_summarized$type_fill)-0.1
      max_value = round(max(table_summarized$type_fill),2) +0.1
      round_dec = ifelse(fill_type == 'Number of Epitopes', 0, 2)
      name_legend = ifelse(fill_type == 'Number of Epitopes', 'Number', 'Density')
      p <- ggplot(table_summarized) + 
        geom_tile(aes(x = variable, 
                      y = ID, 
                      fill = type_fill, 
                      text = paste("Protein:", ID,
                                   "<br>",
                                   "Allele:", variable,
                                   "<br>",
                                   "N_epitopes: ", count,
                                   "<br>",
                                   "Length: ", Length,
                                   "<br>",
                                   "Density: ", round(Density, 2))))+
        theme_classic(base_size = 20)+
        scale_fill_gradientn(name = name_legend,
                             colors = c("white","#1E74A5"),
                             limits = c(min_value, max_value),
                             breaks = c(min_value, mean(c(min_value,max_value)), max_value),
                             labels = c(round(10**min_value,round_dec), round(10**mean(c(min_value,max_value)),round_dec), round(10**max_value,round_dec)),
                             guide = guide_colorbar(barwidth = 1.5,barheight = 8, nbin=50, ticks=FALSE, frame.colour="black"))+
        xlab("") + ylab("") +ylim(y_ids) + xlim(x_ids) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=15),
              axis.text.y = element_text(size = 15))
    }
    
    ggplotly(p, tooltip = c("text"))
  })
  
  # Download button
  output$downloadData3 <- downloadHandler(
    filename = function(){
      paste("Density_Protein.txt")
    },
    content = function(file){
      write.table(d_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })

  ################
  # Viewer tab #
  ################
  
  # Text
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
     </ul>")))})
  
  # Reactive
  reactive_4 <- reactiveValues(doPlot = FALSE)
  observeEvent(input$go_4, {reactive_4$doPlot <- input$go_4})
  allele_4_react <- eventReactive(input$go_4, {input$allele_4})
  cutoff_4_react <- eventReactive(input$go_4, {input$cutoff_4})
  protein_4_react <- eventReactive(input$go_4, {input$protein_4})
  conditional_4_react <- eventReactive(input$go_4, {input$conditional_4})
  
  # Parse
  location_table <- reactive({
    #if (reactive_4$doPlot == FALSE) return()
    
    ep_data5 <- parsed_data()
    
    
    allele_4 <- allele_4_react()
    cutoff_4 <- cutoff_4_react()
    protein_4 <- protein_4_react()
    conditional_4 <- conditional_4_react()
    
    # Parse
    peptide_counts = parsed_data() %>% 
      filter(ID == protein_4) %>% 
      select(Peptide, Pos, Length, allele_4)
    
    if (length(allele_4) != 1){
      if ( conditional_4 == "Intersection"){
        peptide_counts$score <- apply(peptide_counts[, allele_4], 1, max)
      } else if(conditional_4 == "Union"){
        peptide_counts$score <- apply(peptide_counts[, allele_4], 1, min)
      }
    }
    else{
      peptide_counts$score <- peptide_counts[,allele_4]
    }
    peptide_counts <- peptide_counts %>% filter(score <= cutoff_4) %>% 
      dplyr::select(Peptide, Pos, Length, allele_4) %>% 
      melt(id.vars = c("Peptide", "Pos", "Length")) %>% 
      group_by(Peptide) %>% 
      filter(value <= 2) %>% 
      mutate(count = sum(value <= cutoff_4),
             alleles = paste(variable, collapse = " - ")) %>% 
      select(Peptide, count, alleles) %>% as.data.frame()
    
    if(nrow(peptide_counts)==0) return()
    
    ep_data2 = parsed_data() %>%
      filter(ID == protein_4) %>% 
      select(Peptide, Pos, Length)
    
    index = match(ep_data2$Peptide, peptide_counts$Peptide)
    
    ep_data2$Value = peptide_counts$count[index]
    ep_data2$alleles = peptide_counts$alleles[index]
    ep_data2 <- ep_data2 %>% filter(Value > 0)
    
    length_protein <- ep_data2 %>% dplyr::select(Length) %>% unlist() %>% unique()
    local_table <-data.frame(type = c("Protein", rep("Epitope", nrow(ep_data2))),
                             start = c(1, ep_data2$Pos+1),
                             end = c(length_protein+1, ep_data2$Pos+1+nchar(ep_data2$Peptide)),
                             sequence = c(protein_4, ep_data2$Peptide),
                             nivel = c(1, rep(0, nrow(ep_data2))),
                             id = seq(0, nrow(ep_data2)),
                             alleles = c("-", ep_data2$alleles),
                             count = c("NA", ep_data2$Value))
    
    nivel = 0
    while(nivel %in% local_table$nivel){
      subdata <- local_table[local_table$type == "Epitope" & local_table$nivel==nivel,]
      changed = TRUE
      while (changed){
        number_rows = nrow(subdata)
        if(number_rows == 1){
          changed = FALSE
          nivel = nivel - 1
          break
        }
        conteo = 0
        for (i in 1:number_rows){
          # Static row
          static_row <- subdata[i,]
          # Static vector
          v_s <- seq(static_row$start, static_row$end)
          new_data <- subdata[-i,]
          # Compare all vs first vector
          conteo = conteo + 1
          for (row in 1:nrow(new_data)){
            row_data <- new_data[row,]
            v_x <- seq(row_data$start, row_data$end)
            id <- row_data$id
            if (length(intersect(v_s, v_x)) != 0){
              local_table[local_table$id == id,]$nivel <- row_data$nivel - 1
              subdata[subdata$id == id, ]$nivel <- row_data$nivel - 1
            }
          }
          subdata <- subdata[subdata$nivel ==nivel, ]
          final_number_rows = nrow(subdata)
          if (number_rows != final_number_rows){
            break
          }
          if(number_rows == final_number_rows & conteo == number_rows){
            changed = FALSE
          }
        }
      }
      nivel = nivel - 1
    }
    
    return(local_table)
  })
  
  # Table
  output$table_4 <- DT::renderDataTable({
    
    # Updating click button
    #if (reactive_4$doPlot == FALSE) return()
    local_table <- location_table()
    if (is.null(local_table)) return(data.frame(sequence = character(),
                                                start = integer(),
                                                end = integer(),
                                                alleles = character(),
                                                count = integer()))
    local_table <- local_table[-1,c(4,2,3,7,8)]
    DT::datatable(local_table, rownames = FALSE,options= list(scrollX = "300px"))})
  
  # Plot
  output$plot_4 <- renderPlotly({
    #if (reactive_4$doPlot == FALSE) return()
    
    ep_data5 <- parsed_data()
    protein_4 <- protein_4_react()
    
    local_table <- location_table()
    if (is.null(local_table)) return()
    
    local_table$count <- as.numeric(local_table$count) 
    
    if (nrow(local_table) == 1){
      p <- ggplot()+
        geom_rect(data = local_table, aes(xmin = start,
                                          xmax = end,
                                          ymin = nivel -0.3,
                                          ymax = nivel + 0.3,
                                          fill = count,
                                          sequence = "-",
                                          text = paste("Class:", type,"<br>","Seq:", sequence)),
                  alpha = 0.9,
                  size = 0.3,
                  colour = "black")+ 
        geom_text(data = local_table[local_table$type == 'Protein',],
                  aes(x = start/2 + end/2, y = (nivel-0.3)/2 + (nivel+0.3)/2), label = 'protein_4')+
        scale_fill_gradientn(colors = c("yellow", "red"),
                             limits = c(0, max(local_table$count)),
                             na.value = "lightblue")+
        theme_bw(base_size = 12) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")+
        scale_y_continuous(expand = c(0,0), limits = c(-5, 1.5))
      ggplotly(p, tooltip = c("text")) %>% layout(dragmode = "zoom",showlegend = FALSE)
    }
    else{
      
      steps = c(5, 10, 25, 50, 100, 150, 2000, 250, 500, 1000, 1500, 2500, 5000)
      n_splits = 15
      max_x = local_table[local_table$type == "Protein",]$end
      step_x_1 = max_x / n_splits
      step_x = steps[which(abs(steps - step_x_1) == min(abs(steps - step_x_1)))]
      
      p <- ggplot()+
        geom_rect(data = local_table, aes(xmin = start,
                                          xmax = end,
                                          ymin = nivel-0.3,
                                          ymax = nivel+0.3,
                                          fill = count,
                                          text = paste("Seq: ", sequence, 
                                                       "<br>",
                                                       "Position: ", start, "-", end-1,
                                                       "<br>",
                                                       "Alleles: ", alleles)),
                  alpha = 0.9,
                  size = 0.3,
                  colour = "black")+
        geom_text(data = local_table[local_table$type == 'Protein',],
                  aes(x = start/2 + end/2, y = (nivel-0.3)/2 + (nivel+0.3)/2), label = protein_4)+
        scale_fill_gradientn(colors = c("yellow", "red"),
                             limits = c(0, max(local_table$count)),
                             na.value = "lightblue")+
        scale_y_continuous(expand = c(0,0), limits = c(min(local_table$nivel)-1, 1.5)) +
        scale_x_continuous(limits = c(-1, max_x) ,breaks = seq(0, max_x, by = step_x))+
        theme_bw(base_size = 12) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")
      ggplotly(p, tooltip = c("text")) %>% layout(dragmode = "zoom",showlegend = FALSE)
    }
    
  })
  
  # Download button
  output$downloadData4 <- downloadHandler(
    filename = function(){
      paste("Epitopes_Plot_Location.txt")
    },
    content = function(file){
      write.table(location_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  ###################
  # Promiscuity tab #
  ###################
  
  # Text
  output$text_5 <- renderUI({
    list(p(HTML("Epitope Promiscuity shows epitopes predicted to bind a minimum
    number of MHC alleles. The tool depicts a heatmap where alleles are on the x-axis
    and epitopes are on the y-axis. Users should set the strong and weak binding 
    cutoff %rank to filter the strong binder and weak binder epitopes. Strong binder epitopes (SB) are the 
    epitopes with a %rank less or equal than the strong binding cutoff %rank. 
    Weak binder epitopes (WB) are the epitopes with a %rank greater than the strong binding
    cutoff %rank but less or equal than the weak binding cutoff %rank. In the heatmap, 
    SB and WB epitopes are colored in red and orange, respectively. Users also 
    need to select a minimum number of MHC alleles to filter and show only epitopes 
    that are predicted to bind to at least this number of alleles. This tool also shows
    a table containing the amino acid sequence, the position within the protein, the ID 
    of the protein, and the number of alleles to which they are predicted to bind.
    <br>
    <h4> Parameters </h4>
    <ul>
     <li> <b>Minimum number of alleles:</b> This indicates the minimum number of MHC alleles that an epitope needs to bind to be shown in the heatmap.</li>    
     <li> <b>Strong Binding Cutoff % rank:</b> Epitopes with a %rank lower or equal than this cutoff are considered as strong binder epitopes (SB). </li>      
     <li> <b>Weak Binding Cutoff %rank:</b> Epitopes with a %rank lower or equal than this cutoff but greater than Strong Binding Cutoff are considered weak binder epitopes (WB). </li>
     </ul>")))})
  
  # Reactive
  reactive_5 <- reactiveValues(doPlot = FALSE)
  observeEvent(input$go_5, {reactive_5$doPlot <- input$go_5})
  wb_ct_react <- eventReactive(input$go_5, {input$wb_ct})
  sb_ct_react <- eventReactive(input$go_5, {input$sb_ct})
  min_al_react <- eventReactive(input$go_5, {input$min_al})
  
  # Parsing
  parsed_data5 <- shiny::reactive({
    
    
    # Update click button
    #if (reactive_5$doPlot == FALSE) return()
    ep_data5 <- parsed_data()
    wb_ct = wb_ct_react()
    sb_ct = sb_ct_react()
    min_al = min_al_react()
    
    # Parse
    ep_data5$prom <- rowSums(ep_data5[,seq(5, ncol(ep_data5))] < wb_ct)
    ep_data5 <- ep_data5 %>% filter(prom >= min_al) 
    ids <- ep_data5 %>% select(Peptide) %>% unlist()
    ep_data5 <- ep_data5[order(ep_data5[,"prom"], decreasing = TRUE),]
    ep_data5 <- ep_data5 %>% select("Peptide", "Pos", "Length", "ID", "prom")
    
  })
  
  # Table
  output$table_5 <- DT::renderDataTable({
    
    # Update click button
    #if (reactive_5$doPlot == FALSE) return()
    parsed_data5()
    
  }, rownames = FALSE, options = list(scrollX = "300px"))
  
  # Plot
  output$plot_5 <- renderPlotly({
    
    
    if(nrow(parsed_data5()) == 0 ) return()
    
    wb_ct = wb_ct_react()
    sb_ct = sb_ct_react()
    min_al = min_al_react()
    
    # Parse
    ep_data5 <- parsed_data()
    ep_data5$prom <- rowSums(ep_data5[,seq(5, ncol(ep_data5))] < wb_ct)
    ep_data5 <- ep_data5 %>% filter(prom >= min_al) 
    ids <- ep_data5 %>% select(Peptide) %>% unlist()
    ep_data5 <- ep_data5[order(ep_data5[,"prom"], decreasing = TRUE),]
    asd<- factor(levels= ep_data5$Peptide,
                 labels = ep_data5$Peptide)
    mdata <- parsed_data() %>% filter( Peptide %in% ids)
    mdata <- mdata[order(mdata[,'Peptide']),]
    mdata <- mdata[!duplicated(mdata$Peptide),]
    mdata <- melt(mdata, id=c("Peptide", "Pos", "Length", "ID"))
    
    mdata$value2[mdata$value <= wb_ct & mdata$value > sb_ct] = "WB"
    mdata$value2[mdata$value <= sb_ct] = "SB"
    mdata$value2[mdata$value > wb_ct] = "NA"
    mdata$value2 <- as.factor(mdata$value2)
    
    p <-ggplot(mdata) + geom_tile(aes(x=variable, y=Peptide, fill=value2, cutoff = value), width = 0.9, height=0.9, color ="black", size= 0.3) +
      scale_fill_manual(values = c("NA" = "white", "SB"= "red", "WB" = "orange"))+
      ylim(levels(asd))+ xlab("Alleles") + ylab("")+
      theme_bw(base_size = 18)+
      xlab('')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.title =element_blank())+ xlab("Alleles")
    ggplotly(p, tooltip = c("cutoff", "y"))
    
  })
  
  # Download button
  output$downloadData5 <- downloadHandler(
    filename = function(){
      paste("Epitope_promiscuity.txt")
    },
    content = function(file){
      write.table(parsed_data5(), file, sep = "\t",row.names = FALSE, quote = FALSE)
      
    })
  
  
  ####################
  # Conservation tab #
  ####################
  
  # Text
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

     </ul>")))})
  
  # Reactive
  reactive_7 <- reactiveValues(doPlot = FALSE)
  observeEvent(input$go_7, {reactive_7$doPlot <- input$go_7})
  
  allele_7_react <- eventReactive(input$go_7, {input$allele_7})
  cutoff_7_react <- eventReactive(input$go_7, {input$cutoff_7})
  protein_7_react <- eventReactive(input$go_7, {input$protein_7})
  plot_type_7_react <- eventReactive(input$go_7, {input$plot_type_7})
  conditional_7_react <- eventReactive(input$go_7, {input$conditional_7})
  
  
  table_epitopes <- reactive({
    parsed = parsed_data() 
    allele_7 <- allele_7_react()
    cutoff_7 <- cutoff_7_react()
    protein_7 <- protein_7_react()
    plot_type_7 <- plot_type_7_react()
    conditional_7 <- conditional_7_react()
    
    parsed = parsed %>% filter(ID %in% protein_7)
    parsed_1 = parsed %>% select(Peptide, allele_7) %>% unique()
    
    if (length(allele_7) != 1){
      if (conditional_7 == 'Intersection'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, max)
      } else if (conditional_7 == 'Union'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, min) 
      }
    }else{
      parsed_1$plot_score <- parsed_1[, allele_7]
    }
    
    candidate_peptides = parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
    parsed_2 = parsed %>% select(Peptide, ID, allele_7) %>% filter(Peptide %in% candidate_peptides)
    peptide_info = dcast(parsed_2[, c(1,2)], formula = Peptide~ID) %>% as.data.frame() 
    peptide_info$Proteins = apply(peptide_info[,-1], 1, function(x) return(paste(x[!is.na(x)], collapse = '-')))
    peptide_info = peptide_info %>% 
      group_by(Proteins) %>%
      summarize(sequences = paste(Peptide, collapse = '-')) %>%
      as.data.frame()
    peptide_info$N_sequences = apply(peptide_info %>% select(sequences), 1, function(x) return(length(unlist(strsplit(x, split = '-')))))
    return(peptide_info)
  })
  
  # Table
  output$table_7 <- DT::renderDataTable({
    # Updating Click button
    #if (reactive_7$doPlot == FALSE) return()
    
    table_epitopes()
  }, rownames = FALSE, options = list(scrollX = '300px'))
  
  output$plot_7_VD <- renderPlot({
    
    #if (reactive_7$doPlot == FALSE) return()
    
    allele_7 <- allele_7_react()
    cutoff_7 <- cutoff_7_react()
    protein_7 <- protein_7_react()
    plot_type_7 <- plot_type_7_react()
    conditional_7 <- conditional_7_react()
    
    if(length(protein_7)> 6){
      return(NULL)
    }
    
    parsed = parsed_data() %>% filter(ID %in% protein_7)
    parsed_1 = parsed %>% select(Peptide, allele_7) %>% unique()
    
    if (length(allele_7) != 1){
      if (conditional_7 == 'Intersection'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, max)
      } else if (conditional_7 == 'Union'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, min) 
      }
    }else{
      parsed_1$plot_score <- parsed_1[, allele_7]
    }
    
    candidate_peptides = parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
    parsed_2 = parsed %>% select(Peptide, ID, allele_7) %>% filter(Peptide %in% candidate_peptides)
    
    myList = split(parsed_2$Peptide, parsed_2$ID)
    venn <- Venn(myList)
    data <- process_data(venn)
    
    return(ggplot() +
             geom_sf(aes(fill=count), data = venn_region(data)) +
             geom_sf(size = 1.2, color = "black", data = venn_setedge(data), show.legend = T) +
             geom_sf_text(aes(label = name), data = venn_setlabel(data), nudge_y = 0.1, size = 6.5) +
             geom_sf_label(aes(label=count), fontface = "bold", data = venn_region(data), label.size = NA, fill = NA, size = 7) +
             theme_void(base_size = 20) + scale_fill_gradient(low = 'white', high = 'red') +
             theme(legend.position = "none") +
             scale_x_continuous(expand = expansion(mult = 0.5)))
    
    
  })
  
  output$plot_7 <- renderPlotly({
    #if (reactive_7$doPlot == FALSE) return()
    
    allele_7 <- allele_7_react()
    cutoff_7 <- cutoff_7_react()
    protein_7 <- protein_7_react()
    plot_type_7 <- plot_type_7_react()
    conditional_7 <- conditional_7_react()
    
    parsed = parsed_data() %>% filter(ID %in% protein_7)
    parsed_1 = parsed %>% select(Peptide, allele_7) %>% unique()
    
    if (length(allele_7) != 1){
      if (conditional_7 == 'Intersection'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, max)
      } else if (conditional_7 == 'Union'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, min) 
      }
    }else{
      parsed_1$plot_score <- parsed_1[, allele_7]
    }
    
    candidate_peptides = parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
    parsed_2 = parsed %>% select(Peptide, ID, allele_7) %>% filter(Peptide %in% candidate_peptides)
    
    new_data2 = parsed_2 %>% dcast(Peptide~ID)
    new_data2[,-1] = apply(new_data2[,-1], 2, function(x) return(as.integer(ifelse(is.na(x), 0, 1)))) %>% as.data.frame()
    colnames(new_data2) = gsub("\\|", "_", colnames(new_data2))
    
    protein = colnames(new_data2[,-1])
    
    i <- 1
    x = vector("list", 2**length(protein)-1)
    y = vector("list", 2**length(protein)-1)
    
    
    for(j in 1:length(protein)){
      combinations <- combn(protein,j)
      for (k in 1: ncol(combinations)){
        x[[i]] <- sum(colSums(get_intersect_members(new_data2, combinations[,k])))/j
        y[[i]] <- paste(combinations[,k], collapse = "-")
        i = i + 1
      }
    }
    z <- do.call(rbind, Map(data.frame, names = y, values= x))
    z <-z[z$values!=0,]
    
    nintersections = nrow(z)
    nalleles = length(protein)
    set_size <- data.frame(alleles = protein, size = colSums(new_data2[,2:ncol(new_data2)]))
    set_size <- set_size[order(set_size[, 'size']),]
    
    ordered_alleles <- as.character(set_size$alleles)
    dict <- vector("list")
    i = length(ordered_alleles)
    for (allele in ordered_alleles){
      dict[[allele]] <- i
      i = i -1
    }
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
    
    zz <- z[order(z[,"values"], decreasing = TRUE),]
    
    #### Union
    intersect_size_chart<-plot_ly(showlegend = FALSE) %>% add_trace(
      x = 1:nintersections, # Number of bars
      y = zz$values, # Height bar
      type = "bar",
      marker = list(color = "black",
                    hoverinfo = "none"))  %>% 
      layout(xaxis = list(tickfont = list(size = 15)),
             yaxis = list(tickfont = list(size = 15))) %>% 
      add_trace(type = "scatter",
                mode = "text",
                x = 1:nintersections,
                y = zz$values+ max(zz$values) * 0.05,
                text = zz$values,
                textfont = list(color = "black"))
    
    ### Grid
    names_intersection <- as.character(zz$names)
    n_grids_rows = 0
    for (name in names_intersection){
      name_splitted <- unlist(strsplit(name, ","))
      n_grids_rows = n_grids_rows + length(name_splitted)
    }
    
    y_axis <- vector("numeric", n_grids_rows)
    x_axis <- vector("numeric", n_grids_rows)
    combo_axis <- vector("numeric", n_grids_rows)
    name_axis <- vector("numeric", n_grids_rows)
    j = 1
    x = 1
    for(name in names_intersection){
      name_splitted <- unlist(strsplit(name, "-"))
      for (i in 1:length(name_splitted)){
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
    
    grid <- plot_ly(
      type = "scatter",
      mode = "markers",
      marker = list(color = "lightgrey", size = 8)
    ) %>% add_trace(
      type = "scatter",
      x = rep(1:nintersections,
              length(protein)), 
      y = unlist(lapply(1:length(protein), function(x)
        rep(x - 0.5, nintersections))),
      hoverinfo = "none"
    ) %>% add_trace(
      data = group_by(lines2, combo),
      type = "scatter",
      mode = "lines+markers",
      x = lines2$x,
      y = lines2$y-0.5,
      hoverinfo = "text",
      text = ~name,
      line = list(color = "black", width = 3),
      marker = list(color = "black",
                    size = 10))%>% 
      layout(xaxis = list(showticklabels = FALSE,
                          showgrid = FALSE,
                          zeroline = FALSE),
             yaxis = list(showticklabels = FALSE,
                          showgrid = TRUE,
                          range = c(0, length(protein)),
                          zeroline = FALSE,
                          range = 1:10),
             margin = list(t = 0, b = 40))
    
    intersect_size_chart <-
      intersect_size_chart %>% layout(yaxis = list(title = "Intersections size"))
    
    s1 <-
      subplot(
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        set_size_chart,
        nrows = 2,
        widths = c(0.6, 0.4)
      )
    s2 <-
      subplot(intersect_size_chart,
              grid,
              nrows = 2,
              shareX = TRUE) %>% layout(showlegend = FALSE)
    subplot(s1, s2, widths = c(0.3, 0.7))
    
    
    
  })
  
  # Showing plots
  observeEvent(input$plot_type_7, {
    req(input$plot_type_7)
    if (input$plot_type_7 == "Venn Diagram") {
      hide("plot_7")
      show("plot_7_VD")
    } else {
      hide("plot_7_VD")
      show("plot_7")
    }
  })
  
  # Showing plots
  observeEvent(input$go_7, {
    protein_7_react = input$protein_7
    plot_type_7 = input$plot_type_7
    if (length(protein_7_react)<=6 & plot_type_7 == 'Venn Diagram') {
      hide("plot_7")
      show("plot_7_VD")
    } else {
      hide("plot_7_VD")
      show("plot_7")
    }
  })
  
  
  
  
  ###### Example
  
  # Transform Table
  e_data <- reactive({
    example_data = Parse_NetMHCPAN('example.xls','example.fasta', 'Rank')
    example_data[c(5:ncol(example_data))] <- sapply(example_data[c(5:ncol(example_data))], as.numeric)  
    example_data[c(5:ncol(example_data))] <- round(example_data[c(5:ncol(example_data))], 2)
    return(example_data)
    })
  
  output$e_contents_1 <- DT::renderDataTable({
    e_data()}, options = list(scrollX = '300px'))
  
  # Download button
  output$e_downloadData_example <- downloadHandler(
    filename = function(){
      paste("Ouput_example.txt")
    },
    content = function(file){
      write.table(e_data(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  
  output$e_text_1 <- renderUI({
    fasta_file <- htmltools::a("(Download file)", href="https://github.com/SotoLF/Epitope-Evaluator/blob/main/Examples/example.fasta")
    prediction_file <- htmltools::a("(Download file)", href="https://github.com/SotoLF/Epitope-Evaluator/blob/main/Examples/example.xls")
    list(p(HTML(paste("The example data shows the <b> MHC class I epitopes</b> from the SARS-CoV-2 proteome using the NetMHCPan4.1.
    <br>The input files (below) and the parsed table (right) can be downloaded:
    <ol> 
    <li>Fasta File:", fasta_file, "</li>
    <li>Output Prediction File:", prediction_file, "</li>
    </ol>
                      
    The parameters for the prediction file were the following:
                      <ul>
                      <li>PEPTIDE LENGTH: 9mer peptides</li>
                      <li>ALLELES: HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01 </li>
                      <li>Strong Binder: 0.5 </li>
                      <li>Weak Binder: 2 </li>
                      <li>Save prediction to XLS file: Yes </li>
                      </ul>
                      "))))
    
  })
  
  e_parsed_data <- reactive({  
    
    parsed_data <-e_data()
    
    ## Rename columns
    colnames(parsed_data)[1:4] = c("Peptide", "Pos", "Length", "ID")
    return(parsed_data)
  })
  
  ####################
  # Distribution tab #
  ####################
  
  output$e_text_distribution <- renderUI({
    list(p(HTML("Epitope Distribution shows a histogram representing the number of epitopes 
    within a percentile rank or binding score range. The histogram can be shown per MHC
    allele or for the union or intersection of different MHC alleles, which can represent
    the number of epitopes that can be recognized by a heterozygote individual or the 
    number of epitopes that are recognized by the alleles present in a population, 
    respectively.  Users can select between plotting a histogram or a cumulative histogram.
    Hovering over any of the histogram bars provides the number of predicted epitopes with
    a percentile rank lower than a selected cutoff. In addition, the tool shows a heatmap
    indicating the number of epitopes predicted to bind to each allele with a percentile
    rank lower than a defined cutoff.

    <br>
    <h4> Parameters </h4>
    <ul>
    <li> <b>MHC alleles:</b> List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = First allele) </li>
    <li> <b>Shared epitopes:</b> When multiple alleles are selected, users must indicate if they are interested in epitopes that bind to all
    the selected MHC alleles (Intersection) or epitopes that bind to at least one of the MHC alleles selected (Union). (Default: Intersection)  </li>
    <li> <b>Min % rank:</b> Minimum percentile rank to filter the epitopes. Peptides with a rank lower
    than the cutoff selected will be ignored. To consider all epitopes, choose 0. If you are interested
    in only weak binder epitopes, you could set 0.5 as cutoff for MHC Class I and 2 for MHC Class II epitopes. (Default = 0)</li>
    <li> <b>Max % rank:</b> Maximum percentile rank to consider a peptide as predicted epitope. Peptides
    with a percentile rank higher than the cutoff selected will be ignored. (Default MHC Class I = 2, MHC Class II = 10).</li>      
    <li> <b>Step:</b> Width of the bins in the histogram/cumulative histogram plot and in the heatmap. (Default = 0.1)</li>
    <li> <b>Plot type:</b> The distribution of epitopes can be shown as a histogram or a cumulative histogram. (Default = histogram)</li>  
    </ul>")))})

  e_reactive_2 <- reactiveValues(doTable = FALSE, doPlot = FALSE)
  
  observeEvent(input$e_go_2_1, {
    e_reactive_2$doTable <- input$e_go_2_1
    e_reactive_2$doPlot <- input$e_go_2_1})
  
  # Parsing
  e_ep_data.score <- reactive({
    
    # Updating click button
    #allele <- e_allele_react()
    allele <- input$e_allele
    #conditional <- e_conditional_react()
    conditional <- input$e_conditional
    
    ep_data.score <- e_parsed_data() %>% 
      dplyr::select(Peptide, allele) %>% 
      unique() 
    
    if (length(allele) != 1){
      if ( conditional == "Intersection"){
        ep_data.score$plot_score <- apply(ep_data.score[, allele], 1, max)
      } else if(conditional == "Union"){
        ep_data.score$plot_score <- apply(ep_data.score[, allele], 1, min)
      }
    } else{
      ep_data.score$plot_score <- ep_data.score[,allele]}
    
    return(ep_data.score)
  })
  
  
  # Parsing
  e_distribution_table <- eventReactive(c(input$e_go_2_1), {
    
    # Updating click button
    allele <- input$e_allele
    xmin <- input$e_xmin
    xmax <- input$e_xmax
    
    ep_data.score <- e_ep_data.score()
    tmp_data <- e_parsed_data()
    
    keep_peptides = tmp_data$Peptide %in% ep_data.score$Peptide[ep_data.score$plot_score > xmin & ep_data.score$plot_score <= xmax]
    distribution_table = tmp_data[keep_peptides, c(1,2,4)]
    distribution_table$Peptide_length = nchar(distribution_table[,1])
    rownames(distribution_table) <- NULL
    return(distribution_table)
  })
  
  # Table
  output$e_distribution_table <- DT::renderDataTable({
    
    # Updating click button
    
    # Required to update alleles information
    e_parsed_data()
    
    
    e_distribution_table()}, rownames = FALSE, options = list(scrollX = "300px"))
  
  
  e_plotly <- eventReactive(c(input$e_go_2_1), {
    # Updating click button
    allele <- input$e_allele
    xmin <- input$e_xmin
    xmax <- input$e_xmax
    step <- input$e_step
    conditional_plot <- input$e_conditional_plot
    plot_values <- e_ep_data.score()$plot_score
    
    cumulative = ifelse(conditional_plot == 'Histogram', FALSE, TRUE)
    
    plot_ly(x = plot_values,
            type = "histogram",
            opacity = 0.45,
            hoverinfo = "y",
            marker = list(color = "ligthblue", line = list(cauto = FALSE, width= 1, color = "black", cmid = 2)),
            xbins = list(start = xmin, end = xmax, size = step),
            cumulative = list(enabled = cumulative)) %>% 
      layout(xaxis = list(title = "% Rank Predictor", tickfont = list(size = 20), titlefont = list(size = 20)),
             yaxis = list(title = "Number of epitopes", tickfont = list(size = 20), titlefont = list(size = 20)))
  })
  
  # Plot
  output$e_distribution_plot <- renderPlotly({
    
    e_plotly()
    
    
  })
  
  
  e_heatmap_data = eventReactive(c(input$e_go_2_1), {
    
    allele <- input$e_allele
    xmin <- input$e_xmin
    xmax <- input$e_xmax
    step <- input$e_step
    e_parsed_data <- e_parsed_data()
    
    cutoffs = seq(xmin, xmax, by = step)
    cutoffs = cutoffs[cutoffs!=0]
    
    heatmap_data = data.frame(alleles = c(),
                              cutoff = c(),
                              number = c())
    
    for(cutoff  in cutoffs){
      tmp_data = apply(e_parsed_data %>% select(allele), 2, function(x) return(sum(x <= cutoff))) %>% as.data.frame()
      tmp_data$alleles = rownames(tmp_data)
      rownames(tmp_data) = NULL
      tmp_data$cutoff = cutoff
      tmp_data =  tmp_data[, c(2,3,1)]
      colnames(tmp_data) = c('alleles', 'cutoff', 'number')
      heatmap_data = rbind(heatmap_data, tmp_data)
    }
    
    return(heatmap_data)
  })
  
  output$e_heatmap_plot = renderPlot({
    xmin <- input$e_xmin
    xmax <- input$e_xmax
    step <- input$e_step
    e_parsed_data <- e_parsed_data()
    
    cutoffs = seq(xmin, xmax, by = step)
    cutoffs = cutoffs[cutoffs!=0]
    
    heatmap_data = e_heatmap_data()
    
    p = ggplot(heatmap_data, aes( x= cutoff, y = alleles, fill = number)) + 
      geom_tile(color = 'black')+
      scale_fill_gradient(low = 'white', high = '#0C6394')+
      geom_text(aes(label = number), size = 5)+
      theme_bw(base_size = 22) + 
      scale_x_continuous(breaks = cutoffs) + 
      theme(legend.position = 'none') +
      xlab('') + ylab('')
    return(p)
  })
  
  
  
  # Download button
  output$e_downloadData2 <- downloadHandler(
    filename = function(){
      paste("Epitope_Distribution_example.txt")
    },
    content = function(file){
      write.table(e_distribution_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  
  
  
  ####################
  # Intersection tab #
  ####################
  
  # Text
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
     <li> <b>Cutoff %rank:</b> aximum percentile rank to consider a peptide
     as a predicted epitope. Peptides with a percentile rank higher
     than the cutoff selected will be ignored. (Default MHC Class I = 2, MHC Class II = 10)</li>    
     <li> <b>MHC alleles:</b>  List of MHC alleles obtained from
     the input file. Users can choose more than one allele. (Default = First 3 alleles) </li>      
     </ul>")))})
  
  # Reactive
  e_reactive_6 <- reactiveValues(doPlot = FALSE)
  observeEvent(input$e_go_6, {e_reactive_6$doPlot <- input$e_go_6})
  
  e_allele_6_react <- eventReactive(input$e_go_6, {input$e_allele_6})
  e_cutoff_6_react <- eventReactive(input$e_go_6, {input$e_cutoff_6})
  e_plot_type_6_react <- eventReactive(input$e_go_6, {input$e_plot_type_6})
  
  
  # Parse
  e_ep_data6 <- eventReactive(c(input$e_go_6), {
    
    
    # Updating click button
    #if (reactive_6$doPlot == FALSE) return()
    allele_6 <- input$e_allele_6
    cutoff_6 <- input$e_cutoff_6
    
    ep_data5 <- e_parsed_data()
    
    new_data <- data.frame(matrix(ncol = length(allele_6)+1, nrow = nrow(ep_data5)))
    colnames(new_data) <- c("Peptide", allele_6)
    new_data$Peptide<- ep_data5$Peptide
    for (allele in allele_6){
      new_data[allele] <- ifelse(ep_data5[allele] <= cutoff_6, allele, NA)
    }
    
    new_data <- new_data[!duplicated(new_data$Peptide),]
    new_data$Alleles = apply(new_data[,-1], 1, function(x) return(paste(x[!is.na(x)], collapse = '-')))
    new_data = new_data %>% filter(Alleles != '')
    new_data = new_data %>% group_by(Alleles) %>% 
      summarize(Sequences = paste(Peptide, collapse = '-')) %>% 
      as.data.frame()
    colnames(new_data)
    new_data %>% select(Alleles) %>% unlist() %>% unique()
    new_data$N_epitopes = apply(new_data %>% select(Sequences), 1, function(x) return(length(unlist(strsplit(x, split = '-')))))
    new_data$N_alleles = apply(new_data %>% select(Alleles), 1, function(x) return(length(unlist(strsplit(x, split = '-')))))
    new_data = new_data[,c(4,1,2,3)]
    return(new_data)
  })
  
  # Table
  output$e_table_6 <- DT::renderDataTable({
    
    e_ep_data6() %>% arrange(desc(N_alleles))
    
  }, rownames = FALSE, options = list(scrollX = "300px"))
  
  # Plot (Venn Diagram)
  
  t_plot_6_VD <- eventReactive(c(input$e_go_6),{
    ep_data5 <- e_parsed_data()
    
    
    allele_6 <- input$e_allele_6
    cutoff_6 <- input$e_cutoff_6
    
    if(length(allele_6)> 6){
      return(NULL)
    }
    
    ep_data6 = ep_data5 %>% select(1, allele_6) ## Selecting alleles 
    ep_data6 = melt(ep_data6, id.vars = 'Peptide') %>% as.data.frame()
    ep_data6 = ep_data6 %>% filter(value <= cutoff_6) ## Filtering based on cutoff
    
    myList = split(ep_data6$Peptide, ep_data6$variable)
    venn <- Venn(myList)
    data <- process_data(venn)
    
    ggplot() +
      geom_sf(aes(fill=count), data = venn_region(data)) +
      geom_sf(size = 1.2, color = "black", data = venn_setedge(data), show.legend = T) +
      geom_sf_text(aes(label = name), data = venn_setlabel(data), size = 6.5) +
      geom_sf_label(aes(label=count), fontface = "bold", data = venn_region(data), label.size = NA, fill = NA, size = 7) +
      theme_void(base_size = 20) + 
      scale_fill_gradient(low = 'white', high = 'red') +
      theme(legend.position = "none")+
      scale_x_continuous(expand = expansion(mult = 0.7))
  })
  
  output$e_plot_6_VD <-renderPlot({
    
    t_plot_6_VD()
    
  })
  
  # Plot (Up-set)
  
  t_plot_6 <- eventReactive(c(input$e_go_6), {
    
    ep_data5 <- e_parsed_data()
    
    
    allele_6 <- input$e_allele_6
    cutoff_6 <- input$e_cutoff_6
    
    # Parsing
    new_data <- data.frame(matrix(ncol = length(allele_6)+1, nrow = nrow(ep_data5)))
    colnames(new_data) <- c("Peptide", allele_6)
    new_data$Peptide<- ep_data5$Peptide
    for (allele in allele_6){
      new_data[allele] <- ep_data5[allele] <= cutoff_6
      new_data[,allele] <- as.integer(new_data[,allele])}
    
    new_data <- new_data[order(new_data[,'Peptide']),]
    new_data <- new_data[!duplicated(new_data$Peptide),]
    i <- 1
    x = vector("list", 2**length(allele_6)-1)
    y = vector("list", 2**length(allele_6)-1)
    
    
    for(j in 1:length(allele_6)){
      combinations <- combn(allele_6,j)
      for (k in 1: ncol(combinations)){
        x[[i]] <- sum(colSums(get_intersect_members(new_data, combinations[,k])))/j
        y[[i]] <- paste(combinations[,k], collapse = ",")
        i = i + 1
      }
    }
    z <- do.call(rbind, Map(data.frame, names = y, values= x))
    z <-z[z$values!=0,]
    nintersections = nrow(z)
    nalleles = length(allele_6)
    set_size <- data.frame(alleles = allele_6, size = colSums(new_data[,2:ncol(new_data)]))
    set_size <- set_size[order(set_size[, 'size']),]
    
    ordered_alleles <- as.character(set_size$alleles)
    dict <- vector("list")
    i = length(ordered_alleles)
    for (allele in ordered_alleles){
      dict[[allele]] <- i
      i = i -1
    }
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
    
    zz <- z[order(z[,"values"], decreasing = TRUE),]
    
    #### Union
    intersect_size_chart<-plot_ly(showlegend = FALSE) %>% add_trace(
      x = 1:nintersections, # Number of bars
      y = zz$values, # Height bar
      type = "bar",
      marker = list(color = "black",
                    hoverinfo = "none")) %>% 
      layout(xaxis = list(tickfont = list(size = 15)),
             yaxis = list(tickfont = list(size = 15))) %>% 
      add_trace(type = "scatter",
                mode = "text",
                x = 1:nintersections,
                y = zz$values+ max(zz$values) * 0.05,
                text = zz$values,
                textfont = list(color = "black"))
    
    ### Grid
    names_intersection <- as.character(zz$names)
    n_grids_rows = 0
    for (name in names_intersection){
      name_splitted <- unlist(strsplit(name, ","))
      n_grids_rows = n_grids_rows + length(name_splitted)
    }
    
    y_axis <- vector("numeric", n_grids_rows)
    x_axis <- vector("numeric", n_grids_rows)
    combo_axis <- vector("numeric", n_grids_rows)
    name_axis <- vector("numeric", n_grids_rows)
    j = 1
    x = 1
    for(name in names_intersection){
      name_splitted <- unlist(strsplit(name, ","))
      for (i in 1:length(name_splitted)){
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
    
    grid <- plot_ly(
      type = "scatter",
      mode = "markers",
      marker = list(color = "lightgrey", size = 8)
    ) %>% add_trace(
      type = "scatter",
      x = rep(1:nintersections,
              length(allele_6)), 
      y = unlist(lapply(1:length(allele_6), function(x)
        rep(x - 0.5, nintersections))),
      hoverinfo = "none"
    ) %>% add_trace(
      data = group_by(lines2, combo),
      type = "scatter",
      mode = "lines+markers",
      x = lines2$x,
      y = lines2$y-0.5,
      hoverinfo = "text",
      text = ~name,
      line = list(color = "black", width = 3),
      marker = list(color = "black",
                    size = 10))%>% 
      layout(xaxis = list(showticklabels = FALSE,
                          showgrid = FALSE,
                          zeroline = FALSE),
             yaxis = list(showticklabels = FALSE,
                          showgrid = TRUE,
                          range = c(0, length(allele_6)),
                          zeroline = FALSE,
                          range = 1:10),
             margin = list(t = 0, b = 40))
    
    intersect_size_chart <-
      intersect_size_chart %>% layout(yaxis = list(title = "Intersections size"))
    
    s1 <-
      subplot(
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        set_size_chart,
        nrows = 2,
        widths = c(0.6, 0.4)
      )
    s2 <-
      subplot(intersect_size_chart,
              grid,
              nrows = 2,
              shareX = TRUE) %>% layout(showlegend = FALSE)
    subplot(s1, s2, widths = c(0.3, 0.7))
    
  })
  output$e_plot_6 <- renderPlotly({

   t_plot_6()
  })
  
  # Showing plots
  observeEvent(input$e_go_6, {
    e_allele_6_react = input$e_allele_6
    e_plot_type_6 = input$e_plot_type_6
    if (length(e_allele_6_react)<=6 & e_plot_type_6 == 'Venn Diagram') {
      hide("e_plot_6")
      show("e_plot_6_VD")
    } else {
      hide("e_plot_6_VD")
      show("e_plot_6")
    }
  })
  
  # Download button
  output$e_downloadData6 <- downloadHandler(
    filename = function(){
      paste("Up_Set_epitopes_example.txt")
    },
    content = function(file){
      write.table(e_ep_data6(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  
  
  # Density tab #
  ###############
  
  # Text input
  output$e_text_density <- renderUI({
    list(p(HTML("This tool can be used to determine the set of proteins
    containing a high number of predicted epitopes as a first step to
    finding highly immunogenic proteins. The tool displays a scatter 
    plot of protein length versus the number of epitopes predicted to bind 
    an allele or combination of MHC alleles. Hovering over each point shows 
    the name of the protein, number of epitopes, length of the protein, and 
    the epitope density. Also, selecting any protein from the table will
    highlight the respective point in the scatter plot. The tool also shows 
    the absolute number of epitopes within each protein predicted to bind to
    each MHC allele. This visualization can be displayed as a bar plot for a
    small number of proteins, or as a heatmap for several proteins. In both
    cases, users can modify the plots by changing the fill range (i.e., by
    number or by the density of epitopes) and by arranging the set of proteins.
    <br>
     To obtain the scatter plot, users must indicate the maximum cutoff
     %rank to consider a peptide as an epitope and click on the
     'Run analysis' button. This will produce a scatter plot and 
     a table containing the ID of the proteins, their length in amino acids,
     the number of epitopes, and the epitope density. 
    <br>
    To obtain the heatmap, in addition to the maximum cutoff%rank, users must
    select the plot-type (heatmap or bar plot) and the fill-type (by
    number or density of epitopes). When heatmap is selected, this shows
    the alleles on the x-axis and the proteins on the y-axis. The color 
    intensity indicates the log-10 of the number (or the density) of epitopes.
    When 'barplot' is selected, each protein is represented as a bar while
    alleles are on the x-axis and the number (or density) of epitopes are 
    on the y-axis. Hovering over each bar (or cell in the heatmap) will show 
    the protein ID, the allele, the number of epitopes, the length of the 
    protein, and the density of epitopes.
    <br>
    <h4> Parameters </h4>
    <ul>
     <li> <b>Cutoff % rank:</b> Maximum %rank to consider a peptide as a predicted epitope. (Default MHC Class I = 2, MHC Class II = 10) </li>      
     <li> <b>Color by:</b> Fill heatmap by number or density of epitopes. (Default = Epitopes Number)</li>
     <li> <b>Plot type:</b> Heatmap is recommended for several proteins, while a bar plot is recommended for a few proteins. (Default = Heatmap)</li>
     <li> <b>Sort type:</b> A parameter to sort proteins/alleles in a descendent way or as found in the input file. (Default = Descendent)</li>
     </ul>")))})
  
  # Reactive
  e_reactive_3 <- reactiveValues(doTable = FALSE, doPlot = FALSE)
  
  observeEvent(input$e_go_3_1, {
    e_reactive_3$e_doTable <- input$e_go_3_1
    e_reactive_3$e_doPlot <- FALSE})
  
  observeEvent(input$e_go_3_2, {
    e_reactive_3$doTable <- input$e_go_3_1
    e_reactive_3$doPlot <- input$e_go_3_2})

  
  # Parsing data
  e_d_table <- eventReactive(c(input$e_go_3_1), {
    
    # Updating click button
    #if (reactive_3$doTable == FALSE) return()
    cutoff = input$e_cutoff
    ep_data5 <- e_parsed_data()
    
    ep_data5$min <- apply(ep_data5[seq(5, ncol(ep_data5))], 1, min)
    ep_data5$min <- as.numeric(ep_data5$min)
    ep_data5 <- ep_data5 %>% dplyr::filter(min <= cutoff) %>% dplyr::group_by(ID,Length) %>% dplyr::summarize(count=n())
    ep_data5$Density <- round(ep_data5$count/ep_data5$Length,2)
    return(ep_data5)
  })
  
  # Table
  output$e_table_3 <- DT::renderDataTable({
    
    # Updating click button
    #if (reactive_3$doTable == FALSE) return()
    ep_data3 <- e_d_table()
    DT::datatable(ep_data3, rownames = FALSE, options= list(scrollX = "300px"))})
  
  # Plot
  
  output$e_plot_3_1 <- renderPlotly({
    
    # Updating click button
    #if (reactive_3$doTable == FALSE) return()
    ep_data5 <- e_d_table()
    
    s = input$e_table_3_rows_selected
    y_min_value <- min(ep_data5$count)
    y_max_value <- max(ep_data5$count)
    x_min_value <- min(ep_data5$Length)
    x_max_value <- max(ep_data5$Length)
    ids = as.data.frame(ep_data5)[s, ]
    
    plot_ly(type = "scatter") %>% 
      add_trace(data = ep_data5[ep_data5$ID  %in% ids$ID,],
                x = ~as.factor(Length), y = ~as.factor(count),
                marker = list(size = 10,
                              color = "red",
                              line = list(color = "black",
                                          width = 2)),
                showlegend = FALSE,
                
                hoverinfo = "text",
                
                text = paste("Protein: ", ep_data5[ep_data5$ID  %in%  ids$ID,]$ID,
                             "<br>",
                             "N° epitopes: ", ep_data5[ep_data5$ID  %in%  ids$ID,]$count,
                             "<br>",
                             "Length: ", ep_data5[ep_data5$ID  %in%  ids$ID,]$Length, "aminoacids",
                             "<br>",
                             "Density: ", round(ep_data5[ep_data5$ID %in% ids$ID,]$Density, 2))) %>% 
      
      add_trace(data = ep_data5[!(ep_data5$ID  %in%  ids$ID),],
                x = ~Length, y = ~count,
                marker = list(size = 10,
                              color = "rgba(255, 182, 193, .1)",
                              line = list(color = "black",
                                          width = 2)),
                showlegend = FALSE,
                
                hoverinfo = "text",
                
                text = paste("Protein: ", ep_data5[!(ep_data5$ID  %in%  ids$ID),]$ID,
                             "<br>",
                             "N° epitopes: ", ep_data5[!(ep_data5$ID  %in%  ids$ID),]$count,
                             "<br>",
                             "Length: ", ep_data5[!(ep_data5$ID  %in%  ids$ID),]$Length, "aminoacids",
                             "<br>",
                             "Density: ", round(ep_data5[!(ep_data5$ID %in% ids$ID),]$Density, 2))) %>%
      
      layout(title = "Correlation N° Epitopes and Proteins Length",
             plot_bgcolor='#a9c8f5',
             yaxis = list(title = "N of epitopes", zerolinecolor = 'black',range = c(0, y_max_value),
                          zerolinewidth = 2,gridcolor = 'black',
                          tickfont = list(size = 20),  titlefont = list(size = 20)),
             xaxis = list(title = "Protein Length",zerolinecolor = 'black', range = c(0, x_max_value),
                          zerolinewidth = 2, gridcolor = 'black',tickfont = list(size = 20),  titlefont = list(size = 20)))
    
  })
  
  # Plot
  t_plot_3_2 <- eventReactive(c(input$e_go_3_2),{
    
    # Updating click button
    #if (reactive_3$doPlot == FALSE) return()
    cutoff = input$e_cutoff
    table <- e_parsed_data()
    fill_type <- input$e_fill_type
    sort_type <- input$e_sort_type
    plot_type<- input$e_plot_type
    table_melt <- melt(table, id.vars = c("Pos", "Length", "ID", "Peptide"))
    protein_length <- table_melt %>% select(ID, Length) %>% unique()
    
    y_ids <- table$ID %>% unique() %>% unlist() %>% as.character()
    x_ids <- colnames(table)[5:ncol(table)]
    
    table_summarized <- table_melt %>% group_by(ID, variable) %>% filter(value <= cutoff) %>% 
      summarise(count = n()) %>% as.data.frame()
    index = match(table_summarized$ID, protein_length$ID)
    
    if (fill_type == "Number of Epitopes"){
      table_summarized$Length <- protein_length$Length[index] 
      table_summarized$Density <- table_summarized$count/table_summarized$Length
      table_summarized$type_fill <- log10(table_summarized$count)
    }
    else{
      table_summarized$Length <- protein_length$Length[index] 
      table_summarized$Density <- table_summarized$count/table_summarized$Length
      table_summarized$type_fill <- log10(table_summarized$count/table_summarized$Length)
    }
    
    
    if (sort_type != "Default"){
      y_ids <- table_summarized %>% group_by(ID) %>% summarise(accumulated = sum(10**type_fill)) %>% 
        arrange(desc(accumulated)) %>% select(ID) %>% unlist()
      x_ids <- table_summarized %>% group_by(variable) %>% summarise(accumulated = sum(10**type_fill)) %>% 
        arrange(desc(accumulated)) %>% select(variable) %>% unlist() %>% as.character()
    }
    
    if (plot_type == "Bar plot"){
      if(fill_type == "Number of Epitopes"){
        p <-ggplot(table_summarized)+ 
          geom_bar(stat = "identity" ,
                   color ="black",
                   position = "dodge",
                   width = 0.85,
                   aes(x = variable,
                       y = count,
                       fill = factor(ID, levels = y_ids),
                       text = paste("Protein:", ID,
                                    "<br>",
                                    "Allele:", variable,
                                    "<br>",
                                    "N_epitopes: ", count,
                                    "<br>",
                                    "Length: ", Length,
                                    "<br>",
                                    "Density: ", round(Density, 2))))+
          theme_classic(base_size = 20)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15),
                axis.text.y = element_text(size = 15),
                legend.title = element_blank())+
          xlab("") + ylab("") + 
          xlim(x_ids)
      }else{
        p <-ggplot(table_summarized)+ 
          geom_bar(stat = "identity" ,
                   color ="black",
                   position = "dodge",
                   width = 0.85,
                   aes(x = variable,
                       y = Density, 
                       text = paste("Protein:", ID,
                                    "<br>",
                                    "Allele:", variable,
                                    "<br>",
                                    "N_epitopes: ", count,
                                    "<br>",
                                    "Length: ", Length,
                                    "<br>",
                                    "Density: ", round(Density, 2)),
                       fill = factor(ID, levels = y_ids)))+
          theme_classic(base_size = 20)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15),
                axis.text.y = element_text(size = 15),
                legend.title = element_blank())+
          xlab("") + ylab("") + 
          xlim(x_ids)
      }
      
    }else{
      min_value = min(table_summarized$type_fill)-0.1
      max_value = round(max(table_summarized$type_fill),2) +0.1
      round_dec = ifelse(fill_type == 'Number of Epitopes', 0, 2)
      name_legend = ifelse(fill_type == 'Number of Epitopes', 'Number', 'Density')
      p <- ggplot(table_summarized) + 
        geom_tile(aes(x = variable, 
                      y = ID, 
                      fill = type_fill, 
                      text = paste("Protein:", ID,
                                   "<br>",
                                   "Allele:", variable,
                                   "<br>",
                                   "N_epitopes: ", count,
                                   "<br>",
                                   "Length: ", Length,
                                   "<br>",
                                   "Density: ", round(Density, 2))))+
        theme_classic(base_size = 20)+
        scale_fill_gradientn(name = name_legend,
                             colors = c("white","#1E74A5"),
                             limits = c(min_value, max_value),
                             breaks = c(min_value, mean(c(min_value,max_value)), max_value),
                             labels = c(round(10**min_value,round_dec), round(10**mean(c(min_value,max_value)),round_dec), round(10**max_value,round_dec)),
                             guide = guide_colorbar(barwidth = 1.5,barheight = 8, nbin=50, ticks=FALSE, frame.colour="black"))+
        xlab("") + ylab("") +ylim(y_ids) + xlim(x_ids) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=15),
            axis.text.y = element_text(size = 15))
    }
    
    ggplotly(p, tooltip = c("text"))
  })
  output$e_plot_3_2 <- renderPlotly({
    t_plot_3_2()
    
  })
  
  # Download button
  output$e_downloadData3 <- downloadHandler(
    filename = function(){
      paste("Density_Protein_example.txt")
    },
    content = function(file){
      write.table(e_d_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  ################

  # Location tab #
  ################
  
  # Text
  output$e_text_4 <- renderUI({
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
     </ul>")))})
  
  # Reactive
  e_reactive_4 <- reactiveValues(doPlot = FALSE)
  observeEvent(input$e_go_4, {e_reactive_4$doPlot <- input$e_go_4})
  #e_allele_4_react <- eventReactive(input$e_go_4, {input$e_allele_4})
  #e_cutoff_4_react <- eventReactive(input$e_go_4, {input$e_cutoff_4})
  #e_protein_4_react <- eventReactive(input$e_go_4, {input$e_protein_4})
  #e_conditional_4_react <- eventReactive(input$e_go_4, {input$e_conditional_4})
  
  # Parse
  e_location_table <- eventReactive(c(input$e_go_4),{
    #if (reactive_4$doPlot == FALSE) return()
    
    ep_data5 <- e_parsed_data()
    
    
    allele_4 <- input$e_allele_4
    cutoff_4 <- input$e_cutoff_4
    protein_4 <- input$e_protein_4
    conditional_4 <- input$e_conditional_4
    
    # Parse
    peptide_counts = e_parsed_data() %>% 
      filter(ID == protein_4) %>% 
      select(Peptide, Pos, Length, allele_4)
    
    if (length(allele_4) != 1){
      if ( conditional_4 == "Intersection"){
        peptide_counts$score <- apply(peptide_counts[, allele_4], 1, max)
      } else if(conditional_4 == "Union"){
        peptide_counts$score <- apply(peptide_counts[, allele_4], 1, min)
      }
    }
    else{
      peptide_counts$score <- peptide_counts[,allele_4]
    }
    peptide_counts <- peptide_counts %>% filter(score <= cutoff_4) %>% 
      dplyr::select(Peptide, Pos, Length, allele_4) %>% 
      melt(id.vars = c("Peptide", "Pos", "Length")) %>% 
      group_by(Peptide) %>% 
      filter(value <= 2) %>% 
      mutate(count = sum(value <= cutoff_4),
             alleles = paste(variable, collapse = " - ")) %>% 
      select(Peptide, count, alleles) %>% as.data.frame()
    
    if(nrow(peptide_counts)==0) return()
    
    ep_data2 = e_parsed_data() %>%
      filter(ID == protein_4) %>% 
      select(Peptide, Pos, Length)
    
    index = match(ep_data2$Peptide, peptide_counts$Peptide)
    
    ep_data2$Value = peptide_counts$count[index]
    ep_data2$alleles = peptide_counts$alleles[index]
    ep_data2 <- ep_data2 %>% filter(Value > 0)
    
    length_protein <- ep_data2 %>% dplyr::select(Length) %>% unlist() %>% unique()
    local_table <-data.frame(type = c("Protein", rep("Epitope", nrow(ep_data2))),
                             start = c(1, ep_data2$Pos+1),
                             end = c(length_protein+1, ep_data2$Pos+1+nchar(ep_data2$Peptide)),
                             sequence = c(protein_4, ep_data2$Peptide),
                             nivel = c(1, rep(0, nrow(ep_data2))),
                             id = seq(0, nrow(ep_data2)),
                             alleles = c("-", ep_data2$alleles),
                             count = c("NA", ep_data2$Value))
    
    nivel = 0
    while(nivel %in% local_table$nivel){
      subdata <- local_table[local_table$type == "Epitope" & local_table$nivel==nivel,]
      changed = TRUE
      while (changed){
        number_rows = nrow(subdata)
        if(number_rows == 1){
          changed = FALSE
          nivel = nivel - 1
          break
        }
        conteo = 0
        for (i in 1:number_rows){
          # Static row
          static_row <- subdata[i,]
          # Static vector
          v_s <- seq(static_row$start, static_row$end)
          new_data <- subdata[-i,]
          # Compare all vs first vector
          conteo = conteo + 1
          for (row in 1:nrow(new_data)){
            row_data <- new_data[row,]
            v_x <- seq(row_data$start, row_data$end)
            id <- row_data$id
            if (length(intersect(v_s, v_x)) != 0){
              local_table[local_table$id == id,]$nivel <- row_data$nivel - 1
              subdata[subdata$id == id, ]$nivel <- row_data$nivel - 1
            }
          }
          subdata <- subdata[subdata$nivel ==nivel, ]
          final_number_rows = nrow(subdata)
          if (number_rows != final_number_rows){
            break
          }
          if(number_rows == final_number_rows & conteo == number_rows){
            changed = FALSE
          }
        }
      }
      nivel = nivel - 1
    }
    
    return(local_table)
  })
  
  # Table
  output$e_table_4 <- DT::renderDataTable({
    
    # Updating click button
    #if (reactive_4$doPlot == FALSE) return()
    local_table <- e_location_table()
    if (is.null(local_table)) return(data.frame(sequence = character(),
                                                start = integer(),
                                                end = integer(),
                                                alleles = character(),
                                                count = integer()))
    local_table <- local_table[-1,c(4,2,3,7,8)]
    DT::datatable(local_table, rownames = FALSE,options= list(scrollX = "300px"))})
  
  # Plot
  
  t_plot_4 <- eventReactive(c(input$e_go_4),{
    ep_data5 <- e_parsed_data()
    protein_4 <- input$e_protein_4
    
    local_table <- e_location_table()
    if (is.null(local_table)) return()
    
    local_table$count <- as.numeric(local_table$count) 
    
    if (nrow(local_table) == 1){
      p <- ggplot()+
        geom_rect(data = local_table, aes(xmin = start,
                                          xmax = end,
                                          ymin = nivel -0.3,
                                          ymax = nivel + 0.3,
                                          fill = count,
                                          sequence = "-",
                                          text = paste("Class:", type,"<br>","Seq:", sequence)),
                  alpha = 0.9,
                  size = 0.3,
                  colour = "black")+ 
        geom_text(data = local_table[local_table$type == 'Protein',],
                  aes(x = start/2 + end/2, y = (nivel-0.3)/2 + (nivel+0.3)/2), label = 'protein_4')+
        scale_fill_gradientn(colors = c("yellow", "red"),
                             limits = c(0, max(local_table$count)),
                             na.value = "lightblue")+
        theme_bw(base_size = 12) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")+
        scale_y_continuous(expand = c(0,0), limits = c(-5, 1.5))
      ggplotly(p, tooltip = c("text")) %>% layout(dragmode = "zoom",showlegend = FALSE)
    }
    else{
      
      steps = c(5, 10, 25, 50, 100, 150, 2000, 250, 500, 1000, 1500, 2500, 5000)
      n_splits = 15
      max_x = local_table[local_table$type == "Protein",]$end
      step_x_1 = max_x / n_splits
      step_x = steps[which(abs(steps - step_x_1) == min(abs(steps - step_x_1)))]
      
      p <- ggplot()+
        geom_rect(data = local_table, aes(xmin = start,
                                          xmax = end,
                                          ymin = nivel-0.3,
                                          ymax = nivel+0.3,
                                          fill = count,
                                          text = paste("Seq: ", sequence, 
                                                       "<br>",
                                                       "Position: ", start, "-", end-1,
                                                       "<br>",
                                                       "Alleles: ", alleles)),
                  alpha = 0.9,
                  size = 0.3,
                  colour = "black")+
        geom_text(data = local_table[local_table$type == 'Protein',],
                  aes(x = start/2 + end/2, y = (nivel-0.3)/2 + (nivel+0.3)/2), label = protein_4)+
        scale_fill_gradientn(colors = c("yellow", "red"),
                             limits = c(0, max(local_table$count)),
                             na.value = "lightblue")+
        scale_y_continuous(expand = c(0,0), limits = c(min(local_table$nivel)-1, 1.5)) +
        scale_x_continuous(limits = c(-1, max_x) ,breaks = seq(0, max_x, by = step_x))+
        theme_bw(base_size = 12) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")
      ggplotly(p, tooltip = c("text")) %>% layout(dragmode = "zoom",showlegend = FALSE)    }
  })
  
  output$e_plot_4 <- renderPlotly({
    #if (reactive_4$doPlot == FALSE) return()
    t_plot_4()
  })
  
  # Download button
  output$e_downloadData4 <- downloadHandler(
    filename = function(){
      paste("Epitopes_Plot_Location_example.txt")
    },
    content = function(file){
      write.table(e_location_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  ###################
  
  # Promiscuity tab #
  ###################
  
  # Text
  output$e_text_5 <- renderUI({
    list(p(HTML("Epitope Promiscuity shows epitopes predicted to bind a minimum
    number of MHC alleles. The tool depicts a heatmap where alleles are on the x-axis
    and epitopes are on the y-axis. Users should set the strong and weak binding 
    cutoff %rank to filter the strong binder and weak binder epitopes. Strong binder epitopes (SB) are the 
    epitopes with a %rank less or equal than the strong binding cutoff %rank. 
    Weak binder epitopes (WB) are the epitopes with a %rank greater than the strong binding
    cutoff %rank but less or equal than the weak binding cutoff %rank. In the heatmap, 
    SB and WB epitopes are colored in red and orange, respectively. Users also 
    need to select a minimum number of MHC alleles to filter and show only epitopes 
    that are predicted to bind to at least this number of alleles. This tool also shows
    a table containing the amino acid sequence, the position within the protein, the ID 
    of the protein, and the number of alleles to which they are predicted to bind.
    <br>
    <h4> Parameters </h4>
    <ul>
     <li> <b>Minimum number of alleles:</b> This indicates the minimum number of MHC alleles that an epitope needs to bind to be shown in the heatmap.</li>    
     <li> <b>Strong Binding Cutoff % rank:</b> Epitopes with a %rank lower or equal than this cutoff are considered as strong binder epitopes (SB). </li>      
     <li> <b>Weak Binding Cutoff %rank:</b> Epitopes with a %rank lower or equal than this cutoff but greater than Strong Binding Cutoff are considered weak binder epitopes (WB). </li>
     </ul>")))})
  
  # Reactive
  e_reactive_5 <- reactiveValues(doPlot = FALSE)
  observeEvent(input$e_go_5, {e_reactive_5$doPlot <- input$e_go_5})
  e_wb_ct_react <- eventReactive(input$e_go_5, {input$e_wb_ct})
  e_sb_ct_react <- eventReactive(input$e_go_5, {input$e_sb_ct})
  e_min_al_react <- eventReactive(input$e_go_5, {input$e_min_al})
  
  # Parsing
  e_parsed_data5 <- eventReactive(c(input$e_go_5),{
    
    
    # Update click button
    #if (reactive_5$doPlot == FALSE) return()
    ep_data5 <- e_parsed_data()
    wb_ct = input$e_wb_ct
    sb_ct = input$e_sb_ct
    min_al = input$e_min_al
    
    # Parse
    ep_data5$prom <- rowSums(ep_data5[,seq(5, ncol(ep_data5))] < wb_ct)
    ep_data5 <- ep_data5 %>% filter(prom >= min_al) 
    ids <- ep_data5 %>% select(Peptide) %>% unlist()
    ep_data5 <- ep_data5[order(ep_data5[,"prom"], decreasing = TRUE),]
    ep_data5 <- ep_data5 %>% select("Peptide", "Pos", "Length", "ID", "prom")
    
  })
  
  # Table
  output$e_table_5 <- DT::renderDataTable({
    
    # Update click button
    #if (reactive_5$doPlot == FALSE) return()
    e_parsed_data5()
    
  }, rownames = FALSE, options = list(scrollX = "300px"))
  
  # Plot
  t_plot_5 <- eventReactive(c(input$e_go_5),{
    if(nrow(e_parsed_data5()) == 0 ) return()
    
    wb_ct = input$e_wb_ct
    sb_ct = input$e_sb_ct
    min_al = input$e_min_al
    
    # Parse
    ep_data5 <- e_parsed_data()
    ep_data5$prom <- rowSums(ep_data5[,seq(5, ncol(ep_data5))] < wb_ct)
    ep_data5 <- ep_data5 %>% filter(prom >= min_al) 
    ids <- ep_data5 %>% select(Peptide) %>% unlist()
    ep_data5 <- ep_data5[order(ep_data5[,"prom"], decreasing = TRUE),]
    asd<- factor(levels= ep_data5$Peptide,
                 labels = ep_data5$Peptide)
    mdata <- e_parsed_data() %>% filter( Peptide %in% ids)
    mdata <- mdata[order(mdata[,'Peptide']),]
    mdata <- mdata[!duplicated(mdata$Peptide),]
    mdata <- melt(mdata, id=c("Peptide", "Pos", "Length", "ID"))
    
    mdata$value2[mdata$value <= wb_ct & mdata$value > sb_ct] = "WB"
    mdata$value2[mdata$value <= sb_ct] = "SB"
    mdata$value2[mdata$value > wb_ct] = "NA"
    mdata$value2 <- as.factor(mdata$value2)
    
    p <-ggplot(mdata) + geom_tile(aes(x=variable, y=Peptide, fill=value2, cutoff = value), width = 0.9, height=0.9, color ="black", size= 0.3) +
      scale_fill_manual(values = c("NA" = "white", "SB"= "red", "WB" = "orange"))+
      ylim(levels(asd))+ xlab("Alleles") + ylab("")+
      theme_bw()+
      xlab('')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.title =element_blank())+ xlab("Alleles")
    ggplotly(p, tooltip = c("cutoff", "y"))
  })
  output$e_plot_5 <- renderPlotly({
    t_plot_5()
  })
  
  # Download button
  output$e_downloadData5 <- downloadHandler(
    filename = function(){
      paste("Epitope_promiscuity_example.txt")
    },
    content = function(file){
      write.table(e_parsed_data5(), file, sep = "\t",row.names = FALSE, quote = FALSE)
      
    })
  
  
  ####################
  
  ####################
  # Conservation tab #
  ####################
  
  # Text
  output$e_text_7 <- renderUI({
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

     </ul>")))})
  
  # Reactive
  e_reactive_7 <- reactiveValues(doPlot = FALSE)
  observeEvent(input$e_go_7, {e_reactive_7$doPlot <- input$e_go_7})
  
  e_allele_7_react <- eventReactive(input$e_go_7, {input$e_allele_7})
  e_cutoff_7_react <- eventReactive(input$e_go_7, {input$e_cutoff_7})
  e_protein_7_react <- eventReactive(input$e_go_7, {input$e_protein_7})
  e_plot_type_7_react <- eventReactive(input$e_go_7, {input$e_plot_type_7})
  e_conditional_7_react <- eventReactive(input$e_go_7, {input$e_conditional_7})
  
  
  e_table_epitopes <- eventReactive(c(input$e_go_7), {
    parsed = e_parsed_data() 
    allele_7 <- input$e_allele_7
    cutoff_7 <- input$e_cutoff_7
    protein_7 <- input$e_protein_7
    plot_type_7 <- input$e_plot_type_7
    conditional_7 <- input$e_conditional_7
  
    parsed = parsed %>% filter(ID %in% protein_7)
    parsed_1 = parsed %>% select(Peptide, allele_7) %>% unique()
    
    if (length(allele_7) != 1){
      if (conditional_7 == 'Intersection'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, max)
      } else if (conditional_7 == 'Union'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, min) 
      }
    }else{
      parsed_1$plot_score <- parsed_1[, allele_7]
    }
    
    candidate_peptides = parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
    parsed_2 = parsed %>% select(Peptide, ID, allele_7) %>% filter(Peptide %in% candidate_peptides)
    peptide_info = dcast(parsed_2[, c(1,2)], formula = Peptide~ID) %>% as.data.frame() 
    peptide_info$Proteins = apply(peptide_info[,-1], 1, function(x) return(paste(x[!is.na(x)], collapse = '-')))
    peptide_info = peptide_info %>% 
      group_by(Proteins) %>%
      summarize(sequences = paste(Peptide, collapse = '-')) %>%
      as.data.frame()
    peptide_info$N_sequences = apply(peptide_info %>% select(sequences), 1, function(x) return(length(unlist(strsplit(x, split = '-')))))
    return(peptide_info)
  })
  
  # Table
  output$e_table_7 <- DT::renderDataTable({
    # Updating Click button
    #if (reactive_7$doPlot == FALSE) return()
    
    e_table_epitopes()
  }, rownames = FALSE, options = list(scrollX = '300px'))
  
  
  t_plot_7_VD <- eventReactive(c(input$e_go_7), {
    allele_7 <- input$e_allele_7
    cutoff_7 <- input$e_cutoff_7
    protein_7 <- input$e_protein_7
    plot_type_7 <- input$e_plot_type_7
    conditional_7 <- input$e_conditional_7
    
    if(length(protein_7)> 6){
      return(NULL)
    }
    
    parsed = e_parsed_data() %>% filter(ID %in% protein_7)
    parsed_1 = parsed %>% select(Peptide, allele_7) %>% unique()
    
    if (length(allele_7) != 1){
      if (conditional_7 == 'Intersection'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, max)
      } else if (conditional_7 == 'Union'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, min) 
      }
    }else{
      parsed_1$plot_score <- parsed_1[, allele_7]
    }
    
    candidate_peptides = parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
    parsed_2 = parsed %>% select(Peptide, ID, allele_7) %>% filter(Peptide %in% candidate_peptides)
    
    myList = split(parsed_2$Peptide, parsed_2$ID)
    venn <- Venn(myList)
    data <- process_data(venn)
    
    return(ggplot() +
             geom_sf(aes(fill=count), data = venn_region(data)) +
             geom_sf(size = 1.2, color = "black", data = venn_setedge(data), show.legend = T) +
             geom_sf_text(aes(label = name), data = venn_setlabel(data), nudge_y = 0.1, size = 6.5) +
             geom_sf_label(aes(label=count), fontface = "bold", data = venn_region(data), label.size = NA, fill = NA, size = 7) +
             theme_void(base_size = 20) + scale_fill_gradient(low = 'white', high = 'red') +
             theme(legend.position = "none") +
             scale_x_continuous(expand = expansion(mult = 0.5)))
  })
  output$e_plot_7_VD <- renderPlot({
    t_plot_7_VD()
  })
  
  t_plot_7 <- eventReactive(c(input$e_go_7), {
    allele_7 <- input$e_allele_7
    cutoff_7 <- input$e_cutoff_7
    protein_7 <- input$e_protein_7
    plot_type_7 <- input$e_plot_type_7
    conditional_7 <- input$e_conditional_7
    
    parsed = e_parsed_data() %>% filter(ID %in% protein_7)
    parsed_1 = parsed %>% select(Peptide, allele_7) %>% unique()
    
    if (length(allele_7) != 1){
      if (conditional_7 == 'Intersection'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, max)
      } else if (conditional_7 == 'Union'){
        parsed_1$plot_score <- apply(parsed_1[, allele_7], 1, min) 
      }
    }else{
      parsed_1$plot_score <- parsed_1[, allele_7]
    }
    
    candidate_peptides = parsed_1$Peptide[parsed_1$plot_score <= cutoff_7]
    parsed_2 = parsed %>% select(Peptide, ID, allele_7) %>% filter(Peptide %in% candidate_peptides)
    
    new_data2 = parsed_2 %>% dcast(Peptide~ID)
    new_data2[,-1] = apply(new_data2[,-1], 2, function(x) return(as.integer(ifelse(is.na(x), 0, 1)))) %>% as.data.frame()
    colnames(new_data2) = gsub("\\|", "_", colnames(new_data2))
    
    protein = colnames(new_data2[,-1])
    
    i <- 1
    x = vector("list", 2**length(protein)-1)
    y = vector("list", 2**length(protein)-1)
    
    
    for(j in 1:length(protein)){
      combinations <- combn(protein,j)
      for (k in 1: ncol(combinations)){
        x[[i]] <- sum(colSums(get_intersect_members(new_data2, combinations[,k])))/j
        y[[i]] <- paste(combinations[,k], collapse = "-")
        i = i + 1
      }
    }
    z <- do.call(rbind, Map(data.frame, names = y, values= x))
    z <-z[z$values!=0,]
    
    nintersections = nrow(z)
    nalleles = length(protein)
    set_size <- data.frame(alleles = protein, size = colSums(new_data2[,2:ncol(new_data2)]))
    set_size <- set_size[order(set_size[, 'size']),]
    
    ordered_alleles <- as.character(set_size$alleles)
    dict <- vector("list")
    i = length(ordered_alleles)
    for (allele in ordered_alleles){
      dict[[allele]] <- i
      i = i -1
    }
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
    
    zz <- z[order(z[,"values"], decreasing = TRUE),]
    
    #### Union
    intersect_size_chart<-plot_ly(showlegend = FALSE) %>% add_trace(
      x = 1:nintersections, # Number of bars
      y = zz$values, # Height bar
      type = "bar",
      marker = list(color = "black",
                    hoverinfo = "none")) %>% 
      layout(xaxis = list(tickfont = list(size = 15)),
             yaxis = list(tickfont = list(size = 15))) %>% 
      add_trace(type = "scatter",
                mode = "text",
                x = 1:nintersections,
                y = zz$values+ max(zz$values) * 0.05,
                text = zz$values,
                textfont = list(color = "black", size = 15))
    
    ### Grid
    names_intersection <- as.character(zz$names)
    n_grids_rows = 0
    for (name in names_intersection){
      name_splitted <- unlist(strsplit(name, ","))
      n_grids_rows = n_grids_rows + length(name_splitted)
    }
    
    y_axis <- vector("numeric", n_grids_rows)
    x_axis <- vector("numeric", n_grids_rows)
    combo_axis <- vector("numeric", n_grids_rows)
    name_axis <- vector("numeric", n_grids_rows)
    j = 1
    x = 1
    for(name in names_intersection){
      name_splitted <- unlist(strsplit(name, "-"))
      for (i in 1:length(name_splitted)){
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
    
    grid <- plot_ly(
      type = "scatter",
      mode = "markers",
      marker = list(color = "lightgrey", size = 8)
    ) %>% add_trace(
      type = "scatter",
      x = rep(1:nintersections,
              length(protein)), 
      y = unlist(lapply(1:length(protein), function(x)
        rep(x - 0.5, nintersections))),
      hoverinfo = "none"
    ) %>% add_trace(
      data = group_by(lines2, combo),
      type = "scatter",
      mode = "lines+markers",
      x = lines2$x,
      y = lines2$y-0.5,
      hoverinfo = "text",
      text = ~name,
      line = list(color = "black", width = 3),
      marker = list(color = "black",
                    size = 10))%>% 
      layout(xaxis = list(showticklabels = FALSE,
                          showgrid = FALSE,
                          zeroline = FALSE),
             yaxis = list(showticklabels = FALSE,
                          showgrid = TRUE,
                          range = c(0, length(protein)),
                          zeroline = FALSE,
                          range = 1:10),
             margin = list(t = 0, b = 40))
    
    intersect_size_chart <-
      intersect_size_chart %>% layout(yaxis = list(title = "Intersections size"))
    
    s1 <-
      subplot(
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        set_size_chart,
        nrows = 2,
        widths = c(0.6, 0.4)
      )
    s2 <-
      subplot(intersect_size_chart,
              grid,
              nrows = 2,
              shareX = TRUE) %>% layout(showlegend = FALSE)
    subplot(s1, s2, widths = c(0.3, 0.7))
    
    
  })
  output$e_plot_7 <- renderPlotly({
    #if (reactive_7$doPlot == FALSE) return()
    t_plot_7()
    
    
  })
  
  # Showing plots
  observeEvent(input$e_plot_type_7, {
    req(input$e_plot_type_7)
    if (input$e_plot_type_7 == "Venn Diagram") {
      hide("e_plot_7")
      show("e_plot_7_VD")
    } else {
      hide("e_plot_7_VD")
      show("e_plot_7")
    }
  })
  
  # Showing plots
  observeEvent(input$e_go_7, {
    e_protein_7_react = input$e_protein_7
    e_plot_type_7 = input$e_plot_type_7
    if (length(e_protein_7_react)<=6 & e_plot_type_7 == 'Venn Diagram') {
      hide("e_plot_7")
      show("e_plot_7_VD")
    } else {
      hide("e_plot_7_VD")
      show("e_plot_7")
    }
  })
  
}

shinyApp(ui = ui, server = server)
