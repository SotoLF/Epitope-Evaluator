### Installation Packages
#install.packages("tidyselect")
#install.packages("shinydashboard")
#install.packages("rlist")
#install.packages("shinycssloaders")

### Libraries needed
library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(grid)
library(gridExtra)
library(reshape)
library(shinydashboard)
library(plotly)
library(tidyselect)
library(rlist)
library(dplyr)
library(tibble)
library(shinycssloaders)

# Options for Spinner
options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=1)
options(shiny.maxRequestSize = 30*1024^2)


# Function to calculate U-set 
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


ui <- dashboardPage(
  
  # HEADER
  dashboardHeader(title = "Epitope Analyzer",
                  titleWidth = 250),
  # SIDEBAR
  dashboardSidebar(
    conditionalPanel(condition = "input.tabselected==1",
                     column(width = 12,
                            fileInput("input_file", "Upload Input file", accept = ".txt"))), 
    conditionalPanel(condition="input.tabselected==2",
                     column(width = 12,
                            "Distribution Analysis:",
                            checkboxGroupInput('allele','MHC alleles', choices = c(), selected = c()),
                            radioButtons("conditional","Relation among MHC alleles", list("AND", "OR") ),
                            numericInput("xmin", "Min %rank", min = 0, max = 10, value = 0),
                            numericInput("xmax", "Max %rank", min = 0, max = 10, value = 2),
                            numericInput("step", "Step", min = 0, max = 1, value = 0.1),
                            radioButtons("conditional_plot","Plot type", list("Histogram", "Cumulative Histogram") ),
                            actionButton("go_2_1", "Run analysis"),
                            "Analysis of Presence/Absence of epitopes in proteins:",
                            numericInput("n_prots", "Min number of proteins", min = 1, max = 10, value = 2),
                            actionButton("go_2_2", "Plot duplicated epitopes"))),
    
    conditionalPanel(condition="input.tabselected==3",
                     column(width = 12,
                            "Density Analysis:",
                            numericInput("cutoff", "Cutoff %rank", min = 0, max = 10, value = 2),
                            actionButton("go_3_1", "Run analysis"),
                            "Plot customization:",
                            radioButtons("plot_type","Plot type: ", list("Heatmap", "Bar plot")),
                            radioButtons("fill_type","Color by: ", list("Epitopes Number", "Epitopes Density")),
                            actionButton("go_3_2", "Generate Plot"))),

    conditionalPanel(condition="input.tabselected==4",
                     column(width = 12,
                            "Location Analysis",
                            selectInput("protein_4", "Protein", choices = c(), selected =c())),
                     column(width = 12,
                            numericInput("cutoff_4", "Cutoff %rank", min = 0, max = 10, value = 2)),
                     column(width = 12,
                            checkboxGroupInput('allele_4','MHC alleles', choices = c(), selected = c())),
                     column(width = 12,
                            radioButtons("conditional_4","Relation among MHC alleles", list("AND", "OR") )),
                     column(width = 12, 
                            actionButton("go_4", "Run analysis"))),
    
    conditionalPanel(condition ="input.tabselected==5",
                     column(width = 12,
                            "Promiscuity Analysis",
                            numericInput("min_al", "Minimum number of alleles", min = 0, max = 10, value = 6)),
                     column(width = 12,
                            numericInput("sb_ct", "Strong Binding cutoff %rank", min = 0, max = 10, value = 0.5)),
                     column(width = 12,
                            numericInput("wb_ct", "Weak Binding cutoff %rank", min = 0, max = 10, value = 2)),
                     column(width = 12,
                            actionButton("go_5", "Run analysis"))),
    conditionalPanel(condition="input.tabselected==6",
                     column(width = 12,
                            "Up-Set Analysis",
                            numericInput("cutoff_6", "Cutoff %rank", min = 0, max = 10, value = 2)),
                     column(width = 12,
                            checkboxGroupInput("allele_6", "MHC alleles", choices = c(), selected = c())),
                     column(width = 12,
                            actionButton("go_6", "Run analysis")))),
  # BODY
  dashboardBody(
    
    #mainPanel(
    tabsetPanel(
      tabPanel("Input files", value = 1,
               wellPanel(style = "overflow-x:auto; height:200px;",uiOutput('text_1')),
               wellPanel(DT::dataTableOutput("contents_1"))),
      
      tabPanel("Epitope Distribution", value=2,
               h4("Distribution of epitopes by MHC allele"),
               wellPanel(style = "overflow-x:auto; height:200px;",uiOutput('text_2')),
               fluidRow(
                 column(5,withSpinner(plotlyOutput("plot_2_1", height = 500), type = 5)),
                 column(7,withSpinner(DT::dataTableOutput("table_2"), type = 5)),
                 column(2,downloadButton("downloadData2", "Download Table"), offset = 6)),
               fluidRow(column(12, withSpinner(plotlyOutput("plot_2_2"), type = 5)))),
      
      tabPanel("Epitope Density ", value = 3,
               h4("Correlation between number of epitopes and protein lengths"),
               wellPanel(style = "overflow-x:auto; height:200px;",uiOutput('text_3')),
               fluidRow(column(6, withSpinner(plotlyOutput("plot_3_1", height = 500), type = 5)),
                        column(6, withSpinner(DT::dataTableOutput("table_3"), type = 5)),
                        column(2,downloadButton("downloadData3", "Download Table"), offset = 6)),
               fluidRow(column(12, withSpinner(plotlyOutput("plot_3_2"), type = 5)))),
      
      tabPanel("Epitope Location", value=4,
               h4("Amino acid position of the epitopes in proteins"),
               wellPanel(style = "overflow-x:auto; height:200px;",uiOutput('text_4')),
               fluidRow(column(6, withSpinner(plotlyOutput("plot_4", height = 200), type = 5)),
                        column(6, withSpinner(DT::dataTableOutput("table_4"), type = 5)),
                        column(2,downloadButton("downloadData4", "Download Table"), offset = 6))),
      
      tabPanel("Epitope Promiscuity", value=5,
               h4("Number of alleles and binding rank per epitopes"),
               wellPanel(style = "overflow-x:auto; height:200px;",uiOutput('text_5')),
               fluidRow(column(6, withSpinner(plotlyOutput("plot_5", height = 500), type = 5)),
                        column(6, withSpinner(DT::dataTableOutput("table_5"), type = 5))),
               fluidRow(column(2, downloadButton("downloadData5", "Download Table"), offset = 6))),
      
      tabPanel("Epitope Intersection", value = 6,
               h4("Set of epitopes shared between different MHC allele combinations"),
               wellPanel(style = "overflow-x:auto; height:200px;",uiOutput('text_6')),
               withSpinner(plotlyOutput("plot_6", height = "550px"), type = 5),
               fluidRow(column(12, withSpinner(DT::dataTableOutput("table_6"), type = 5))),
               fluidRow(column(2, downloadButton("downloadData6", "Download Table"), offset = 6))),
      
      id = "tabselected"
    )
  )
)

server <- function(session, input, output){
  
  #############
  # Input tab #
  #############
  
  output$text_1 <-renderUI({
    github_code <- a("GitHub", href="https://github.com/SotoLF/Epitope-Analyzer")
    case_1_file <- a("File 1", href="https://raw.githubusercontent.com/SotoLF/Epitope-Analyzer/main/Examples/SARS_cov_CTLs_parsed.txt")
    case_2_file <- a("File 2", href="https://raw.githubusercontent.com/SotoLF/Epitope-Analyzer/main/Examples/SARS_cov_HTLs_parsed.txt")
    case_3_file <- a("case 3", href="https://raw.githubusercontent.com/SotoLF/Epitope-Analyzer/main/Examples/leishmania_HTLs.txt?token=AHHMUKJWPDHSUTX76XWNOTDBFMSLM")
    pmid1 <- a ("32711842", href= "https://pubmed.ncbi.nlm.nih.gov/32711842/")
    pmid2 <- a ("32406916", href= "https://pubmed.ncbi.nlm.nih.gov/32406916/")
    pmid3 <- a ("32308001", href= "https://pubmed.ncbi.nlm.nih.gov/32308001/")
    list(p(HTML(paste( "<b> Epitope Analyzer </b> is a Shiny app aimed to filter and analyze predicted T-cell epitopes. 
    Epitope Analyzer includes 5 tools: (1) Epitope Distribution, (2) Epitope Density, (3) Epitope Location, (4) Epitope Promiscuity, 
    and (5) Epitope Intersection. It needs a txt format file tab-delimited where the first 4 columns are the peptide sequence, the position
    of the peptide within the protein, the length of the protein, and the name of the protein. The following columns should be named
    with the allele containing the %rank score assigned to each peptide for each MHC allele. The <b>%rank score</b> is a number between
    0 and 100. It is defined as the rank of the predicted binding score compared to a set of random natural peptides by several predictors
    of T-cell epitopes (",pmid1,pmid2,pmid3, "). Usually, the output of these predictors requires to be parsed before uploading in the shiny app. Once
    the input file is uploaded, parameters need to be set by the user followed by clicking on the <b>Run analysis</b>  or <b> Plot</b>  buttons. 
    Each of the tools shows (1) graphics that can be downloaded as png files by clicking on the camera icon (<b> Download plot as png</b>), 
    and (2) tables that can be downloaded by clicking on the <b> Download table</b>  button.

                       <p>
                       We have included 2 example files that can be downloaded and used as input files on Epitope Analyzer.
                       <ul>
                       <li>", case_1_file, ": Predicted MHC Class I epitopes from SARS-CoV-2.</li>
                       <li>", case_2_file, ": Predicted MHC Class II epitopes from SARS-CoV-2</li>
                       </ul> 
                       The R scripts of Epitope Analyzer can be freely downloaded from", github_code ,"and launched locally.
                       <p>
                       <h4>How to cite Epitope Analyzer</h4>
                       Paper Citation
                       <p>
                       <p>
                       <h4> Contact </h4>
                       Juan Fuxman Bass - fuxman@bu.edu
                       <br>
                       Luis F Soto - lufesu98@gmail.com
                       "))))
  })
  # Uploading table
  output$contents_1 <- DT::renderDataTable({
    file <- input$input_file
    ext <- tools::file_ext(file$datapath)
    req(file)
    read_delim(file$datapath,"\t",
               escape_double = FALSE,
               trim_ws = TRUE)
    
  }, options = list(scrollX = "300px"))
  
  #################
  # Uploading tab #
  #################
  ep_data <- reactive({
    
    inFile <- input$input_file
    df_raw_file <-data.frame(read_delim(inFile$datapath,"\t",
                                   escape_double = FALSE,
                                   trim_ws = TRUE))
    
    ## Rename columns
    colnames(df_raw_file)[1:4] = c("Peptide", "Pos", "Length", "ID")
    
    ## Updating buttons
    updateCheckboxGroupInput(session = session, inputId = "allele", choices = colnames(df_raw_file)[c(-1, -2, -3, -4)], selected = colnames(df_raw_file)[c(-1, -2, -3, -4)][1])
    allele_example <- colnames(df_raw_file)[5]
    if(startsWith(toupper(allele_example), "HLA")){
      updateNumericInput(session = session, inputId = "xmin", min = 0, max = 2, value = 0)
      updateNumericInput(session = session, inputId = "xmax", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "cutoff", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "cutoff_4", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "sb_ct", min = 0, max = 2, value = 0.5)
      updateNumericInput(session = session, inputId = "wb_ct", min = 0, max = 10, value = 2)
      updateNumericInput(session = session, inputId = "cutoff_6", min = 0, max = 10, value = 2)
    }else{
      updateNumericInput(session = session, inputId = "xmin", min = 0, max = 10, value = 0)
      updateNumericInput(session = session, inputId = "xmax", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutoff", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutoff_4", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "sb_ct", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "wb_ct", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutoff_6", min = 0, max = 10, value = 10)
    }
    updateNumericInput(session = session, inputId = "n_prots", min = 1, max = df_raw_file %>% 
                         dplyr::select(ID) %>% 
                         unlist() %>% 
                         unique() %>% 
                         length(), value = 2)
    updateSelectInput(session = session, inputId = "protein_4", choices = df_raw_file %>% dplyr::select(ID) %>% unlist() %>% unique(), selected = df_raw_file["ID"][1,])
    updateCheckboxGroupInput(session = session, inputId = "allele_4", choices = colnames(df_raw_file)[c(-1, -2, -3, -4)], selected = colnames(df_raw_file)[c(-1, -2, -3, -4)][1])
    updateCheckboxGroupInput(session = session, inputId = "allele_6", choices = colnames(df_raw_file)[c(-1, -2, -3, -4)], selected = colnames(df_raw_file)[c(-1, -2, -3, -4)][1:5])
    updateNumericInput(session = session,
                       inputId = "min_al",
                       max = length(colnames(df_raw_file)[c(-1, -2, -3, -4)]),
                       value = length(colnames(df_raw_file)[c(-1, -2, -3, -4)])-1)
    
    return(df_raw_file)
  })

  
  ####################
  # Distribution tab #
  ####################
  
  # Text
  output$text_2 <- renderUI({
    list(p(HTML("Epitope Distribution shows (1) the distribution of the number of unique epitopes that were predicted to bind
    to selected MHC alleles within an %rank interval and (2) a heatmap showing the presence of epitopes across different proteins. 
    To obtain the distribution of epitopes, users should select one or more MHC alleles, set the min %rank and the max%rank to 
    consider a peptide as an epitope, and select the plot type which could be a histogram or cumulative histogram. 
    The histogram and cumulative histogram plots show the %rank intervals on the x-axis and the number of epitopes
    on the y-axis. The density is represented with a red line in the histogram plot. To show the epitopes present across
    multiple proteins, users need to indicate a minimum number of proteins in which the epitopes are contained.
    The Shiny app will show a heatmap where the epitopes are on the x-axis while proteins are on the y 
    axis. The cells are colored with blue or white indicating the presence or absence of the epitope, respectively.
    <br>
    <h4> Parameters </h4>
    <ul>
     <li> <b>Relation among MHC alleles:</b> When multiple alleles are selected, users must indicate if they are interested in shared epitopes that bind all the selected MHC alleles (AND) or epitopes that bind at least one of the MHC alleles selected (OR). </li>
     <li> <b>MHC alleles:</b> List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = first allele) </li>
    <li> <b>Min % rank:</b> Minimum %rank to filter the epitopes. (Default = 0)</li>
     <li> <b>Max % rank:</b> Maximum %rank to consider a peptide as predicted epitope. (Default MHC Class I = 2, MHC Class II = 10). </li>      
     <li> <b>Step:</b> Width of the bins in the histogram plot. (Default = 0.1)</li>
     <li> <b>Plot type:</b> The distribution of epitopes can be shown as a histogram or a Cumulative histogram. If histogram is selected, the density is represented with a red line.</li>  
     <li> <b>Min number of proteins:</b> Minimum number of proteins where epitopes are conserved (Default = 2). </li>     
     </ul>")))})
  
  # Reactive
  reactive_2 <- reactiveValues(doTable = FALSE, doPlot = FALSE)
  
  observeEvent(input$go_2_1, {
    reactive_2$doTable <- input$go_2_1
    reactive_2$doPlot <- FALSE})
  
  observeEvent(input$go_2_2, {
    reactive_2$doTable <- input$go_2_1
    reactive_2$doPlot <- input$go_2_2})
  
  allele_react <- eventReactive(input$go_2_1, {input$allele})
  xmin_react <- eventReactive(input$go_2_1, {input$xmin})
  xmax_react <- eventReactive(input$go_2_1, {input$xmax})
  step_react <- eventReactive(input$go_2_1, {input$step})
  conditional_plot_react <-eventReactive(input$go_2_1, {input$conditional_plot})
  conditional_react <- eventReactive(input$go_2_1, {input$conditional})
  n_prots_react <- eventReactive(input$go_2_2, {input$n_prots})
  
  # Parsing
  data_values <- reactive({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_2$doTable == FALSE) return()
    allele <- allele_react()
    conditional <- conditional_react()
    
    ep_data3 <- ep_data() %>% 
      dplyr::select(Peptide, allele) %>% 
      unique() 
    
    if (length(allele) != 1){
      if ( conditional == "AND"){
        ep_data3$plot_score <- apply(ep_data3[, allele], 1, max)
      } else if(conditional == "OR"){
        ep_data3$plot_score <- apply(ep_data3[, allele], 1, min)
      }
    } else{
      ep_data3$plot_score <- ep_data3[,allele]}
    
    return(ep_data3)
  })
  
  # Parsing
  f_table <- reactive({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_2$doTable == FALSE) return()
    allele <- allele_react()
    xmin <- xmin_react()
    xmax <- xmax_react()
    
    data_values <- data_values()
    tmp_data <- ep_data()

    keep = tmp_data$Peptide %in% data_values$Peptide[data_values$plot_score > xmin & data_values$plot_score <= xmax]
    data_keep = tmp_data[keep, c(1,2,4)]
    data_keep$Peptide_length = nchar(data_keep[,1])
    rownames(data_keep) <- NULL
    return(data_keep)
  })
  
  # Table
  output$table_2 <- DT::renderDataTable({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    ep_data3 <- ep_data()
    
    # Updating click button
    if (reactive_2$doTable == FALSE) return()
    f_table()}, rownames = FALSE, options = list(scrollX = "300px"))
  
  # Plot
  output$plot_2_1 <- renderPlotly({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_2$doTable == FALSE) return()
    allele <- allele_react()
    xmin <- xmin_react()
    xmax <- xmax_react()
    step <- step_react()
    conditional_plot <- conditional_plot_react()
    plot_values <- data_values()$plot_score
    
    if(conditional_plot == "Histogram"){
      cumulative = FALSE
      fit <- density(plot_values[plot_values <=xmax & plot_values > xmin])
      plot_ly(x = plot_values,
              type = "histogram",
              opacity = 0.8,
              hoverinfo = "y",
              marker = list(color = "ligthblue", line = list(cauto = FALSE,width= 1, color= "black", cmid = 2)),
              xbins = list(start = xmin, end = xmax, size = step),
              cumulative = list(enabled = FALSE)) %>% 
        add_trace(x = fit$x, y = fit$y, type = "scatter", mode = "lines", fill = "tozeroy", yaxis = "y2",
                  fillcolor = "transparent",
                  marker = list(color = "red")) %>% 
        layout(xaxis = list(title = "% Rank Predictor"),
               yaxis = list(title = "Number of epitopes"),
               yaxis2 = list(overlaying = "y", side = "right", position = 0.9),
               showlegend = FALSE)
    }else{
      cumulative = TRUE
      plot_ly(x = plot_values,
              type = "histogram",
              opacity = 0.8,
              hoverinfo = "y",
              marker = list(color = "ligthblue", line = list(cauto = FALSE,width= 1, color= "black", cmid = 2)),
              xbins = list(start = xmin, end = xmax, size = step),
              cumulative = list(enabled = cumulative)) %>%
        layout(xaxis = list(title = "% Rank Predictor"),
               yaxis = list(title = "Number of epitopes"))
    }
  
    })
  
  
  # Plot
  output$plot_2_2 <- renderPlotly({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_2$doPlot == FALSE) return()
    
    n_prots <- n_prots_react()
    
    table <- f_table() 
    x_peptides <- table %>% 
      select(Peptide, ID) %>% 
      unique() %>% 
      group_by(Peptide) %>% 
      summarise(conteo = n()) %>% 
      as.data.frame() %>% arrange(desc(conteo)) %>%
      filter(conteo >= n_prots) %>% 
      select(Peptide)
    
    if (nrow(x_peptides) == 0){
      y_proteins = ""
    }else{
      keep = table$Peptide %in% x_peptides$Peptide
      y_proteins = table[keep, ] %>% dplyr::select(ID) %>% unlist() %>% unique() 
    }
    
    
    p <- ggplot(table) + 
      geom_tile(aes(x = Peptide,
                    y = ID,
                    fill = "blue2",
                    text = paste("Protein:", ID,
                                 "<br>",
                                 "Epitope:", Peptide)), alpha = 0.7 ,width = 0.9, size = 1.2, height = 0.9, color = "black")+
      theme_classic(base_size = 12)+
      scale_fill_manual(values = c("blue2"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
            legend.position = "none") +
      xlab("")+
      ylab("") + 
      xlim(x_peptides$Peptide)+
      ylim(y_proteins)
    
    ggplotly(p, tooltip = c("text")) 
    
  })
  
  # Download button
  output$downloadData2 <- downloadHandler(
    filename = function(){
      paste("Epitope_Distribution.txt")
    },
    content = function(file){
      write.table(f_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  
  ###############
  # Density tab #
  ###############
  
  # Text input
  output$text_3 <- renderUI({
    list(p(HTML("Epitope Density shows a scatter plot where the protein length in amino acids 
    is on the x-axis and the number of epitopes is on the y-axis. Hovering over each point shows the 
    name of the protein, number of epitopes, length of the protein, and the epitope density. The epitope
    density is calculated as the number of epitopes divided by the length of the proteins. It also shows
    a heatmap representing the number (or density) of epitopes in each protein across each allele.
    <br>
    To obtain the scatter plot, users must indicate the maximum cutoff %rank to consider a peptide as an
    epitope and click on the <b>Run analysis</b> button. This will produce a scatter plot and a table containing the ID of the proteins,
    their length in amino acids, the number of epitopes, and the epitope density. Clicking on any row of the table will highlight 
    the point in the scatter plot. 
    <br>
    To obtain the heatmap, in addition to the maximum cutoff%rank, users must select 
    the plot-type (heatmap or bar plot) and the fill-type (by number or density of epitopes). We recommend using the 
    plot-type <b>heatmap</b> when several proteins are analyzed. When heatmap is selected, this shows the alleles on the x-axis 
    and the proteins on the y-axis. The color intensity indicates the log 10 of the number or the density of epitopes. When
    <b>barplot</b> is selected, each protein is represented as a bar while alleles are on the x-axis and the number (or density) 
    of epitopes are on the y-axis. Hovering over each bar (or cell in the heatmap) will show the protein ID, the allele, the 
    number of epitopes, the length of the protein, and the density of epitopes.

    <br>
    <h4> Parameters </h4>
    <ul>
     <li> <b>Cutoff % rank:</b> Maximum %rank to consider a peptide as predicted epitope. (Default MHC Class I = 2, MHC Class II = 10). </li>      
     <li> <b>Color by:</b> Fill heatmap by number or density of epitopes.(Default = Epitopes Number)</li>
     <li> <b>Plot type:</b> Heatmap is recommended for several proteins, while bar plot is recommended for few proteins. (Default = Heatmap)</li>     
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
  
  # Parsing data
  d_table <- reactive({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_3$doTable == FALSE) return()
    cutoff = cutoff_3_react()
    ep_data5 <- ep_data()
    
    
    ep_data5$min <- apply(ep_data5[seq(5, ncol(ep_data5))], 1, min)
    ep_data5$min <- as.numeric(ep_data5$min)
    ep_data5 <- ep_data5 %>% dplyr::filter(min <= cutoff) %>% dplyr::group_by(ID,Length) %>% dplyr::summarize(count=n())
    ep_data5$Density <- round(ep_data5$count/ep_data5$Length,2)
    return(ep_data5)
  })
  
  # Table
  output$table_3 <- DT::renderDataTable({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_3$doTable == FALSE) return()
    ep_data3 <- d_table()
    DT::datatable(ep_data3, rownames = FALSE, options= list(scrollX = "300px"))})
  
  # Plot
  output$plot_3_1 <- renderPlotly({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_3$doTable == FALSE) return()
    ep_data5 <- d_table()
    
    s = input$table_3_rows_selected
    y_min_value <- min(ep_data5$count)
    y_max_value <- max(ep_data5$count)
    x_min_value <- min(ep_data5$Length)
    x_max_value <- max(ep_data5$Length)
    x_step <- ceiling((x_max_value - x_min_value)/10)
    y_step <- ceiling((y_max_value - y_min_value)/10)
    ids = as.data.frame(ep_data5)[s, ]
    plot_ly(type = "scatter") %>% 
      add_trace(data = ep_data5[ep_data5$ID  %in% ids$ID,],
                x = ~as.factor(Length), y = ~as.factor(count),
                marker = list(size = 10,
                              color = "red",
                              line = list(color = "rgba(152, 0, 0, .8)",
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
                              line = list(color = "rgba(152, 0, 0, .8)",
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
             yaxis = list(title = "N of epitopes",tickmode = "linear", dtick = y_step),
             xaxis = list(title = "Protein Length", tickmode = "linear", dtick = x_step))
    
  })
  
  # Plot
  output$plot_3_2 <- renderPlotly({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_3$doPlot == FALSE) return()
    cutoff = cutoff_3_2_react()
    table <- ep_data()
    fill_type <- fill_type_react()
    
    table_melt <- melt(table, id.vars = c("Pos", "Length", "ID", "Peptide"))
    protein_length <- table_melt %>% select(ID, Length) %>% unique()
    table_summarized <- table_melt %>% group_by(ID, variable) %>% filter(value <= cutoff) %>% 
      summarise(count = n()) %>% as.data.frame()
    index = match(table_summarized$ID, protein_length$ID)
    
    if (fill_type == "Epitopes Number"){
      table_summarized$Length <- protein_length$Length[index] 
      table_summarized$Density <- table_summarized$count/table_summarized$Length
      table_summarized$type_fill <- log10(table_summarized$count)
    }
    else{
      table_summarized$Length <- protein_length$Length[index] 
      table_summarized$Density <- table_summarized$count/table_summarized$Length
      table_summarized$type_fill <- log10(table_summarized$count/table_summarized$Length)
    }
    
    y_ids <- table_summarized %>% group_by(ID) %>% summarise(accumulated = sum(10**type_fill)) %>% 
      arrange(desc(accumulated)) %>% select(ID) %>% unlist()
    x_ids <- table_summarized %>% group_by(variable) %>% summarise(accumulated = sum(10**type_fill)) %>% 
      arrange(desc(accumulated)) %>% select(variable) %>% unlist() %>% as.character()
    
    plot_type<- plot_type_react()
    if (plot_type == "Bar plot"){
      if(fill_type == "Epitopes Number"){
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
          theme_classic(base_size = 12)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
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
          theme_classic(base_size = 12)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                legend.title = element_blank())+
          xlab("") + ylab("") + 
          xlim(x_ids)
      }
      
    }else{
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
        theme_classic(base_size = 12)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              legend.position = "none")+
        scale_fill_gradientn(colors = c("white","blue"),
                             limits = c(min(table_summarized$type_fill)-0.1, round(max(table_summarized$type_fill),2) +0.1))+
        xlab("") + ylab("") +ylim(y_ids) + xlim(x_ids) 
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
  # Location tab #
  ################
  
  # Text
  output$text_4 <- renderUI({
    list(p(HTML("Epitope Location shows the location of epitopes within the protein (blue bar).
    Epitopes are represented on a color gradient from yellow to red reflecting the number of MHC
    alleles they are predicted to bind to. Hovering on each epitope will show the amino acid sequence,
    the position within the protein, and the MHC alleles that are predicted to recognize it. Users must
    select the protein, set a cutoff %rank to filter epitopes, and indicate the MHC alleles of interest.
    When multiple alleles are selected, users also need to specify if they are interested in epitopes that
    are predicted to bind all the selected alleles (AND) or at least one of the selected alleles (OR). 
    This tool also shows a table containing the sequence of the epitope, the start and end amino acid position
    within the protein, the alleles, and the number of alleles to which they were predicted to bind. This
    analysis may take a few minutes depending on the number of epitopes and the length of the protein.

    <br>
    <h4> Parameters </h4>
    <ul>
     <li> <b>Protein:</b> List of proteins obtained from the input file. Users must select a protein to be analyzed. (Default = first protein)</li>    
     <li> <b>Cutoff % rank:</b> Maximum %rank to consider a peptide as predicted epitope. (Default MHC Class I = 2, MHC Class II = 10). </li>      
     <li> <b>Condition over MHC alleles:</b> When multiple alleles are selected, users must indicate if they are interested in shared epitopes that bind all the selected MHC alleles (AND) or epitopes that bind at least one of the MHC alleles selected (OR). </li>
     <li> <b>MHC alleles:</b> List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = first allele) </li>
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
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    ep_data5 <- ep_data()
    
    if (reactive_4$doPlot == FALSE) return()
    allele_4 <- allele_4_react()
    cutoff_4 <- cutoff_4_react()
    protein_4 <- protein_4_react()
    conditional_4 <- conditional_4_react()

    # Parse
    peptide_counts = ep_data() %>% 
      filter(ID == protein_4) %>% 
      select(Peptide, Pos, Length, allele_4)
    
    if (length(allele_4) != 1){
      if ( conditional_4 == "AND"){
        peptide_counts$score <- apply(peptide_counts[, allele_4], 1, max)
      } else if(conditional_4 == "OR"){
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
    
    ep_data2 = ep_data() %>%
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
                             sequence = c("NA", ep_data2$Peptide),
                             nivel = c(1, rep(0, nrow(ep_data2))),
                             id = seq(0, nrow(ep_data2)),
                             alleles = c("NA", ep_data2$alleles),
                             count = c("Na", ep_data2$Value))
    
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
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_4$doPlot == FALSE) return()
    local_table <- location_table()
    if (is.null(local_table)) return()
    local_table <- local_table[-1,c(4,2,3,7,8)]
    DT::datatable(local_table, rownames = FALSE,options= list(scrollX = "300px"))})
  
  # Plot
  output$plot_4 <- renderPlotly({
    
    # Validation files
    validate(
      need(input$input_file != "", " ")
    )
    
    ep_data5 <- ep_data()
    
    # Updating click button
    if (reactive_4$doPlot == FALSE) return()
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
                  alpha = 0.7,
                  size = 0.1,
                  colour = "black")+ 
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
      n_splits = 10
      max_x = local_table[local_table$type == "Protein",]$end
      step_x = max_x / n_splits
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
                  alpha = 0.7,
                  size = 0.1,
                  colour = "black")+
        scale_fill_gradientn(colors = c("yellow", "red"),
                             limits = c(0, max(local_table$count)),
                             na.value = "lightblue")+
        scale_y_continuous(expand = c(0,0), limits = c(min(local_table$nivel)-1, 1.5)) +
        scale_x_continuous(limits = c(-1, max_x) ,breaks = seq(0, max_x,round(step_x)))+
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
    
    # Validation of input file
    validate(
      need(input$input_file != "", " ")
    )
    
    # Update click button
    if (reactive_5$doPlot == FALSE) return()
    ep_data5 <- ep_data()
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
    
    # Validation of input file
    validate(
      need(input$input_file != "", " ")
    )
    
    # Update click button
    if (reactive_5$doPlot == FALSE) return()
    parsed_data5()
    
  }, rownames = FALSE, options = list(scrollX = "300px"))
  
  # Plot
  output$plot_5 <- renderPlotly({
    
    # Validation of input file
    validate(
      need(input$input_file != "", " ")
    )
    
    # Update click button
    if (reactive_5$doPlot == FALSE) return()
    wb_ct = wb_ct_react()
    sb_ct = sb_ct_react()
    min_al = min_al_react()
    
    # Parse
    ep_data5 <- ep_data()
    ep_data5$prom <- rowSums(ep_data5[,seq(5, ncol(ep_data5))] < wb_ct)
    ep_data5 <- ep_data5 %>% filter(prom >= min_al) 
    ids <- ep_data5 %>% select(Peptide) %>% unlist()
    ep_data5 <- ep_data5[order(ep_data5[,"prom"], decreasing = TRUE),]
    asd<- factor(levels= ep_data5$Peptide,
                 labels = ep_data5$Peptide)
    mdata <- ep_data() %>% filter( Peptide %in% ids)
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
  
  ##############
  # Up-Set tab #
  ##############
  
  # Text
  output$text_6 <- renderUI({
    list(p(HTML("Epitope Intersection shows an Up-Set plot where the MHC alleles
    are represented as the sets and the epitopes are represented as the elements.
    The selected MHC alleles are represented with bars on the left side
    indicating the number of predicted epitopes for each MHC allele. 
    The number of epitopes for each MHC allele or intersection of MHC alleles is
    represented with bars at the top. Individual points in the grid indicate epitopes binding
    to a specific MHC allele while connected points indicate epitopes that can bind
    to multiple MHC alleles. This tool also provides a table containing the number
    of alleles, the alleles, the epitope sequences within each region, and their respective number of epitopes.

    <br>
    <h4> Parameters </h4>
    <ul>
     <li> <b>Cutoff %rank:</b> Maximum %rank to consider a peptide as predicted epitope. (Default MHC Class I = 2, MHC Class II = 10)</li>    
     <li> <b>MHC alleles:</b> List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = first 5 alleles) </li>      
     </ul>")))})
  
  # Reactive
  reactive_6 <- reactiveValues(doPlot = FALSE)
  observeEvent(input$go_6, {reactive_6$doPlot <- input$go_6})
  allele_6_react <- eventReactive(input$go_6, {input$allele_6})
  cutoff_6_react <- eventReactive(input$go_6, {input$cutoff_6})
  
  # Parse
  ep_data6 <- shiny::reactive({
    
    # Validation of input file
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_6$doPlot == FALSE) return()
    allele_6 <- allele_6_react()
    cutoff_6 <- cutoff_6_react()
    
    ep_data5 <- ep_data()
    new_data <- data.frame(matrix(ncol = length(allele_6)+1, nrow = nrow(ep_data5)))
    colnames(new_data) <- c("Peptide", allele_6)
    new_data$Peptide<- ep_data5$Peptide
    for (allele in allele_6){
      new_data[allele] <- ep_data5[allele] <= cutoff_6
      new_data[,allele] <- as.integer(new_data[,allele])}
    
    new_data <- new_data[!duplicated(new_data$Peptide),]
    i <- 1
    x = vector("list", 2**length(allele_6)-1)
    y = vector("list", 2**length(allele_6)-1)
    
    epitopes <- c()
    
    for(j in 1:length(allele_6)){
      combinations <- combn(allele_6,j)
      for (k in 1: ncol(combinations)){
        x[[i]] <- sum(colSums(get_intersect_members(new_data, combinations[,k])))/j
        epitopes <- c(epitopes, ep_data5[as.numeric(rownames(get_intersect_members(new_data, combinations[,k]))), 1])
        y[[i]] <- paste(combinations[,k], collapse = "-")
        i = i + 1
      }
    }
    
    z <- do.call(rbind, Map(data.frame, names = y, values= x))
    z <-z[z$values!=0,]
    
    epitopes2 <- c()
    counter <- 0
    for (i in z$values){
      epitopes2 <- c(epitopes2, paste(epitopes[seq(counter + 1, counter + i, by = 1)], collapse = "-"))
      print(epitopes2)
      counter <- counter + i
    }
    
    comb <- c()
    for (i in z$names){
      comb <- c(comb, length(unlist(strsplit(i, split = "-"))))
    }
    
    new <- data.frame(comb, z[, 1], epitopes2, z[, 2])
    colnames(new) <- c("N_alleles", "Alleles", "Epitopes", "N_epitopes")
    
    return(new)
  })
  
  # Table
  output$table_6 <- DT::renderDataTable({
    
    # Validation of input file
    validate(
      need(input$input_file != "", " ")
    )
    
    # Updating click button
    if (reactive_6$doPlot == FALSE) return()
    ep_data6()
    
  }, rownames = FALSE, options = list(scrollX = "300px"))
  
  # Plot
  output$plot_6 <- renderPlotly({
    
    # Validation of input file
    validate(
      need(input$input_file!= "", " ")
    )
    
    ep_data5 <- ep_data()
    
    # Updating click button
    if (reactive_6$doPlot == FALSE) return()
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
             yaxis = list(
               categoryarray = rev(set_size$alleles),
               categoryorder = "array"))
    
    zz <- z[order(z[,"values"], decreasing = TRUE),]
    
    #### Intersection
    intersect_size_chart<-plot_ly(showlegend = FALSE) %>% add_trace(
      x = 1:nintersections, # Number of bars
      y = zz$values, # Height bar
      type = "bar",
      marker = list(color = "black",
                    hoverinfo = "none")) %>% 
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
  
  # Download button
  output$downloadData6 <- downloadHandler(
    filename = function(){
      paste("Up_Set_epitopes.txt")
    },
    content = function(file){
      write.table(ep_data6(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
}
shinyApp(ui = ui, server = server)
