### Installation Packages
#install.packages("tidyselect")
#install.packages("shinydashboard")
#install.packages("rlist")
#detach("package:MASS", unload = TRUE)
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

#setwd("Z:/PPP/shiny_app")
#install.packages('rsconnect')
#rsconnect::setAccountInfo(name='lufesu',
#                          token='46F0C29680D4D3291BC4F11F053B61CB',
#                          secret='HXRRBuVER09Qz+FJxXqnQRyhQyKA/KtLYLkCsMe7')
#setwd("Z:/ppp")

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
  dashboardHeader(title="Epitope Analyzer",
                  titleWidth = 250),
  # SIDEBAR
  dashboardSidebar(
    conditionalPanel(condition="input.tabselected==1",
                     column(width = 12,
                            fileInput("file1", "Upload a txt file", accept = ".txt"))),
    conditionalPanel(condition="input.tabselected==2",
                     column(width = 12,
                            checkboxGroupInput('allele','MHC Allele', choices = c(), selected = c()),
                            numericInput("xmin", "Min value", min = 0, max = 10, value = 0),
                            numericInput("xmax", "Max value", min = 0, max = 10, value = 2),
                            numericInput("step", "Step", min = 0, max = 1, value = 0.1),
                            numericInput("n_prots", "Min number of proteins conserved", min = 1, max = 10, value = 2))),
    
    conditionalPanel(condition="input.tabselected==3",
                     column(width = 12,
                            numericInput("cutoff", "Cutoff", min = 0, max = 10, value = 2))),

    conditionalPanel(condition="input.tabselected==4",
                     column(width = 12,
                            selectInput("protein", "Protein", choices = c(), selected =c())),
                     column(width = 12,
                            numericInput("cutpoint", "Cutoff", min = 0, max = 10, value = 2)),
                     column(width = 12,
                            checkboxGroupInput('allele2','MHC Allele', choices = c(), selected = c())),
                     actionButton("run", "Run analysis")),
    
    conditionalPanel(condition ="input.tabselected==5",
                     column(width = 12,
                            numericInput("cutpoint2", "Cutoff", min = 0, max = 10, value = 2)),
                     column(width = 12,
                            numericInput("min_al", "Minimum number of alleles", min = 0, max = 10, value = 6)),
                     column(width = 12,
                            numericInput("sb_ct", "Strong Binding cutoff", min = 0, max = 10, value = 0.5)),
                     column(width = 12,
                            numericInput("wb_ct", "Weak Binding cutoff", min = 0, max = 10, value = 2))),
    conditionalPanel(condition="input.tabselected==6",
                     column(width = 12,
                            numericInput("cutoff_5", "Cutoff", min = 0, max = 10, value = 2)),
                     column(width = 12,
                            checkboxGroupInput("allele3", "MHC Allele", choices = c(), selected = c())),
                     actionButton("go", "Run analysis"))),
  # BODY
  dashboardBody(
    
    #mainPanel(
    tabsetPanel(
      tabPanel("Input files", value = 1,
               h4( "Input files"),
               wellPanel(DT::dataTableOutput("contents"))),
      
      tabPanel("Epitopes Distribution", value=2,
               h4("Distribution of epitopes"),
               wellPanel(uiOutput('text_1')),
               fluidRow(
                 column(5,plotlyOutput("plot1", height = 500)),
                 column(7,DT::dataTableOutput("table_1")),
                 column(2,downloadButton("downloadData2", "Download"), offset = 6)),
               fluidRow(column(12, plotlyOutput("plot1_1")))),
      
      tabPanel("Density Protein", value = 3,
               h4("Epitopes Density per Protein"),
               wellPanel(uiOutput('text_2')),
               fluidRow(column(6, plotlyOutput("plot2", height = 500)),
                        column(6, DT::dataTableOutput("table_2")),
                        column(2,downloadButton("downloadData3", "Download"), offset = 6)),
               fluidRow(column(12, plotlyOutput("plot4_4")))),
      
      tabPanel("Epitopes Plot Location", value=4,
               h4("Location of all epitopes over the protein"),
               wellPanel(uiOutput('text_3')),
               fluidRow(column(6, plotlyOutput("plot3", height = 200)),
                        column(6, DT::dataTableOutput("table_3")),
                        column(2,downloadButton("downloadData4", "Download"), offset = 6))),
      
      tabPanel("Epitope promiscuity", value=5,
               h4("Promiscuity of epitopes"),
               wellPanel(uiOutput('text_4')),
               fluidRow(column(6, plotlyOutput("plot4", height = 500)),
                        column(6, DT::dataTableOutput("table_4"))),
               fluidRow(column(2, downloadButton("downloadData5", "Download"), offset = 6))),
      
      tabPanel("Uset", value = 6,
               h4("Tab content"),
               plotlyOutput("plot5", height = "550px"),
               fluidRow(column(12, DT::dataTableOutput("table_6"))),
               fluidRow(column(2, downloadButton("downloadData6", "Download"), offset = 6))),
      
      id = "tabselected"
    )
  )
)

server <- function(session, input, output){
  
  #### Input Tab
  output$contents <- DT::renderDataTable({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    req(file)
    read_delim(file$datapath,"\t",
               escape_double = FALSE,
               trim_ws = TRUE)
    
  }, options= list(scrollX = "300px"))
  
  #### Epitope data table (First Tab)
  ep_data <- reactive({
    inFile <- input$file1
    mydata <-data.frame(read_delim(inFile$datapath,"\t",
                                   escape_double = FALSE,
                                   trim_ws = TRUE))
    updateCheckboxGroupInput(session = session, inputId = "allele", choices = colnames(mydata)[c(-1, -2, -3, -4)], selected = colnames(mydata)[c(-1, -2, -3, -4)][1])
    allele_example <- colnames(mydata)[5]
    if(startsWith(toupper(allele_example), "HLA")){
      updateNumericInput(session = session, inputId = "xmin", min = 0, max = 2, value = 0)
      updateNumericInput(session = session, inputId = "xmax", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "cutoff", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "cutpoint", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "cutpoint2", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "sb_ct", min = 0, max = 2, value = 0.5)
      updateNumericInput(session = session, inputId = "wb_ct", min = 0, max = 10, value = 2)
      updateNumericInput(session = session, inputId = "cutoff_5", min = 0, max = 10, value = 2)
    }else{
      updateNumericInput(session = session, inputId = "xmin", min = 0, max = 10, value = 0)
      updateNumericInput(session = session, inputId = "xmax", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutoff", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutpoint", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutpoint2", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "sb_ct", min = 0, max = 2, value = 2)
      updateNumericInput(session = session, inputId = "wb_ct", min = 0, max = 10, value = 10)
      updateNumericInput(session = session, inputId = "cutoff_5", min = 0, max = 10, value = 10)
    }
    updateNumericInput(session = session, inputId = "n_prots", min = 1, max = mydata %>% 
                         dplyr::select(ID) %>% 
                         unlist() %>% 
                         unique() %>% 
                         length(), value = 2)
    updateSelectInput(session = session, inputId = "protein", choices = mydata %>% dplyr::select(ID) %>% unlist() %>% unique(), selected = mydata["ID"][1,])
    updateCheckboxGroupInput(session = session, inputId = "allele2", choices = colnames(mydata)[c(-1, -2, -3, -4)], selected = colnames(mydata)[c(-1, -2, -3, -4)][1])
    updateCheckboxGroupInput(session = session, inputId = "allele3", choices = colnames(mydata)[c(-1, -2, -3, -4)], selected = colnames(mydata)[c(-1, -2, -3, -4)][1:5])
    updateNumericInput(session = session,
                       inputId = "min_al",
                       max = length(colnames(mydata)[c(-1, -2, -3, -4)]),
                       value = length(colnames(mydata)[c(-1, -2, -3, -4)])-1)
    return(mydata)
  })
  
  # Second tab
  sub_data <- reactive({
    
    if (length(input$allele != 1)){
      ep_data3 <- ep_data()
      ep_data3 <- ep_data3[order(ep_data3[,'Peptide']),]
      ep_data3 <- ep_data3[!duplicated(ep_data3$Peptide),]
      ep_data3 <- ep_data3 %>% select(input$allele)
      apply(ep_data3,1, max)
    }else{
      ep_data3 <- ep_data()
      ep_data3 <- ep_data3[order(ep_data3[,'Peptide']),]
      ep_data3 <- ep_data3[!duplicated(ep_data3$Peptide),]
      ep_data3$input$allele
    }
  })
  
  f_table <- reactive({
    if (length(input$allele != 1)){
      ep_data2 <- ep_data() %>% select(input$allele)
      ep_data4 <- ep_data()
      ep_data4$Max_score <- apply(ep_data2,1,max)
    }else{
      ep_data4 <- ep_data() %>% select(input$allele)
    }
    ep_data4 <- ep_data4 %>% filter(Max_score <= input$xmax & Max_score >= input$xmin)
    ep_data4$Length <- nchar(ep_data4$Peptide)
    ep_data4 <- ep_data4[,c(1,3,2,4)]
    return(ep_data4)
  })
  
  output$downloadData2 <- downloadHandler(
    filename = function(){
      paste("Epitope_Distribution.txt")
    },
    content = function(file){
      write.table(f_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  output$plot1 <- renderPlotly({
    sub_data <- sub_data()
    plot_ly(x = sub_data,
            type = "histogram",
            opacity = 0.8,
            hoverinfo = "y",
            marker = list(color = "red", line = list(cauto = FALSE,width= 1, color= "black", cmid = 2)),
            xbins = list(start = input$xmin, end = input$xmax, size = input$step),
            cumulative = list(enabled = TRUE)) %>% 
      layout(xaxis = list(title = "% Rank Predictor"),
             yaxis = list(title = "Number of epitopes"))})
  
  output$plot1_1 <- renderPlotly({
    table <- f_table() 
    x_peptides <- table %>% 
      select(Peptide, ID) %>% 
      unique() %>% 
      group_by(Peptide) %>% 
      summarise(conteo = n()) %>% 
      as.data.frame() %>% arrange(desc(conteo)) %>% 
      filter(conteo >= input$n_prots) %>% 
      select(Peptide) %>% unlist()
    
    p <- ggplot(table) + 
      geom_tile(aes(x = Peptide, y = ID, fill = "black"), alpha = 0.8 ,width = 0.9, height = 0.9, color = "black")+
      theme_classic(base_size = 12)+
      scale_fill_manual(values = c("black"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
            legend.position = "none") +
      xlab("")+
      ylab("") + 
      xlim(x_peptides)
    
    ggplotly(p)
    
    })
  
  output$table_1 <- DT::renderDataTable({
    f_table()},options= list(scrollX = "300px"))
  
  output$text_1 <- renderUI({
    list(p(HTML("Epitopes Distribution allows to see the distribution of the number of epitopes with less % binding affinity
                than <b> Max value </b> but greater than <b> Min value </b>. Users can also set the <b> MHC Allele</b>" )),
         p(HTML("<b>MHC Allele:</b> MHC alleles which were used in the prediction of epitopes. First allele is set by default. Multiple MHC Alleles can be chosen.
                <p>
                <b> Min value:</b> The lowest % binding affinity to be considered in the filter of epitopes. 0 by default.
                <p>
                <b> Max value:</b> The highest % binding affinity to be considered in the filter of epitopes. 2 by default.
                <p>
                <b> Step:</b> Value to set the width of bars. 0.1 by default.
                " )))})
  
  # Third tab
  
  d_table <- reactive({
    ep_data5 <- ep_data()
    ep_data5$min <- apply(ep_data5[seq(5, ncol(ep_data5))], 1, min)
    ep_data5$min <- as.numeric(ep_data5$min)
    ep_data5 <- ep_data5 %>% dplyr::filter(min <= input$cutoff) %>% dplyr::group_by(ID,Length) %>% dplyr::summarize(count=n())
    ep_data5$Density <- round(ep_data5$count/ep_data5$Length,2)
    return(ep_data5)
  })
  
  output$table_2 <- DT::renderDataTable({
    ep_data5 <- d_table()
    DT::datatable(ep_data5,options= list(scrollX = "300px"))})
  
  
  output$plot2 <- renderPlotly({
    ep_data5 <- d_table()
    s = input$table_2_rows_selected
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
                             "Density: ", round(ep_data5[!(ep_data5$ID %in% ids$ID),]$Density, 2)))%>%
      
      layout(title = "Correlation N° Epitopes and Proteins Length",
             yaxis = list(title = "N of epitopes",tickmode = "linear", dtick = y_step),
             xaxis = list(title = "Protein Length", tickmode = "linear", dtick = x_step))

  })
  
  output$downloadData3 <- downloadHandler(
    filename = function(){
      paste("Density_Protein.txt")
    },
    content = function(file){
      write.table(d_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  output$text_2 <- renderUI({
    list(p(HTML("Density Protein allows the know which proteins have more proportion of epitopes and see if there is a relatrion between the length of the protein and the number of predicted epitopes.
                Click on rows will highlight the protein." )),
         p(HTML("<b>Cutoff:</b> The highest % Binding Affinity to be considered as predicted epitope. 2 by default.
                " )))})
  
  # Fourth tab
  v_reactive <- reactiveValues(doPlot = FALSE)
  
  observeEvent(input$run, {
    v_reactive$doPlot <- input$run
  }, once = TRUE)
  
  alleles4 <- eventReactive(input$run, {
    input$allele2
  })
  cutoff4 <- eventReactive(input$run, {
    input$cutpoint
  })
  protein4 <- eventReactive(input$run, {
    input$protein
  })
  
  
  location_table <- reactive({
    
    peptide_counts = 
      ep_data() %>% 
      filter(ID == protein4()) %>% 
      select(Peptide, Pos, Length, all_of(alleles4())) %>% 
      melt(id.vars = c("Peptide", "Pos", "Length")) %>% 
      filter(value <= cutoff4()) %>% 
      group_by(Peptide) %>% 
      mutate(count = sum(value <= cutoff4()),
             alleles = paste(variable, collapse = " - ")) %>% 
      select(Peptide, count, alleles) %>% as.data.frame()

    ep_data2 = ep_data() %>%
      filter(ID == protein4()) %>% 
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

  output$downloadData4 <- downloadHandler(
    filename = function(){
      paste("Epitopes_Plot_Location.txt")
    },
    content = function(file){
      write.table(location_table(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  output$table_3 <- DT::renderDataTable({
    if (v_reactive$doPlot == FALSE) 
      return(DT::datatable(ep_data()[,3:4] %>% unique(), rownames = FALSE))
    
    local_table <- location_table()
    local_table <- local_table[-1,c(4,2,3,7,8)]
    DT::datatable(local_table, rownames = FALSE,options= list(scrollX = "300px"))})
  
  output$text_3 <- renderUI({
    list(p(HTML("Location of all epitopes are plotted along the full protein. Each epitope are represented as red bar while the protein is indicated by a blue bar.
                Additionally, some filters could be applied such as selecting the protein, the cutoff and the MHC Allele.
                The output is described with a table which contains information of the filtered epitopes and the plot of locaation epitopes")))})
  
  output$plot3 <- renderPlotly({
    if (v_reactive$doPlot == FALSE) return()
    
    local_table <- location_table()
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
  
  # Fifth tab
  ep_data4 <- shiny::reactive({
    ep_data5 <- ep_data()
    ep_data5$prom <- rowSums(ep_data5[,seq(5, ncol(ep_data5))] < input$cutpoint2)
    ep_data5 <- ep_data5 %>% filter(prom >=input$min_al) 
    ids <- ep_data5 %>% select(Peptide) %>% unlist()
    ep_data5 <- ep_data5[order(ep_data5[,"prom"], decreasing = TRUE),]
    ep_data5 <- ep_data5 %>% select("Peptide", "Pos", "Length", "ID", "prom")
    
  })
  output$table_4 <- DT::renderDataTable({
    ep_data4()
  },options= list(scrollX = "300px"))
  
  output$downloadData5 <- downloadHandler(
    filename = function(){
      paste("Epitope_promiscuity.txt")
    },
    content = function(file){
      write.table(ep_data4(), file, sep = "\t",row.names = FALSE, quote = FALSE)
      
    })
  
  output$plot4_4 <- renderPlotly({
    table <- ep_data()
    table_melt <- melt(table, id.vars = c("Pos", "Length", "ID", "Peptide"))
    protein_length <- table_melt %>% select(ID, Length) %>% unique()
    table_summarized <- table_melt %>% group_by(ID, variable) %>% filter(value <= input$cutpoint2) %>% 
      summarise(count = n()) %>% as.data.frame()
    index = match(table_summarized$ID, protein_length$ID)
    table_summarized$Length <- protein_length$Length[index]
    table_summarized$Density = table_summarized$count/table_summarized$Length
    p <- ggplot(table_summarized) + 
      geom_tile(aes(x = variable, 
                    y = ID, 
                    fill = Density, 
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
      scale_fill_gradientn(colors = c("white","black"),
                           limits = c(0, round(max(table_summarized$Density),2) +0.01))+
      xlab("") + ylab("")
    
    ggplotly(p, tooltip = c("text")) 
  })
  output$plot4 <- renderPlotly({
    ep_data5 <- ep_data()
    ep_data5$prom <- rowSums(ep_data5[,seq(5, ncol(ep_data5))] < input$cutpoint2)
    ep_data5 <- ep_data5 %>% filter(prom >= input$min_al) 
    ids <- ep_data5 %>% select(Peptide) %>% unlist()
    ep_data5 <- ep_data5[order(ep_data5[,"prom"], decreasing = TRUE),]
    asd<- factor(levels= ep_data5$Peptide,
                 labels = ep_data5$Peptide)
    mdata <- ep_data() %>% filter( Peptide %in% ids)
    mdata <- mdata[order(mdata[,'Peptide']),]
    mdata <- mdata[!duplicated(mdata$Peptide),]
    mdata <- melt(mdata, id=c("Peptide", "Pos", "Length", "ID"))
    
    mdata$value2[mdata$value <= input$wb_ct & mdata$value > input$sb_ct] = "WB"
    mdata$value2[mdata$value <= input$sb_ct] = "SB"
    mdata$value2[mdata$value > input$wb_ct] = "NA"
    mdata$value2 <- as.factor(mdata$value2)
    
    p <-ggplot(mdata) + geom_tile(aes(x=variable, y=Peptide, fill=value2, cutoff = value), width = 0.9, height=0.9, color ="black", size= 0.3) +
      scale_fill_manual(values = c("NA" = "white", "SB"= "red", "WB" = "orange"))+
      ylim(levels(asd))+ xlab("Alleles") + ylab("")+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.title =element_blank())+ xlab("Alleles")
    ggplotly(p, tooltip = c("cutoff", "y"))
    
  })
  output$text_4 <- renderUI({
    list(p(HTML("Epitope Promiscuity shows a table where the epitopes are filtered based on a <b>cutoff</b> and the number of alleles to be bound is indicated.
                It also shows a heatmap where Strong Binding epitopes <b>(SB)</b> and Weak Binding epitopes <b>(WB)</b> are colored in red and orange, respectively.
                SB and WB are set based on the <b>Strong Binder cutoff</b>  and <b>Weak Binder cutoff </b>." )),
         p(HTML("<b>Cutoff:</b> The highest % Binding Affinity to be considered as predicted epitope.
                <p>
                <b>Minimum number of alleles</b>: A minimum number of alleles can be chosen to not show unnecessary epitopes.
                <p>
                <b> Strong Binding Cutoff </b>: Epitopes with % Rank Binding affinity lower or equal than this cutoff are considered as SB. 
                <p>
                <b> Weak Binding Cutoff </b>: Epitopes with % Rank Binding affinity lower or equal than this cutoff but higher than <b> Strong Binding Cutoff </b> are considered as WB.
                " )))})
  
  # Sixth tab
  v <- reactiveValues(doPlot = FALSE)
  
  observeEvent(input$go, {
    v$doPlot <- input$go
  })
  
  alleles <- eventReactive(input$go, {
    ep_data5 <- ep_data()
    input$allele3
  })
  cutoff5 <- eventReactive(input$go, {
    input$cutoff_5
  })
  
  ep_data6 <- shiny::reactive({
    ep_data5 <- ep_data()
    alleles <- alleles()
    new_data <- data.frame(matrix(ncol = length(alleles)+1, nrow = nrow(ep_data5)))
    colnames(new_data) <- c("Peptide", alleles)
    new_data$Peptide<- ep_data5$Peptide
    for (allele in alleles){
      new_data[allele] <- ep_data5[allele] <= cutoff5()
      new_data[,allele] <- as.integer(new_data[,allele])}
    
    new_data <- new_data[!duplicated(new_data$Peptide),]
    i <- 1
    x = vector("list", 2**length(alleles)-1)
    y = vector("list", 2**length(alleles)-1)
    
    epitopes <- c()
    
    for(j in 1:length(alleles)){
      combinations <- combn(alleles,j)
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
  
  output$table_6 <- DT::renderDataTable({
    ep_data6()
  },options= list(scrollX = "300px"))
  
  output$downloadData6 <- downloadHandler(
    filename = function(){
      paste("Unput.txt")
    },
    content = function(file){
      write.table(ep_data6(), file, sep = "\t",row.names = FALSE, quote = FALSE)
    })
  
  output$plot5 <- renderPlotly({
    ep_data5 <- ep_data()
    alleles <- alleles()
    new_data <- data.frame(matrix(ncol = length(alleles)+1, nrow = nrow(ep_data5)))
    colnames(new_data) <- c("Peptide", alleles)
    new_data$Peptide<- ep_data5$Peptide
    for (allele in alleles){
      new_data[allele] <- ep_data5[allele] <= cutoff5()
      new_data[,allele] <- as.integer(new_data[,allele])}
    
    new_data <- new_data[order(new_data[,'Peptide']),]
    new_data <- new_data[!duplicated(new_data$Peptide),]
    i <- 1
    x = vector("list", 2**length(alleles)-1)
    y = vector("list", 2**length(alleles)-1)
    
    if (v$doPlot == FALSE) return()
    
    for(j in 1:length(alleles)){
      combinations <- combn(alleles,j)
      for (k in 1: ncol(combinations)){
        x[[i]] <- sum(colSums(get_intersect_members(new_data, combinations[,k])))/j
        y[[i]] <- paste(combinations[,k], collapse = ",")
        i = i + 1
      }
    }
    z <- do.call(rbind, Map(data.frame, names = y, values= x))
    z <-z[z$values!=0,]
    nintersections = nrow(z)
    nalleles = length(alleles)
    set_size <- data.frame(alleles = alleles, size = colSums(new_data[,2:ncol(new_data)]))
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
              length(alleles)), 
      y = unlist(lapply(1:length(alleles), function(x)
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
                          range = c(0, length(alleles)),
                          zeroline = FALSE,
                          range = 1:10),
             margin = list(t = 0, b = 40))
    
    intersect_size_chart <-
      intersect_size_chart %>% layout(yaxis = list(title = "Intersections size"))
    
    # The y axis labels of the
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
}
shinyApp(ui = ui, server = server)

