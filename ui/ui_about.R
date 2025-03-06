# ui_about.R - UI para la pestaña About
#
# Este archivo contiene la definición de la pestaña "About" que proporciona
# información general sobre Epitope-Evaluator.

#' UI para la pestaña About
#'
#' @return UI para la pestaña About
about_ui <- function() {
  tabPanel(title = span("About", style = 'font-size:130%'),
           sidebarLayout(
             sidebarPanel(
               width = 4,
               h3('Epitope-Evaluator'),
               wellPanel(
                 style = 'overflow-x:auto; height:70;', 
                 uiOutput('text_1_left')
               )
             ), 
             mainPanel(
               width = 8,
               h1('What is Epitope-Evaluator'),
               wellPanel(
                 uiOutput('text_1'),
                 imageOutput('image_1', width = '100%', height = '100%')
               )
             )
           )
  )
}

#' Servidor para la pestaña About
#'
#' @param output Objeto output de Shiny
#' @return NULL
about_server <- function(output) {
  # Texto del panel lateral
  output$text_1_left <- renderUI({
    "Epitope-Evaluator: an interactive web application to study predicted T-cell epitopes"
  })
  
  # Texto principal
  output$text_1 <- renderUI({
    github_code <- htmltools::a("GitHub", href = "https://github.com/SotoLF/Epitope-Evaluator")
    
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
    The R scripts of Epitope-Evaluator can be freely downloaded from", github_code, "and launched locally.
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
                       </ul>"))))
  })
  
  # Imagen
  output$image_1 <- renderImage({
    list(src = 'www/images/Github1.png', 
         contentType = 'image/png', 
         res = 300, 
         width = '100%',
         alt = "Epitope-Evaluator GitHub repository")
  }, deleteFile = FALSE)
}
