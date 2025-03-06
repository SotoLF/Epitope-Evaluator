# ui_tutorial.R - UI para la pestaña Tutorial
#
# Este archivo contiene la definición de la pestaña "Tutorial" que proporciona
# videos tutoriales sobre el uso de la aplicación.

#' UI para la pestaña Tutorial
#'
#' @return UI para la pestaña Tutorial
tutorial_ui <- function() {
  tabPanel(title = span("Tutorial", style = 'font-size:130%'), 
           sidebarLayout(
             sidebarPanel(width = 4),
             mainPanel(
               width = 8,
               bsCollapse(
                 bsCollapsePanel(
                   'Epitope-Evaluator Sections',
                   wellPanel(
                     tags$video(
                       id = "video1", 
                       type = "video/mp4",
                       src = "videos/Sections.mp4", 
                       controls = "controls", 
                       width = "800px", 
                       height = "450px"
                     )
                   ),
                   style = 'primary'
                 )
               ),
               bsCollapse(
                 bsCollapsePanel(
                   'Upload Files',
                   wellPanel(
                     tags$video(
                       id = "video2", 
                       type = "video/mp4",
                       src = "videos/UploadFiles.mp4", 
                       controls = "controls", 
                       width = "800px", 
                       height = "450px"
                     )
                   ),
                   style = 'primary'
                 )
               ),
               bsCollapse(
                 bsCollapsePanel(
                   'Epitope-Distribution',
                   wellPanel(
                     tags$video(
                       id = "video3", 
                       type = "video/mp4",
                       src = "videos/EpitopeDistribution.mp4", 
                       controls = "controls", 
                       width = "800px", 
                       height = "450px"
                     )
                   ),
                   style = 'primary'
                 )
               ),
               bsCollapse(
                 bsCollapsePanel(
                   'Epitope-Intersection',
                   wellPanel(
                     tags$video(
                       id = "video4", 
                       type = "video/mp4",
                       src = "videos/EpitopeIntersection.mp4", 
                       controls = "controls", 
                       width = "800px", 
                       height = "450px"
                     )
                   ),
                   style = 'primary'
                 )
               ),
               bsCollapse(
                 bsCollapsePanel(
                   'Epitope-Density',
                   wellPanel(
                     tags$video(
                       id = "video5", 
                       type = "video/mp4",
                       src = "videos/EpitopeDensity.mp4", 
                       controls = "controls", 
                       width = "800px", 
                       height = "450px"
                     )
                   ),
                   style = 'primary'
                 )
               ),
               bsCollapse(
                 bsCollapsePanel(
                   'Epitope-Viewer',
                   wellPanel(
                     tags$video(
                       id = "video6", 
                       type = "video/mp4",
                       src = "videos/EpitopeLocation.mp4", 
                       controls = "controls", 
                       width = "800px", 
                       height = "450px"
                     )
                   ),
                   style = 'primary'
                 )
               ),
               bsCollapse(
                 bsCollapsePanel(
                   'Epitope-Promiscuity',
                   wellPanel(
                     tags$video(
                       id = "video7", 
                       type = "video/mp4",
                       src = "videos/EpitopePromiscuity.mp4", 
                       controls = "controls", 
                       width = "800px", 
                       height = "450px"
                     )
                   ),
                   style = 'primary'
                 )
               ),
               bsCollapse(
                 bsCollapsePanel(
                   'Epitope-Conservation',
                   wellPanel(
                     tags$video(
                       id = "video8", 
                       type = "video/mp4",
                       src = "videos/EpitopeConservation.mp4", 
                       controls = "controls", 
                       width = "800px", 
                       height = "450px"
                     )
                   ),
                   style = 'primary'
                 )
               )
             )
           ))
}
