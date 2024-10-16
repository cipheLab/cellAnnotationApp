
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinybusy)
library(bslib)
library(xgboost)
library(dplyr)
library(tidyr)
library(flowCore)
library(DT)
#library('FlowCIPHE',lib.loc = 'FlowCIPHE')


source("functionsShinyApp-copy.R")

options(expressions = 5e5, shiny.maxRequestSize = 100 * 1024 ^ 3)


ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = 'Cell Annotation App'),
  dashboardSidebar(
    sidebarMenu(
      id = "tab",
      menuItem(" Infos", tabName = "infos", icon = icon("info")),
      menuItem(" Upload FCS", tabName = "upload", icon = icon("upload")),
      menuItem(" Preprocessing", tabName = "preprocessing", icon = icon("adjust")),
      menuItem("Algorithms", tabName = "annotation_algorithms", icon = icon("pencil"),
               menuSubItem("XGboost", tabName = "XGboost"),
               menuSubItem("Scyan", tabName = "Scyan"),
               menuSubItem("Scaffold", tabName = "Scaffold")),
      menuItem("Results", tabName = "annotationResults", icon = icon("area-chart"),
               menuSubItem("Stats", tabName = "Stats"),
               menuSubItem("Visualizations", tabName = "Visualizations")
               
               ),
      
      actionButton("downloadFCS", "Download enriched FCS", icon = icon("download"))
    )
  ),
  dashboardBody(
    
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css?v=1")
  ,
    
    add_busy_spinner(spin = "atom", onstart = TRUE,position="full-page" ,timeout = 50),
    
  
    
    conditionalPanel(
      condition = "input.tab == 'infos'",
     
      box(title="About this tool",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          id = 'Infos',
        
          HTML("
          <p><strong>This tool allows you to : </strong></p>
          <ol>
            <li>Load one or more cytometry data files.</li>
            <li>Apply the necessary transformations.</li>
            <li>Annotate each cell in the dataset.</li>
          </ol>
          <p><strong>Steps :</strong></p>
          <ol>
            <li>Load your FCS file(s) that you want to annotate.</li>
            <li>If your data is not compensated or transformed, proceed to the <em>Preprocessing</em> section.</li>
            <li>Finally, choose the algorithm you wish to use for annotation and visualization of the annotation results.</li>
            <li>You can then export the results.</li>
          </ol>
        ")
      )
    ),
    
    ### UPLOAD
  
    
    conditionalPanel(
      
      condition = "input.tab == 'upload'",
      
      tabItems(
        
        tabItem(
          
          tabName = "upload",
          
         # Choose files to annotate
          box(title="Choose one or multiple FCS files",
              
               status = "primary",
              
               solidHeader = TRUE,
              
               collapsible = TRUE,
              
              fileInput("fcsFile", NULL, multiple = TRUE, buttonLabel = "Upload...", accept = c('.fcs')),
              
              textOutput("notifLoading"),
              
              DTOutput("ViewFCS")
          )
        )
      ),
      
    ),
    
    
    #### PREPROCESS
    
    
    conditionalPanel(
      
      condition = "input.tab == 'preprocessing'",

      uiOutput("fcsTransformed"),
      
      fluidPage(
        
        box(HTML("<h4> If you want to compensate and/or transform your dataset, else ignore this part </h4>" )
 
          ),
        fluidRow(
          column(
            width = 12,
            box(width = 9,
              title = "Preprocessing Settings",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              tabItems(
                tabItem(
                  tabName = "preprocessing",
                  
                  # Set if you want to compensate your data or not
                  checkboxInput("compensation", "Apply spillover matrix", value = TRUE),
                  
                  uiOutput("compensation_message"),
                  
                  # Set markers to transform
                  selectizeInput("marker_untrans", "Markers to transform", choices = c(""), multiple = TRUE, width=700),
                  
                  fluidRow(
                    
                    column(
                      2, 
                      tags$br(),
                      actionButton("clear_all_transform", "Clear All", style = "background-color: #B2DFEE; border-width: 2px; border-color: black;"),
                  
                    )
                  ),
                  
                  # Select transformation
                  
                  fluidRow(
                    
                    column(
                      2, 
                      tags$br(),
                  selectInput("transfo", "Transformation", choices = list("None", "logicle", "arcsinh"), selected="logicle", width=700))),
                  
                  # Set transformation argument
                  numericInput("trans_args", "Arg", value = 500, width=700),
                  
                  # Apply transformation
                  actionButton(inputId = "submit", label = "Apply", style = " background-color: #B2DFEE;border-width: 2px; border-color: black;"),
                  
                  textOutput("message")
                )
              )
            )
          )
        )
      )
      
    )
    
    ,
    
    #### ANNOTATION XGBOOST
    
    conditionalPanel(
      
      condition = "input.tab == 'XGboost' ", 
      
      # Description of XGBoost
      box(HTML("<h4> XGBoost is a machine learning algorithm used in this case to annotate the cell types present in the loaded files, by being trained on the basis of marker expression. An XGBoost model has been trained with a particular set of markers (features) that must match the markers in the fcs file to be annotated. </h4>" )),
     
      box( width = 12,  title = "Annotation XGboost ",
           
            status = "primary",
           
            solidHeader = TRUE,
           
            collapsible = TRUE,
           
           # Load or build XGBoost model
           box(width = 12 , status="info", title = 'STEP 1 - Load or choose XGBoost model and search commons markers',
            
            # Upload model
            fileInput("modelXGBoost", "Upload model RDS file :", multiple = TRUE, buttonLabel = "Upload ...", accept = c('.rds','RDS')),
            
            # Choose if you want a model XGBoost already existing
            radioButtons("oldModelXGBoost", "Or load existing model :", choices = list(""), selected=NULL),
            
            HTML("<h4> Choose markers for annotation </h4>"),
           
            layout_columns(
              
            # Markers used for prediction  
            selectizeInput("markersUsedForPredictions", "Markers in fcs files", choices = c(""), multiple = TRUE),
            
            # Markers used for train the model
            selectizeInput("markersUsedForModelTraining", "Markers in XGBoostmodel", choices = c(""), multiple = TRUE),
     
            ))
      
      ,
         # XGBoost predictions
          box(width = 12 , status="info", title = 'STEP 2 - Predictions',
            
           
          # Launch XGBoost annotation
          
          actionButton("annotate_data_with_XGboost", "Annotate Selected Files",style = "background-color: #B2DFEE ; border-width: 2px")
      )),
      
    ),
    
    #### ANNOTATION SCYAN
    
    conditionalPanel(
      
      condition = "input.tab == 'Scyan'",
      
      box(HTML("<h4> Scyan algorithm need a knowledge table to annotate a ungated fcs file. The knowledge table contains well-known marker per expression. One row by cell type and one column by marker.
               Each expression should be beetween -1 and 1. -1 for negative expression and 1 for positive exression. </h4>" )
          ),
      
      box(width = 10,status = "primary",
          
          title="Scyan Annotation :",
          
          solidHeader = TRUE,
          
          collapsible = TRUE,
        
        # Load or build Scyan knowledge table  
        box(width = 10 , status="info",title='STEP 1 - Load or build scyan knowledge table',
            
            
            fileInput("knowledgeTable", "Upload knowledge Table (CSV or Excel format) :", multiple = FALSE, buttonLabel = "Upload...", accept = c('.xlx','.xlsx','.csv')),
          
                     checkboxInput("oldKnowledgeTable", "OR use knowledge table of DYADEM project", value=FALSE)),
    
            conditionalPanel(
              condition = "input.oldKnowledgeTable == true",
              img(src = 'scyanImage.png', align = "center")
            ),
           
         
            #fileInput("buildKnowledgeTable", "   OR build knowledge table :", multiple = TRUE, buttonLabel = "Upload...", accept = c('.csv','.txt')),
        
        box(width = 10 , status="info",title='STEP 2 - Choose markers for annotation', 
            
            layout_columns(
              
              # Markers used for prediction  
              selectizeInput("markersUsedForPredictionsScyan", "Markers in fcs files", choices = c(""), multiple = TRUE),
              
              # Markers used for train the model
              selectizeInput("markersPresentInKnowledgeTable", "Markers in knowledge table", choices = c(""), multiple = TRUE),
     
            )
         
          ),
        
        # Annotation part 
        box(width = 12 , status="info",title='STEP 3 - Annotation',
            
            # Button to run Scyan annotation data
            actionButton("annotate_data_with_Scyan", "Run Scyan annotation",style = "background-color: #B2DFEE ; border-width: 2px")),
          
      )
    ),
    
    #### ANNOTATION SCAFFOLD
    
    conditionalPanel(
      condition = "input.tab == 'Scaffold'",
      
     
      box(width = 12,  title = "Annotation Scaffold  ",
           status = "primary",
           solidHeader = TRUE,
           collapsible = TRUE,
          
          # Clara clustering
           box(width = 12 , status="info", title = 'STEP 1 - Clustering',
               
               # Button to choose markers used for clustering
               selectizeInput("marker_clustering", "Markers used for clustering", choices = c(""), multiple = TRUE),
               
               # Set k parameter
               numericInput("clustering_parameter", "k_parameter", value = 300),
               
               # Action button to run clusteringp^ù
               actionButton("clustering", "Run CLARA clustering",style = "background-color: #B2DFEE ; border-width: 2px"),
               
               # Choose if you want a model XGBoost already existing
               checkboxInput("alreadyClustered", "My file(s) is(are) already clustered :", value=FALSE)),
              
               conditionalPanel(
                 
                 condition = "input.alreadyClustered == true",
                 
                 selectInput("clusteringColumn", "Choose name of column that contain clusters", choices = c(""), multiple = FALSE, width=500)

               ),
               
              
               
          
          
          # Scaffold Map
           box(width = 12 , status="info", title = 'STEP 2 - Build or load scaffold Map',
               
               # Add landmarks if you want to buils scaffold map
               selectizeInput("landmark", "Populations landmarks used for build the map", choices = list("DYADEM project landmarks", ""), multiple=FALSE),
               
               # Action button to build scaffod map
               actionButton("scaffoldMap", "Build Scaffold Map",style = "background-color: #B2DFEE ; border-width: 2px"),
               
               textOutput("messageBuildScaffoldMap"),
               
               # Upload model
               fileInput("oldScaffoldMap", " OR upload Scaffold Map :", multiple = FALSE, buttonLabel = "Upload ...", accept = c('.scaffold')),
               
               
               
               

           ),
           
          # Scaffold predictions
           box(width = 12 , status="info", title = 'STEP 3 - Annotation',
               
               # Action Button to run scaffold predictions
               actionButton("RunScaffoldAnnotation", "Run Scaffold Annotation", style="background-color: #B2DFEE ; border-width: 2px"),
               
               textOutput("messageScaffoldAnnotation")
               
    
           )

      )
    ),
    
  conditionalPanel(
    condition = "input.tab == 'Stats'",
    titlePanel(""),  # Si le titre n'est pas nécessaire, vous pouvez le supprimer.
    
    fluidRow(
      div(class = "right-aligned",
          selectInput("enrichedFile", "Choose file ", choices = list(""), selected = NULL, width = 800)
      )
    ),
    
    fluidRow(
      box(width = 4, title = 'XGBoost results', status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          DTOutput("resultsXGBoost")
      ),
      
      box(width = 4, title = 'Scyan results', status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          DTOutput("resultsScyan")
      ),
      
      box(width = 4, title = 'Scaffold results', status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          DTOutput("resultsScaffold")
      )
    ),
    
    fluidRow(
      column(6, 
             downloadButton('downloadAll', 'Télécharger tous les résultats en ZIP'),  
             textOutput("download")  
      )
    )
  )
  
  ,
  conditionalPanel(
    
    condition = "input.tab == 'Visualizations'",
    
    box(width=10, title = "Run algorithm",status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(div(class = "right-aligned",
                     selectInput("enrichedFile", "Choose file ", choices = list(""), selected = NULL, width=800)
        ),
          column(2,checkboxGroupInput("step","", choices=c("PCA", "MultiOptSNE","Rphenograph","PARC","UMAP","vaevictis", "phateR"))),
        
        ),
        column(2,tags$br(),actionButton("run","Run Pipeline")))
    
    
    
    
    )
  )
  
  
  
)
