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

source("functionsShinyApp.R")

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
    fluidRow(
    infoBoxOutput("progressBox1"),
    infoBoxOutput("progressBox2")
    ),
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css?v=1")
  ,
    
    add_busy_spinner(spin = "atom", onstart = TRUE,position="full-page" ,timeout = 50),
      tags$style(HTML("
        .content-wrapper, .main-footer {
          margin-bottom: 0 !important;
        }
        .content {
          height: calc(100vh - 50px);
          overflow-y: auto;
        }
        .wrapper {
          height: 100vh;
          overflow: hidden;
        }
      ")),
   
    
    conditionalPanel(
      condition = "input.tab == 'infos'",
     
      box(width = 10,title="About this tool",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          id = 'Infos',
        
          HTML("
          <p><strong>Steps :</strong></p>
          <ol>
            <li>Load your FCS file(s) that you want to annotate.</li>
            <li>If your data is not compensated or transformed, proceed to the <em>Preprocessing</em> section.</li>
            <li>Finally, choose the algorithm you wish to use for annotation and visualization of the annotation results.</li>
            <li>Export the results.</li>
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
          box(width = 10,title="Choose one or multiple FCS files",
              
               status = "primary",
              
               solidHeader = TRUE,
              
               collapsible = TRUE,
              
              fileInput("fcsFile", NULL, multiple = TRUE, buttonLabel = "Upload...", accept = c('.fcs')),
              
              textOutput("notifLoading"),
             
          )
        )
      ),
      
    ),
    
    
    #### PREPROCESS
    
    
    conditionalPanel(
      
      condition = "input.tab == 'preprocessing'",

      uiOutput("fcsTransformed"),
      
      fluidPage(
        
        box(width = 10,HTML("<h4> If you want to apply compensation matrix and/or transform your dataset, else just ignore this part </h4>" )
 
          ),
        fluidRow(
          column(
            width = 12,
            box(width = 10,
              title = p("Preprocessing Settings",actionButton("HELP_PREPROCESSING","", icon("question-circle"), class="btn-xs")),
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              tabItems(
                tabItem(
                  tabName = "preprocessing",
                  
                  # Checkbox to apply compensation
                  checkboxInput("compensation", "Apply spillover matrix", value = TRUE),
                  uiOutput("compensation_message"),
                  
                  # Select markers to transform
                  selectizeInput("marker_untrans", "Markers to transform", choices = c(""), multiple = TRUE, width = 700),
                  
                  fluidRow(
                    column(
                      2,
                      tags$br(),
                      actionButton("clear_all_transform", "Clear All", style = "border-width: 2px; border-color: black;")
                    )
                  ),
                  
                  # Select transformation
                  fluidRow(
                    column(
                      2,
                      tags$br(),
                      selectInput("transfo", "Transformation", choices = list("logicle", "arcsinh"), selected = "logicle", width = 700)
                    )
                  ),
                  
                  # Set transformation argument
                  numericInput("trans_args", "Argument", value = 500, width = 700),
                  
                  # Apply transformation
                  actionButton(inputId = "submit", label = "Apply", style = " background-color: #B2DFEE;border-width: 2px; border-color: black;"),
                  
                  
                  # Apply transformation
                  actionButton(inputId = "unTransf", label = "Untransform all", style = "border-width: 2px; border-color: black;"),
                  
                  # Display messages
                  textOutput("message"),
                  # Select markers for visualization
                  selectInput("marker_x", "X-axis Marker", choices = c(""), selected = NULL, width = 700),
                  selectInput("marker_y", "Y-axis Marker", choices = c(""), selected = NULL, width = 700),
                  
                 
                  
                  # Display the plot
                  plotOutput("plot_density", height = "500px", width = "500px")
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
      box(width = 12,HTML("<h4> XGBoost is a machine learning algorithm used in this case to annotate the cell types present in the loaded files, by being trained on the basis of marker expression. An XGBoost model has been trained with a particular set of markers (features) that must match the markers in the fcs file to be annotated. </h4>" )),
     
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
          textInput("nameEnrichXGBOOST", "Please, enter the column name that will contain the XGBOOST enrichment :", "XGBoostPopID"),
           
          # Launch XGBoost annotation
          
          actionButton("annotate_data_with_XGboost", "Annotate Selected Files",style = "background-color: #B2DFEE ; border-width: 2px")
      )),
      
    ),
    
    #### ANNOTATION SCYAN
    
    conditionalPanel(
      
      condition = "input.tab == 'Scyan'",
      
      box(width = 10,HTML("<h4> Scyan algorithm need a knowledge table to annotate a ungated fcs file. The knowledge table contains well-known marker per expression. One row by cell type and one column by marker.
               Each expression should be beetween -1 and 1. -1 for negative expression and 1 for positive exression. </h4>
                   <h10> Quentin Blampey et al (2023) Briefings in Bioinformatics, Volume 24, Issue 5 <h10>" )
          ),
      
      box(width = 10,status = "primary",
          
          title="Scyan Annotation :",
          
          solidHeader = TRUE,
          
          collapsible = TRUE,
        
        # Load or build Scyan knowledge table  
        box(width = 10 , status="info",title=p('STEP 1 - Load or build scyan knowledge table' ,actionButton("HELP_KNOWLEDGE_TABLE","", icon("question-circle"), class="btn-xs")),
            
            
            fileInput("knowledgeTable", "Upload knowledge Table (csv or excel format) :", multiple = FALSE, buttonLabel = "Upload...", accept = c('.xlx','.xlsx','.csv')),
          
                     checkboxInput("oldKnowledgeTable", "OR use knowledge table of DYADEM project", value=FALSE)),
    
           
         
            #fileInput("buildKnowledgeTable", "   OR build knowledge table :", multiple = TRUE, buttonLabel = "Upload...", accept = c('.csv','.txt')),
        
        box(width = 10 , status="info",title=p('STEP 2 - Choose markers for annotation',actionButton("HELP_SCYAN_MARKERS","", icon("question-circle"), class="btn-xs")), 
            
            layout_columns(
              
              # Markers used for prediction  
              selectizeInput("markersUsedForPredictionsScyan", "Markers in fcs files", choices = c(""), multiple = TRUE),
              
              # Markers used for train the model
              selectizeInput("markersPresentInKnowledgeTable", "Markers in knowledge table", choices = c(""), multiple = TRUE),
     
            )
         
          ),
        
        # Annotation part 
        box(width = 12 , status="info",title=p('STEP 3 - Run algorithm',actionButton("HELP_SCYAN_RUN","", icon("question-circle"), class="btn-xs")),
            layout_columns(
            numericInput("std", "std", value = 0.25, width=700),
            numericInput("lr", "lr", value = 0.0001, width=700),
        ),
            
        textInput("nameEnrichScyan", "Please, enter the column name that will contain the SCYAN enrichment, default is:", "scyanPopID"),
        # Button to run Scyan annotation data
            actionButton("annotate_data_with_Scyan", "Run Scyan annotation",style = "background-color: #B2DFEE ; border-width: 2px")),
          
      )
    ),
    
    #### ANNOTATION SCAFFOLD
    
    conditionalPanel(
      condition = "input.tab == 'Scaffold'",
      
     
      box(
        width = 12, 
        title = "Annotation Scaffold",
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        
        # CLARA clustering
        box(
          width = 12, 
          status = "info", 
          title = p('STEP 1 - Clustering (CLARA)', 
                    actionButton("HELP_SCAFFOLD_CLUSTERING", "", icon("question-circle"), class = "btn-xs")
          ),
          
          # Checkbox to determine if files are already clustered
          checkboxInput("alreadyClustered", "Select if your file(s) is(are) already clustered", value = FALSE),
          
          # Show these inputs if "alreadyClustered" is unchecked
          conditionalPanel(
            condition = "input.alreadyClustered == false",
            
            # Button to choose markers used for clustering
            selectizeInput("marker_clustering", "Markers used for clustering", choices = c(""), multiple = TRUE),
            
            # Set k parameter
            numericInput("clustering_parameter", "k_parameter", value = 300),
            
            # Action button to run clustering
            actionButton("clustering", "Run CLARA clustering", style = "background-color: #B2DFEE; border-width: 2px")
          ),
          
          # Show this input if "alreadyClustered" is checked
          conditionalPanel(
            condition = "input.alreadyClustered == true",
            
            # Select input for the name of the column containing clusters
            selectInput("clusteringColumn", "Choose name of column that contains clusters", choices = c(""), multiple = FALSE, width = 500),
            
            # Button to choose markers used for scaffold map
            selectizeInput("marker_scaffoldmap", "Select markers for scaffold map", choices = c(""), multiple = TRUE),
            
          )
        ),
      
      
               
              
               
          
          
          # Scaffold Map
           box(width = 12 , status="info", title = p('STEP 2 - Build or load scaffold Map',actionButton("HELP_BUILD_MAP","", icon("question-circle"), class="btn-xs")),
               fileInput("landmarks", "Add landmarks FCS", multiple = TRUE, buttonLabel = "Upload...", accept = c('.fcs')),
               # Action button to build scaffod map
               actionButton("scaffoldMap", "Build Scaffold Map", style = "background-color: #B2DFEE; border-width: 2px;"),
               textOutput("messageBuildScaffoldMap"),
               tags$div(style = "margin-top: 20px;",  # Ajout d'un espacement
                        fileInput("oldScaffoldMap", " OR upload Scaffold Map :", 
                                  multiple = FALSE, 
                                  buttonLabel = "Upload ...", 
                                  accept = c('.scaffold'))
               )
               
               

           ),
           
          # Scaffold predictions
           box(width = 12 , status="info", title = p('STEP 3 - Annotation',actionButton("HELP_SCAFFOLD_ANNOTATION","", icon("question-circle"), class="btn-xs")),
               textInput("nameEnrichScaffold", "Please, enter the column name that will contain the SCAFFOLD enrichment:", "scaffoldPopID"),
               # Action Button to run scaffold predictions
               actionButton("RunScaffoldAnnotation", "Run Scaffold Annotation", style="background-color: #B2DFEE ; border-width: 2px"),
               
               textOutput("messageScaffoldAnnotation")
               
    
           )

      )
    ),

  # RESULTS
    
  conditionalPanel(
    condition = "input.tab == 'Stats'",
    titlePanel(""),  
    
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
             downloadButton('downloadAll', 'Download stats'),  
             textOutput("download")  
      )
    )
  )
  
  ,
  # Visualization Panel with Legend
  conditionalPanel(
    condition = "input.tab == 'Visualizations'",
    # Dropdown to select FCS file
    selectInput("enrichedFileViz", "Choose a File", choices = list(), selected = NULL),
    box(width = 8 , 
      title = "Run Algorithm",
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      
     
        column(
          width = 3,
          
          
          # Input for selecting markers (columns)
          selectizeInput("marker_viz", "Select markers", choices = c(), multiple = TRUE, width = "100%"),
          
          # Dropdown to select a column for coloring points
          selectInput("color_by", "Color by", choices = c(), selected = NULL),
          
          # Checkbox group to select steps (PCA, UMAP, t-SNE)
          checkboxGroupInput("step", "Run : ", choices = c("PCA", "UMAP", "t-SNE")),
          
          # Parameter settings for PCA
          conditionalPanel(
            condition = "input.step.includes('PCA')"
          ),
          
          # Parameter settings for UMAP
          conditionalPanel(
            condition = "input.step.includes('UMAP')",
            numericInput("umap_neighbors", "Number of Neighbors (UMAP)", value = 15, min = 2, max = 100),
            numericInput("umap_min_dist", "Minimum Distance (UMAP)", value = 0.1, min = 0.01, max = 1, step = 0.01)
          ),
          
          # Parameter settings for t-SNE
          conditionalPanel(
            condition = "input.step.includes('t-SNE')",
            numericInput("tsne_perplexity", "Perplexity (t-SNE)", value = 30, min = 5, max = 50),
            numericInput("tsne_theta", "Theta (Barnes-Hut Approximation) (t-SNE)", value = 0.5, min = 0, max = 1, step = 0.01),
            numericInput("tsne_iterations", "Max Iterations (t-SNE)", value = 500, min = 250, max = 2000)
          ),
          
          # Button to run the analysis
          actionButton("run", "Run", style = "background-color: #B2DFEE; border-width: 2px; border-color: black;")
        ),
        
        column(
          width = 9,
          
          # Tabs to display results for PCA, UMAP, and t-SNE
          tabsetPanel(
            tabPanel("PCA", 
                     sliderInput("sizePoint", "Points size", value = 0.5, min = 0.1, max = 2,step=0.1),
                     plotOutput("pcaPlot", height = "500px", width = "500px"),
                  
                     plotOutput("pcaDens", height = "500px", width = "500px"), # Density plot for PCA
                     uiOutput("pcaLegend")),
            tabPanel("UMAP", 
                     sliderInput("sizePoint", "Points size", value = 0.5, min = 0.1, max = 2,step=0.1),
                     plotOutput("umapPlot", height = "500px", width = "500px"),
                  
                     plotOutput("umapDens", height = "500px", width = "500px"), # Density plot for UMAP
                  
                     uiOutput("umapLegend")),
            tabPanel("t-SNE", 
                     sliderInput("sizePoint", "Points size", value = 0.5, min = 0.1, max = 2,step=0.1),
                     plotOutput("tsnePlot", height = "500px", width = "500px"),
                   
                     plotOutput("tsneDens", height = "500px", width = "500px"), # Density plot for t-SNE
                     
                     uiOutput("tsneLegend"))
          )
          
        )
      
    )
  )
  
  
  )
  
  
  
)
