# Load libraries
library(shiny)
library(xgboost)
library(plyr)
library(dplyr)
library(tidyr)
library(flowCore)
#library(FlowCIPHE)
library(shinydashboard)
library(shinyjs)
library(shinybusy)
library(bslib)
library(shinycssloaders)
library(gtools)
library(flowCore)
library(igraph)
library(Rcpp)
install.packages("reticulate")
library(reticulate)
library(openxlsx)
library(parallel)
library(stringr)
library(zip)



source_python("scyanFunctions.py")
s <- import("scyan")
source("functionsShinyApp.R")
source("functionsScaffold.R")


options(expressions = 5e5, shiny.maxRequestSize = 100 * 1024^3)

server <- function(input, output, session) {
  # Initialize objects
  listObject <- reactiveValues(
    listFCS = NULL,
    flow.frames = NULL,
    flow.frames.transformed = NULL,
    marker_untrans = NULL,
    markerXGBoost = NULL,
    model = NULL,
    models = NULL,
    flow.frames.enriched = NULL,
    flow.frames.tab = NULL,
    listLandmark = NULL,
    knowledgeTable = NULL,
    table=NULL,
    Map = NULL,
    resultsScyan = NULL,
    resultsScaffold = NULL,
    resultsXGboost = NULL,
    clusteringColumn = NULL
  )

  # Upload files

  observeEvent(input$fcsFile, {
    progress <- Progress$new()


    files <- input$fcsFile

    newfilesnames <- files$datapath


    filesname <- c(names(listObject$flow.frames), as.vector(newfilesnames))

    listObject$listFCS <- newfilesnames

    # Load and read fcs file one by one

    i <- 0

    new.flow.frames <- lapply(as.vector(listObject$listFCS), function(x) {
      i <<- i + 1

      progress$set(message = paste0("Reading file ... ", i, "/", length(listObject$listFCS), "."), value = i / length(as.vector(listObject$listFCS)))

      return(read.FCS(x, emptyValue = FALSE))
    })

    listObject$flow.frames <- c(listObject$flow.frames, new.flow.frames)

    # names of flow.frames
    names(listObject$flow.frames) <- files$name

    # Print the fcs names
    output$notifLoading <- renderText(paste0(names(unlist(listObject$flow.frames))))

    data <- as.data.frame(listObject$flow.frames[[1]]@exprs)

    results <- head(data)

  
  })



  observeEvent(input$clear_all_transform, {
    
    # Markers to transform
    updateSelectInput(session, "marker_untrans", selected = c("FSC-A", "SSC-A", "UV-BUV395-A", "UV-BUV661-A", "UV-BUV737-A", "V-BV421-A", "V-V500-A", "V-BV650-A", "V-BV711-A", "B-FITC-A", "B-PE-Cy5-5-A", "G-PE-Cy5-A", "G-PE-Cy7-A", "R-APC-A", "R-Alexa700-A", "R APC-Cy7-A"), choices = listObject$marker_untrans)

    shinyjs::show("markers_analyse")
  })



  # Launch Preprocessing

  observeEvent(input$submit, {
    
    # If fcs file is loaded

    if (!is.null(listObject$flow.frames)) {
      # Apply compensation or not

      if (input$compensation) {
        compens <- TRUE
      } else {
        compens <- FALSE
      }

      # Apply preprocessing
      listObject$flow.frames.transformed <- pree_process_fcs(listObject$flow.frames, input$trans_args, input$transfo, compens, input$marker_untrans)

      output$message <- renderText("Preprocessing applied successfully!")
    }
  })


  # Extract markers

  observe({
    # If fcs file is loaded
    if (!is.null(listObject$flow.frames)) {
      # Extract markers from fcs file
      listObject$marker_untrans <- extract_markers(listObject$flow.frames[[1]], NULL)

      # List of markers to transform
      updateSelectInput(session, "marker_untrans", choices = listObject$marker_untrans)
    }
  })

  # Extract markers for XGBoost

  observe({
    # If fcs file is loaded
    if (!is.null(listObject$flow.frames)) {
      # Load model trained on DYADEM project
      listObject$model <- readRDS("model_xgboost_landmark_internship_dyadem_project.rds")

      # Extract markers in the fcs file
      listObject$markersXGBoost <- extract_markers(listObject$flow.frames[[1]], NULL)

      # List of markers that match with model features
      updateSelectInput(session, "markersUsedForPredictions", choices = listObject$markersXGBoost)

      if (!is.null(listObject$model)) {
        # List of markers used by the model for training
        updateSelectInput(session, "markersUsedForModelTraining", choices = listObject$model$feature_names, selected = listObject$model$feature_names)
      }

      # List of markers that match with model features
      updateSelectInput(session, "markersUsedForPredictionsScyan", choices = listObject$markersXGBoost)

    }
  })


  # updateSelectInput(session,"exclude_markers_XGBoost",choices=listObject$markersXGBoost)

  observeEvent(input$annotate_data_with_XGboost, {
    progress <- Progress$new()

    # print(listObject$flow.frames.transformed)

    i <- 0
  
    if (is.null(listObject$flow.frames.transformed)) {
        listObject$flow.frames.transformed <- listObject$flow.frames
    }
    if (is.null(listObject$flow.frames.enriched)) {
        listObject$flow.frames.enriched <- listObject$flow.frames.transformed
      }
    


    # Predict annotation with XGBoost on all fcs file
    listObject$flow.frames.enriched <- lapply(listObject$flow.frames.enriched, function(x) {
      i <<- i + 1
      data <- preparMatrix(x, input$markersUsedForPredictions, input$markersUsedForModelTraining, NULL)

      # Scaling
      data <- scale(data)

      # Transform to XGBoost object
      data <- xgb.DMatrix(as.matrix(data))


      progress$set(message = paste0("Predict XGBoost annotations ... ", i, "/", length(listObject$flow.frames.transformed), "."), value = i / length(as.vector(listObject$flow.frames.transformed)))


      #  Prediction (annotation)
      predictions <- predict(listObject$model, data)

      # Transform to matrix
      predictions <- as.matrix(predictions)

      # Decaler les pr??dictions pour que ca matche avec les annotations scaffold
      predictions <- ifelse(predictions >= 10, predictions + 1, predictions)

      # Add 1 to all predictions
      predictions <- predictions + 1

      colnames(predictions) <- "popIDXGBoost"

      # Ajouter la nouvelle colonne dans le fichier fcs
      x <- enrich.FCS.CIPHE(x, predictions)

      return(x)
    })


    ## Display enrichment

    # After cell annotation :
    if (!is.null(listObject$flow.frames.enriched)) {
      # Initialize the result list


      # Show results for all files
      listObject$resultsXGBoost <- lapply(listObject$flow.frames.enriched, function(x) {
        # Convert to a dataframe
        x <- as.data.frame(x@exprs)

        # Add informations of all pop
        x <- table(x$popIDXGBoost)

        x <- as.data.frame(x)

        # Call the function to match label with popID
        popLabel <- matchPopIDD()

        labels <- c()

        for (i in x$Var1) {
          for (j in 1:29) {
            if (i == j) {
              labels <- c(labels, popLabel[[j]])
            }
          }
        }

        # Add new label column
        x$label <- labels

        # Add percentage column
        x <- x %>%
          mutate(Percentage = round(Freq / sum(Freq) * 100, 3))


        # dataframe that contains xgboost results
        results <- data.frame(popID = x$Var1, label = x$label, count = x$Freq, Percentage = x$Percentage)



        return(results)
      })

      # Add files names to each result
      names(listObject$resultsXGBoost) <- names(unlist(listObject$flow.frames))



      updateSelectInput(session, "enrichedFile", choices = names(unlist(listObject$flow.frames)))
    }
  })


  # Make list of XGBoost model

  observe({
    # List of model already saved in app
    models <- readRDS("model_xgboost_landmark_internship_dyadem_project.rds")


    models <- basename("model_xgboost_landmark_internship_dyadem_project.rds")

    # List of models already present
    listObject$models <- models

    # Create a list of models
    updateRadioButtons(session, "oldModelXGBoost", choices = listObject$models, selected = NULL)

    if (!is.null(listObject$flow.frames)) {
      listObject$model <- readRDS("model_xgboost_landmark_internship_dyadem_project.rds")

      # Extract marker in the fcs file
      listObject$markersXGBoost <- extract_markers(listObject$flow.frames[[1]], listObject$model)

      # Choose markers used for XGBoost Annotation
      updateSelectInput(session, "exclude_markers_XGBoost", choices = listObject$markersXGBoost, selected = c("FSC-A", "SSC-A", "UV-BUV395-A", "UV-BUV661-A", "UV-BUV737-A", "V-BV421-A", "V-V500-A", "V-BV650-A", "V-BV711-A", "B-FITC-A", "B-PE-Cy5-5-A", "G-PE-Cy5-A", "G-PE-Cy7-A", "R-APC-A", "R-Alexa700-A", "R APC-Cy7-A"))


      # Files to annotate
      files <- c("All", names(unlist(listObject$flow.frames)))

      # Update files to annotate
      updateRadioButtons(session, "FilesToAnnotate", choices = files, selected = "All")
    }
  })





  ## Clustering CLARA

  observeEvent(input$clustering, {
    # If there are file uploaded
    if (!is.null(listObject$flow.frames)) {
      if (is.null(listObject$flow.frames.transformed)) {
        # To avoid bug if no transformation is applied
        listObject$flow.frames.transformed <- listObject$flow.frames
      }

      # Apply clustering on all files

      # Apply clusterisation

      listObject$flow.frames.transformed <- claraClustering(listObject$flow.frames.transformed, input$clustering_parameter)
    }
  })

  # Display markers present in FCS file

  observe({
    if (!is.null(listObject$flow.frames)) {

      listObject$marker_untrans <- extract_markers(listObject$flow.frames[[1]], NULL)

      # Modify selectInput with the markers presents in the fcs file

      updateSelectInput(session, "marker_clustering", choices = listObject$marker_untrans, selected = c("FSC-A", "SSC-A", "UV-BUV395-A", "UV-BUV661-A", "UV-BUV737-A", "V-BV421-A", "V-V500-A", "V-BV650-A", "V-BV711-A", "B-FITC-A", "B-PE-Cy5-5-A", "G-PE-Cy5-A", "G-PE-Cy7-A", "R-APC-A", "R-Alexa700-A", "R APC-Cy7-A"))
    }
  })

  # Scaffold Map

  observeEvent(input$scaffoldMap, {
    # For each fcs file build a csv tab

    if (!is.null(listObject$flow.frames)) {
      if (is.null(listObject$flow.frames.transformed)) {
        listObject$flow.frames.transformed <- listObject$flow.frames
      }

      if (input$alreadyClustered == TRUE) {
        listObject$clusteringColumn <- input$clusteringColumn
      }
      # Apply on all transformed files

      listObject$flow.frames.tab <- lapply(listObject$flow.frames.transformed, function(x) {
        # Build csv tab

        return(builddCSVTab(x, listObject$clusteringColumn))
      })

      # Initialize gated.flow.frames
      gated.flow.frames <- NULL

      # List of landmark used by scaffold

      listLandmark <- list.files("/home/maelleWorkspace/GatedCleanLandmarkFSCTransformed", full.names = TRUE) 

      # Name of pop that you want to annotate

      namesLandmarks <- basename(listLandmark)

      # List of files you want to annotate

      listFCS <- as.character(listLandmark)[mixedorder(namesLandmarks)]

      # Set landmarks names
      temp.names <- namesLandmarks

      # Add names for each landmark
      filesname <- c(names(gated.flow.frames), as.vector(temp.names))

      # Loop for all files
      i <- 0

      new.flow.frames <- lapply(as.vector(listFCS), function(x) {
        i <<- i + 1

        fcs <- read.FCS(x, emptyValue = FALSE)

        if (dim(fcs)[1] == 1) {
          fcs@exprs <- fcs@exprs[c(1, 1, 1), ]
        }

        return(fcs)
      })

      # Landmarks
      gated.flow.frames <- c(gated.flow.frames, new.flow.frames)

      # Add names to each landmark
      names(gated.flow.frames) <- filesname

      # Set markers used for scaffold annotation
      col.names <- input$marker_clustering

      # Initialize ungated file
      clusteredFiles <- NULL

      # Contain output of builddCSVTab function
      clusteredFiles <- c(listObject$flow.frames.tab)

      # Assign names of ungated files to the list of files
      names(clusteredFiles) <- names(listObject$flow.frames)

      # Assign NA
      map.clusteredFiles.names <- names(listObject$flow.frames)

      names(gated.flow.frames) <- filesname

      clustedFiles <- clusteredFiles


      # To increase C stack limit
      sourceCpp("forceatlas2.cpp")

      # Run scaffold analysis
      result <- runn_analysis_gated(gated.flow.frames, # Landmarks
        clustedFiles, # Clustered files to annotate
        outputDir = getwd(), # Output directory
        map.clusteredFiles.names, # Clustered files names
        FALSE, # Boolean
        col.names, # Markers used for scaffold annotation
        col.names, # Markers used for annotation
        col.names.inter_cluster = NULL,
        ew_influence = NULL,
        inter_cluster.weight_factor = 0.7,
        inter.cluster.connections = TRUE,
        overlap_method = "repel"
      )

      # Show notification chen map is created
      output$messageBuildScaffoldMap <- renderText("Scaffold map created successfully !")

      # Assign scaffold map to the variable listObject$Map
      listObject$Map <- result
    }
  })


  # Run scaffold annotation
  observeEvent(input$RunScaffoldAnnotation, {
    # If scaffold Map is created
    if (!is.null(listObject$Map)) {
      if (is.null(listObject$flow.frames.enriched)) {
        listObject$flow.frames.enriched <- listObject$flow.frames.transformed
      }


      # Assign scaffold Map to scaffold variable
      scaffold <- listObject$Map


      # Name of file ungated
      list2 <- names(unlist(listObject$flow.frames))

      list1 <- names(unlist(listObject$flow.frames))

      # Initialize list that contains list of files ungated
      listFiles <- NULL
      print(listObject$flow.frames.enriched)
      # List of files ungated

      # Run scaffold annotation
      print("listFiles")
      print(listFiles)

      print("scaffold")
      print(scaffold)

      res <- scaffold_events_export(list1, list2, listObject$flow.frames.enriched, scaffold, "CLARABIS")


      listObject$flow.frames.enriched <- res

      listObject$resultsScaffold <- lapply(listObject$flow.frames.enriched, function(x) {
        # Convert new enriched fcs expression matrix to a dataframe
        resultScaffold <- as.data.frame(x@exprs)

        # Convert popScaffoldID results to a dataframe
        resultScaffold <- as.data.frame(table(resultScaffold$popIDscaffoldBis))

        # Add percentage column
        resultScaffold <- resultScaffold %>%
          mutate(Percentage = round(Freq / sum(Freq) * 100, 3))


        # Call the function to match label with popID

        popLabel <- matchPopIDD()

        labels <- c()

        for (i in resultScaffold$Var1) {
          for (j in 1:29) {
            if (i == j) {
              labels <- c(labels, popLabel[[j]])
            }
          }
        }

        # Add new label column
        resultScaffold$label <- labels

        # dataframe that contains xgboost results
        results <- data.frame(popID = as.vector(resultScaffold$Var1), label = as.vector(resultScaffold$label), count = resultScaffold$Freq, Percentage = resultScaffold$Percentage)

        # Sort by popID number
        results <- results[order(as.numeric(results$popID)), ]

        return(results)
      })

      names(listObject$resultsScaffold) <- names(unlist(listObject$flow.frames))

      output$messageScaffoldAnnotation <- renderText("Scaffold annotation is done, see Visualization section !")

      updateSelectInput(session, "enrichedFile", choices = names(unlist(listObject$flow.frames)))
    }
  })

  # SCYAN annotation (python module)
  
  observeEvent(input$knowledgeTable, {
    # Read the knowledge table
    
    listObject$table <- readKnowledgeTable(input$knowledgeTable$datapath)
   
    # List of markers used by the model for training
    updateSelectInput(session, "markersPresentInKnowledgeTable", choices = colnames(listObject$table), selected = colnames(listObject$table))
    
  })
  
  observe({
    if (input$oldKnowledgeTable) {
    
      listObject$table <- readKnowledgeTable("knowledgeTableScyanWithNI.csv")
     
      updateSelectInput(session, "markersPresentInKnowledgeTable", 
                        choices = colnames(listObject$table), 
                        selected = colnames(listObject$table))
    } else {
      updateSelectInput(session, "markersPresentInKnowledgeTable", 
                        choices = NULL, 
                        selected = NULL)  
    }
  })
  observeEvent(input$annotate_data_with_Scyan, {
    progress <- Progress$new()

    progress$set(message = "Run Scyan algorithm ... ")

    # Knowledge table :

    # If you want to load a knowledge table (csv format)

# 
#     # If you want to build a knowledge table
# 
#     observeEvent(input$buildKnowledgeTable, {
# 
# 
# 
# 
#     })

    if (is.null(listObject$flow.frames.transformed)) {
      
      listObject$flow.frames.transformed <- listObject$flow.frames
    }




    listObject$resultsScyan <- lapply(listObject$flow.frames.transformed, function(x) {
      # Prepare matrix with the good marker
      data <- preparMatrix(x, input$markersUsedForPredictionsScyan, input$markersPresentInKnowledgeTable, "scyan")
      
      # FCS to annotate

      adata <- builAnnDataObject(data)

      # Scale the data
 
      s$preprocess$scale(adata)
      
     
      # Build the model
      model <- runModel(adata, listObject$table)
      
      # Predict the model on adata
      test(model)

      # Obs for adata object
      resultsScyan <- table(as.data.frame(adata$obs$scyan_pop))

      # Convert to dataframe
      resultsScyan <- as.data.frame(resultsScyan)


      # Extract popID
      popID <- as.vector(sapply(as.vector(resultsScyan$adata.obs.scyan_pop), function(x) {
        return(strsplit(x, "_")[[1]][1])
      }))

      # Extract label
      label <- as.vector(sapply(as.vector(resultsScyan$adata.obs.scyan_pop), function(x) {
        return(strsplit(x, "_")[[1]][2])
      }))

      # Add Scyan popID
      resultsScyan$popID <- popID

      # Add labels
      resultsScyan$label <- label

      # Add percentage column
      resultsScyan <- resultsScyan %>%
        mutate(Percentage = round(Freq / sum(Freq) * 100, 3))


      # dataframe that contains xgboost results
      results <- data.frame(popID = as.vector(resultsScyan$popID), label = as.vector(resultsScyan$label), count = resultsScyan$Freq, Percentage = resultsScyan$Percentage)

      # Sort by popID number
      results <- results[order(as.numeric(results$popID)), ]

      return(results)
    })
    names(listObject$resultsScyan) <- names(unlist(listObject$flow.frames))

    updateSelectInput(session, "enrichedFile", choices = names(unlist(listObject$flow.frames)))
  })


  # Show result of selected files
  observe({
    

    if (!is.null(listObject$resultsXGBoost)) {
      
      # Convert result into datatable format
      resultsXGBoost <- DT::datatable(listObject$resultsXGBoost[[input$enrichedFile]], rownames = FALSE)
      output$resultsXGBoost <- renderDT(resultsXGBoost)
    }
    if (!is.null(listObject$resultsScyan)) {
      resultsScyan <- DT::datatable(listObject$resultsScyan[[input$enrichedFile]], rownames = FALSE)
      output$resultsScyan <- renderDT(resultsScyan)
    }
    if (!is.null(listObject$resultsScaffold)) {
      resultsScaffold <- DT::datatable(listObject$resultsScaffold[[input$enrichedFile]], rownames = FALSE)
      output$resultsScaffold <- renderDT(resultsScaffold)
    }
  })
  
  
    
  output$downloadAll <- downloadHandler(
    filename = function() {
      paste0("All_Statistics_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Create temporary repertory
      temp_dir <- tempdir()
      temp_files <- c()
      

      for (flow_frame_name in names(unlist(listObject$flow.frames))) {
   
        wb <- createWorkbook()
        if (!is.null(listObject$resultsScyan)) {
          addWorksheet(wb, "Scyan")
          Scyan <- listObject$resultsScyan[[flow_frame_name]]
          writeData(wb, sheet = "Scyan", Scyan)
        }
        
        if (!is.null(listObject$resultsXGBoost)) {
          addWorksheet(wb, "XGBoost")
          XGBoost <- listObject$resultsXGBoost[[flow_frame_name]]
          writeData(wb, sheet = "XGBoost", XGBoost)
        }
        
        if (!is.null(listObject$resultsScaffold)) {
          addWorksheet(wb, "Scaffold")
          Scaffold <- listObject$resultsScaffold[[flow_frame_name]]
          writeData(wb, sheet = "Scaffold", Scaffold)
        }
        
        
        temp_file <- file.path(temp_dir, paste0(flow_frame_name, ".xlsx"))
        saveWorkbook(wb, temp_file, overwrite = TRUE)
        

        temp_files <- c(temp_files, temp_file)
      }
      
      zip::zipr(zipfile = file, files = temp_files)
    }
  )
    
  
  
  
  
  

  observe({
    
    # If file is not clustered
    listObject$clusteringColumn <- "CLARA"

    if (input$alreadyClustered == TRUE) {
      
      listObject$marker_untrans <- extract_markers(listObject$flow.frames[[1]], NULL)

      updateSelectInput(session, "clusteringColumn", choices = listObject$marker_untrans)
    }
  })

  
  # VIZUALISATION
  
  observeEvent(input$run, {
    roots.output <- "/home/maelleWorkspace/OUTPUT/"
    roots.temp <- "/home/maelleWorkspace/TEMP/"
    
    # MultiOptSNE

    if ("MultiOptSNE" %in% input$step) {
      progress$set(message = "Write CSV", value = 0.1)

      write.csv(listObject$flow.frames.enriched@exprs[, unlist(input$marker_clustering)], paste0(roots.temp, rep, ".csv"), quote = FALSE, row.names = FALSE)
      progress$set(message = "Command line", value = 0.2)
      a <- paste0(roots.temp, rep, ".csv")
      b <- paste0(roots.temp, rep, "_tsne.csv")
      temp <- paste0(roots.temp, rep, ".log")
      iter <- 1000
      run <- "/media/data/cyto/MultiOptTSNE/Multicore-opt-SNE/MulticoreTSNE/run/run_optsne.py"
      cmd <- paste0(
        "python2 ", run, " --optsne --data ", a,
        " --outfile ", b, " --n_threads 35 --perp 30 --early_exaggeration 12 --n_iter ", iter, " ",
        "> ", temp, " 2>&1"
      )
      progress$set(message = "Run MultiOptTsne", value = 0.3)

      # Initialize output results
      system(cmd)
      roots.temp <- "/media/data/cyto/MultiOptTSNE/TEMP/"
      c <- paste0(roots.temp, rep, "_tsne.csv")
      mat <- read.csv(c, header = FALSE)
      # tester 2-3 essais
      if ("TSNE1" %in% colnames(fcs)) {
        fcs@exprs[, "TSNE1"] <- as.vector(unlist(mat[, 1]))
        
      } else {
        fcs <- FlowCIPHE::enrich.FCS.CIPHE(fcs, mat[, 1], "TSNE1")
      }
      if ("TSNE2" %in% colnames(fcs)) {
        fcs@exprs[, "TSNE2"] <- as.vector(unlist(mat[, 2]))
      } else {
        fcs <- FlowCIPHE::enrich.FCS.CIPHE(fcs, mat[, 2], "TSNE2")
      }
      values$data.concat <- fcs
      auto_save_r_data()
    }
  })
  
  observeEvent(input$downloadFCS, {

       if (!is.null(listObject$flow.frames.enriched)) {

         lapply(listObject$flow.frames.enriched, function(x){
         write.FCS(x,paste0(str_replace(names(unlist(listObject$flow.frames)),".fcs",""),"_corrected.fcs"))
          })}




  })

  
}
