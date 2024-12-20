##### Functions used by Cell Annotation App

#############################################
############  PREPROCESS ####################
#############################################


# Function that prepare the expression matrix with the good markers

preparMatrix<-function(flowFrame,markersInFCSFile, markersInModel, algo){
  # Extract expression matrix
  matrix<- as.data.frame(flowFrame@exprs)
  
  # Conditionnal loop for scyan algorithm
  if (!is.null(algo)){
    
  markersInFCSFile<-c(markersInFCSFile,"Time")
  
  markersInModel<-c(markersInModel,"Time")
  
  }else{
    
    markersInFCSFile<-c(markersInFCSFile)
    
    markersInModel<-c(markersInModel)
    
  }
  
  # Keep only markers you want
  df <- matrix[, markersInFCSFile]

  # Give features names of model to our expression matrix
  colnames(df)<-markersInModel
  
  # Return dataFrame
  return(df)
}


# Extract markers presents in a fcs file

extract_markers<-function(fcs_file, model){
  
  # Extract the part that contains the markers
  data<-as.data.frame(fcs_file@parameters@data)

  # Initialize marker vector
  markers<-c()
   
  # Initialize fluo vector 
  fluo <-c()
 
  # Loop 
  for (marker in 1:nrow(data)){
     
    # If  marker not contain description exemple for FSC-A/SSC-A
    if (is.na(data[marker,"desc"]) ){
      
      newMarker<-data[marker,"name"]
      
      newFluo<-data[marker,"name"]
     
    }else{
      
      newMarker<-data[marker,"desc"]
      
      newFluo<-data[marker,"name"]
    }
    
    markers<-c(markers,newMarker)
    
    fluo<-c(fluo,newFluo)
   
  }
 
  # Add markers 
  markers<-as.matrix(markers) 
  
  markers<-as.vector(markers[,1])
  
  fluo<-as.matrix(fluo)
  
  fluo<-as.vector(fluo[,1])
  
  names(fluo)<-markers
 
  return(fluo)
}


# MAIN function for preprocess 

pree_process_fcs <- function(flow.frames, arg, transformation, compens, markers_untrans){
    if(is.na(arg) || arg == ""){arg <- NULL}
    
    progress <- Progress$new()
    i <- 0
    
    # Apply preprocess on all fcs files
    
    flow.frames <- lapply(flow.frames, function(flow.frame){
      
      i <<- i + 1
     progress$set(message = paste0("Preprocessing ... ", i, "/", length(flow.frames), "."), value = i / length(as.vector(flow.frames)))
      
      # Compensation
      if (compens == TRUE) {
        
        # Found compensation matrix in the fcs file
        
        comp <- found.spill.CIPHE(flow.frames[[1]])
        
        # If there is a compensation matrix
        
        if(length(comp)>0){
        
          print("Found compensation matrix, applying...")
          
          comp <- flow.frame@description[comp][[1]]
       
          if(is.character(comp)){
            
            comp <- strsplit(comp, ",")[[1]]
            
            num.channels <- as.numeric(comp[1])
            
            # New matrix 
            m <- matrix(nrow = num.channels, byrow = T, data = as.numeric(comp[(num.channels + 2):length(comp)]))
            colnames(m) <- comp[2:(1 + num.channels)]
            comp <- m
          }
          
          # Apply compensation
          flow.frame <- compensate.CIPHE(flow.frame)
        }
      }

      
      if(transformation == "arcsinh"){
        
        # ARCSINH
        marker_trans <- colnames(flow.frame)[-which(colnames(flow.frame) %in% markers_untrans)]
        
        # Apply arcsinh function
        flow.frame <- arcsinh.CIPHE(flow.frame, markers_untrans, arg)
        
      } else if (transformation == "logicle") {
      
        # LOGICLE 
        
        flow.frame <- logicle.CIPHE(flow.frame, arg , markers_untrans) 
      }
      
      return(flow.frame)

    
   })
    return(flow.frames)
  }



# Logicle transformation 

logicle.CIPHE<-function (flow.frame, value = NULL, markers) 
  
{
  
  markers.transform <- markers
  #print(markers.transform)

  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description == x)))))

  list.index <- gsub("N", "", list.index)
  list.index <- gsub("\\$P", "", list.index)
  
  if (is.null(value)) {
    if (!is.null(flow.frame@description[[paste0("$P", list.index[1], 
                                                "MS")]])) {
      r.values <- unlist(lapply(list.index, function(x) as.integer(flow.frame@description[[paste0("$P", 
                                                                                                  x, "MS")]])))
    }
    else if (!is.null(flow.frame@description[[paste0("P", 
                                                     list.index[1], "MS")]])) {
      r.values <- unlist(lapply(list.index, function(x) as.integer(flow.frame@description[[paste0("P", 
                                                                                                  x, "MS")]])))
    }
    else {
      r.values <- rep(90, length(list.index))
    }
  }
  else {
    r.values <- rep(value, length(list.index))
  }
  w.values <- (4.5 - log10(262143/abs(r.values)))/2
  
  w.values[which(w.values < 0)] <- 0.5
  
  w.values[which(is.infinite(w.values))] <- 0.5
  for (t in 1:length(markers.transform)) {
    lgcl <- flowCore::logicleTransform(w = w.values[t])
    flow.frame <- flowCore::transform(flow.frame, transformList(markers.transform[t], 
                                                                lgcl))
  }
  return(flow.frame)
}



# Arcsinh transformation

arcsinh.CIPHE<- function (flow.frame, marker = NULL, arg) 
{
  raw <- flow.frame@exprs
  mat <- flow.frame@exprs
  if (is.null(marker) || length(marker) < 1) {
    marker <- c("FSC-A","SSC-A",colnames(flow.frame))
  }
  mat <- mat[, marker]
  colnames(mat) <- marker
  mat <- asinh(mat/arg)
  raw[, marker] <- mat[, marker]
  flow.frame@exprs <- raw
  return(flow.frame)
}



# Function to create a csv file 

builddCSVTab<-function(file, enrichName){ 
     
    # Extract colnames of expression matrix
    markers<- colnames(as.data.frame(file@exprs))
  
    flow.frame<-file
    
    # Extract paramters
    labels <- pData(flow.frame@parameters)[,2]
    
    params.names <- pData(parameters(flow.frame))[,1]
    
    # Give marker names to column that have NA 
    labels[which(is.na(labels))] <- colnames(flow.frame)[c(which(is.na(labels)))]
    
    labels[which(labels=="<NA>")] <- colnames(flow.frame)[c(which(labels=="<NA>"))]
    
    names(params.names) <- labels
    
    # over.clustering
    flow.frame.e <- flow.frame
    
    tab <- as.data.frame(flow.frame.e@exprs)
    
    data<-as.data.frame(file@parameters@data)
    
    fluo<-c()
    
    for (i in markers){
      
      for (marker in 1:nrow(data)){
        
        if (i %in% data[marker,"desc"]  | i %in% data[marker,"name"]){
          
          newFluo<-data[marker,"desc"]
          
        }else{
          newFluo<- i
        }
        
      }
      fluo<-c(fluo,newFluo)
    }
    colnames(tab) <- fluo
    
    tab.me<- ddply(.data = tab, .variables = enrichName, .fun = colwise(median, is.numeric))
    
    pop.size <- ddply(.data = tab, .variables = enrichName, .fun = nrow)
    
    names(pop.size) <- gsub("V1", "popsize", names(pop.size))
    
    res <- base::merge(tab.me, pop.size, by = enrichName)

  return(res)
}




############################################
############ ANNOTATION ####################
############################################


###### SCAFFOLD ANNOTATION

# CLARA CLUSTERING
claraClustering <- function(listFCS, markersClustering,k) {
  progress <- Progress$new()
  
  # Detect the operating system
  is_windows <- .Platform$OS.type == "windows"
  
  if (is_windows) {
    # Sequential execution for Windows
    message("Windows OS detected: Running CLARA clustering sequentially...")
    fcs_results <- lapply(seq_along(listFCS), function(i) {
      # Update progress
      progress$set(
        message = paste0("Processing file ", i, " of ", length(listFCS), "..."),
        value = i / length(listFCS)
      )
      
      # Prepare matrix to have only selected markers
      m <- listFCS[[i]]@exprs[,markersClustering]
    
      # Perform CLARA clustering
      cl <- cluster::clara(m, k, samples = 50)
      
      # Extract clustering results
      groups <- cl$clustering
      
      # Convert results to a matrix
      new_col <- as.matrix(groups)
      
      # Name the column that contains clustering results
      marker <- "CLARA"
      
      # Add column name
      colnames(new_col) <- marker
      
      # Add the new column with clustering results
      enrich.FCS.CIPHE(listFCS[[i]], new_col)
    })
  } else {
    # Parallel execution for non-Windows systems
    num_cores <- min(detectCores() - 1, 35) # Use up to 35 cores or available cores - 1
    message("Non-Windows OS detected: Running CLARA clustering in parallel...")
    
    fcs_results <- mclapply(seq_along(listFCS), mc.cores = num_cores, FUN = function(i) {
      # Prepare matrix to have only selected markers
      m <- listFCS[[i]]@exprs[,markersClustering]
      
      # Perform CLARA clustering
      cl <- cluster::clara(m, k, samples = 50)
      
      # Extract clustering results
      groups <- cl$clustering
      
      # Convert results to a matrix
      new_col <- as.matrix(groups)
      
      # Name the column that contains clustering results
      marker <- "CLARA"
      
      # Add column name
      colnames(new_col) <- marker
      
      # Add the new column with clustering results
      enrich.FCS.CIPHE(listFCS[[i]], new_col)
    })
  }
  
  progress$close()
  return(fcs_results)
}


# Function that build a scaffold map 

runn_analysis_gated<-function(flow.frames, clusteredFiles, map.clustedFiles.names, outputDir, 
         groups.clustering, col.names.gated, col.names.matrix, inter.cluster.connections, 
         col.names.inter_cluster, inter_cluster.weight_factor, overlap_method, ew_influence){
  

  print(sprintf("Markers used for SCAFFoLD: %s", paste(col.names.gated, collapse = ", ")))
  print(paste("Using as reference ", map.clustedFiles.names, sep = " "))
 
  gated_data <- load_attractors_from_gated_data(flow.frames, col.names.gated)
  
 
  tab.attractors <- gated_data$tab.attractors
  att.labels <- gated_data$cellType_key$population
  
  G.attractors <- NULL
  
  ret <- process_files(clusteredFiles, map.clustedFiles.names, G.attractors, tab.attractors, 
                       att.labels,col.names.gated=col.names.gated, col.names.matrix=col.names.matrix, 
                       scaffold.mode="gated", inter.cluster.connections=inter.cluster.connections, 
                       col.names.inter_cluster=col.names.inter_cluster
  )
  
  ret <- c(list(scaffold.col.names = col.names.gated, landmarks.data = gated_data$downsampled.data), ret)
  print(paste("Map file created:", sprintf("%s.scaffold", map.clustedFiles.names, sep = "/")))
  my_save(ret, paste0(outputDir, sprintf("%s.scaffold", map.clustedFiles.names), sep = ""))
  return(ret)
}



##########################################
########### VISUALIZATION ################
##########################################



# Function to match cell annotation with pop ID 
matchPopIDD<- function(){
  popLabel<-list()
  popLabel[[1]]<- "cDC2"
  popLabel[[2]]<- "CSTHi"
  popLabel[[3]]<- "CSTLo"
  popLabel[[4]]<- "DNDc"
  popLabel[[5]]<- "Eosinos"
  popLabel[[6]]<- "Granulos"
  popLabel[[7]]<- "Monos"
  popLabel[[8]]<- "NKcCD11b-"
  popLabel[[9]]<- "NKcCD11b+"
  popLabel[[10]]<- "NKTc"
  popLabel[[11]]<- "NI"
  popLabel[[12]]<- "pDCs"
  popLabel[[13]]<- "RPMac"
  popLabel[[14]]<- "TDN"
  popLabel[[15]]<- "TDP"
  popLabel[[16]] <- "TRegCM"
  popLabel[[17]]<- "B1B"
  popLabel[[18]]<-"B2B"
  popLabel[[19]]<- "CD4CM"
  popLabel[[20]]<- "CD4DN"
  popLabel[[21]]<-"CD4EM"
  popLabel[[22]]<- "CD4Naives"
  popLabel[[23]]<- "CD8CM"
  popLabel[[24]]<- "CD8DN"
  popLabel[[25]]<- "CD8EM"
  popLabel[[26]]<- "CD8Naives"
  popLabel[[27]]<-"cDC1"
  popLabel[[28]]<-"TRegEM"
  popLabel[[29]]<- "MZM"
  
  return(popLabel)
  
}



# Function to add a new column to a fcs file

enrich.FCS.CIPHE <- function(original, new.column)
{
  new_p <- parameters(original)[1,]
  
  ## Now, let's change it's name from $P1 to $P26 (or whatever the next new number is)
  new_p_number <- as.integer(dim(original)[2]+1)
  rownames(new_p) <- c(paste0("$P", new_p_number))
  
  ## Now, let's combine the original parameter with the new parameter
  
  library('BiocGenerics') ## for the combine function
  allPars <- BiocGenerics::combine(parameters(original), new_p)
  
  ## Fix the name and description of the newly added parameter, say we want to be calling it cluster_id
  new_p_name <- colnames(new.column)
  allPars@data$name[new_p_number] <- new_p_name
  allPars@data$desc[new_p_number] <- new_p_name
  
  # combine expression matrix 
  new_exprs <- cbind(original@exprs, new.column)
  
  new_kw <- original@description
  new_kw["$PAR"] <- as.character(new_p_number)
  new_kw[paste0("$P",as.character(new_p_number),"N")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"S")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"E")] <- "0,0"
  new_kw[paste0("$P",as.character(new_p_number),"G")] <- "1"
  new_kw[paste0("$P",as.character(new_p_number),"B")] <- new_kw["$P1B"]
  new_kw[paste0("$P",as.character(new_p_number),"R")] <- new_kw["$P1R"]
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- new_kw["flowCore_$P1Rmin"]
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- new_kw["flowCore_$P1Rmax"]
  
  ## Now, let's just combine it into a new flowFrame
  new_fcs <- new("flowFrame", exprs=new_exprs, parameters=allPars, description=new_kw)
  
  return(new_fcs)
}








