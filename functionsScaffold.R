##This file contains all the function used in the shiny application, grouped in sections.
##Many functions remain untouched from the original code found on this link:
##https://github.com/nolanlab/scaffold/tree/multiFilesClustering
##They are not commented.
##The new functions are on the top of their respective sections.

## OPTIONS
options(stringsAsFactors = F)

######## 0. GLOBAL FUNCTION ####################################################
################################################################################

my_load <- function(f_name){
  con <- file(f_name, "rb")
  retval <- unserialize(con)
  close(con)
  return(retval)
}


enrich.FCS.CIPHE <- function(original, new.column){
  
  new_p <- parameters(original)[1,]
  
  ## Now, let's change it's name from $P1 to $P26 (or whatever the next new number is)
  new_p_number <- as.integer(dim(original)[2]+1)
  rownames(new_p) <- c(paste0("$P", new_p_number))
  
  ## Now, let's combine the original parameter with the new parameter
  library('BiocGenerics') ## for the combine function
  allPars <- combine(parameters(original), new_p)
  
  ## Fix the name and description of the newly added parameter, say we want to be calling it cluster_id
  new_p_name <- colnames(new.column)
  allPars@data$name[new_p_number] <- new_p_name
  allPars@data$desc[new_p_number] <- new_p_name
  
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




write.FCS.CIPHE <- function(fcs, fcs.path, values.range = 262144, flow.data = T){
  if(is.null(flow.data)){
    if(is.null(target.fcs@description[["SPILL"]])){
      flow.data <- FALSE
    } else {
      flow.fata <- TRUE
    }
  }
  
  fcs.out <- flowFrame(fcs@exprs)
  
  descR <- description(fcs.out)
  lapply(1:ncol(fcs@exprs),function(x)
  {
    if(flow.data)
    {
      descR[[paste0("$P",x,"R")]] <<- values.range
    }
    else
    {
      descR[[paste0("$P",x,"R")]] <<- description(fcs)[[paste0("$P",x,"R")]]
    }
    descR[[paste0("$P",x,"S")]] <<- description(fcs)[[paste0("$P",x,"S")]]
  })
  pd <- pData(parameters(fcs))
  for(p in seq_along(pd[["name"]]))
  {
    descR[[sprintf("flowCore_$P%sRmax", p)]] <- pd[p,"minRange"]
    descR[[sprintf("flowCore_$P%sRmin", p)]] <- pd[p,"maxRange"]
  }
  levels(fcs.out@parameters@data[[2]]) <- fcs@parameters@data[[2]]
  fcs.out@parameters@data[[2]] <- fcs@parameters@data[[2]]
  if(flow.data)
  {
    fcs.out@parameters@data[[3]] <- values.range
    fcs.out@parameters@data[[4]] <- 0
    fcs.out@parameters@data[[5]] <- values.range-1
  }
  else
  {
    fcs.out@parameters@data[[3]] <- c(fcs@parameters@data[[3]], max(as.numeric(new.column)))
    fcs.out@parameters@data[[4]] <- c(fcs@parameters@data[[4]], 0)
    fcs.outs@parameters@data[[5]] <- c(fcs@parameters@data[[5]], max(as.numeric(new.column))-1)
  }
  
  fcs.out <- flowFrame(fcs@exprs, description = descR, parameters = fcs.out@parameters)
  if(flow.data)
  {
    if(!is.null(fcs@description[["$TIMESTEP"]]))
    {
      fcs.out@description[["$TIMESTEP"]] <- fcs@description[["$TIMESTEP"]]
    }
    if(!is.null(fcs@description[["APPLY COMPENSATION"]]))
    {
      fcs.out@description[["APPLY COMPENSATION"]] <- fcs@description[["APPLY COMPENSATION"]]
    }
  }
  if(!is.null(fcs@description[["SPILL"]]))
  {
    fcs.out@description[["SPILL"]] <- fcs@description[["SPILL"]]
  }
  
  write.FCS(fcs.out, fcs.path)
} 

updateFlowFrameKeywordsCIPHE <- function(flowFrame){
  
  row.names(flowFrame@parameters) <- paste0("$P",c(1:length(row.names(flowFrame@parameters))))
  params = parameters(flowFrame)
  pdata = pData(params)
  for (i in 1:ncol(flowFrame)){
    
    s = paste("$P",i,"S",sep="");
    n = paste("$P",i,"N",sep="");
    r = paste("$P",i,"R",sep="");
    b = paste("$P",i,"B",sep="");
    e = paste("$P",i,"E",sep="");
    fcmax1 <- paste("flowCore_$P",i,"Rmax",sep="");
    fcmin1 <- paste("flowCore_$P",i,"Rmin",sep="");
    fcmax <- paste("flowCore_P",i,"Rmax",sep="");
    fcmin <- paste("flowCore_P",i,"Rmin",sep="");
    display = paste0("P",i,"DISPLAY")
    bs = paste0("P",i,"BS")
    ms = paste0("P",i,"MS")
    
    keyval=list();
    label <- pData(flowFrame@parameters)[,"desc"][i]
    if(is.na(label)) {label <- colnames(flowFrame)[i] }
    keyval[[s]] = label
    keyval[[n]] = colnames(flowFrame)[i]         
    keyval[[r]] = ceiling(max(exprs(flowFrame)[,i])-min(exprs(flowFrame)[,i]))
    keyval[[b]] = 32;
    keyval[[e]] = "0,0";
    keyval[[fcmax1]] <- ceiling(max(exprs(flowFrame)[,i])-min(exprs(flowFrame)[,i]))
    keyval[[fcmin1]] <- ceiling(min(exprs(flowFrame)[,i]))
    keyval[[fcmax]] <- ceiling(max(exprs(flowFrame)[,i])-min(exprs(flowFrame)[,i]))
    keyval[[fcmin]] <- ceiling(min(exprs(flowFrame)[,i]))
    keyval[[bs]] <- 0
    keyval[[display]] <- "LOG"
    keyval[[ms]] <- 0
    
    keyword(flowFrame) = keyval;
    pdata[i,"minRange"]=min(exprs(flowFrame)[,i])
    pdata[i,"maxRange"]=max(exprs(flowFrame)[,i])
    
  }
  pData(params)=pdata
  parameters(flowFrame)=params
  row.names(flowFrame@parameters) <- paste0("$P",c(1:length(row.names(flowFrame@parameters))))
  return(flowFrame)
}

#Depending on the OS different functions are to choose a directory because no one is currently crossplatform.
chooseDir <-  function() {
  OS <- Sys.info()["sysname"]
  if (OS=="Windows") {
    Dir <-
      choose.dir(default = "", caption = "Select a Folder for saving:")
  }
  else if (OS=="Linux") {
    Dir <- tk_choose.dir(default = "", caption = "Select a Folder for saving:")
  }
  else {
    Dir <- choose.mac.dir()
  }
  return(Dir)
}

#Function used to select a folder via interface on a mac OS system.
choose.mac.dir <- function() {
  system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Select a Folder for saving:\")' > /tmp/R_folder",
         intern = FALSE, ignore.stderr = TRUE)
  p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  return(ifelse(length(p), p, NA))
}

get_cluster_groups_table <- function(v, key) {
  tags$table(class = "table table-hover table-striped",
             tags$tr(tags$td(
               v[1],
               tags$button(class = "btn btn-xs btn-warning pull-right", onClick = sprintf("Shiny.onInputChange('clusteringui_remove_clustering_group', {'key':'%s', 'x':Math.random()})", key),
                           tags$span(class = "glyphicon glyphicon-trash")
               )
             )),
             ifelse(length(v > 1),
                    tagList(lapply(tail(v, n = -1), function(x) {tags$tr(tags$td(x))})),
                    tagList()
             )
  )
}

returnOrder <- function(inputId, vars) {
  tagList(
    shiny::singleton(tags$head(tags$script(src = 'sort.js'))),
    shiny::singleton(tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'sort.css'))),
    shiny::HTML(html_list(vars, inputId)),
    tags$script(paste0("$(function() {$( '#",inputId,"' ).sortable({placeholder: 'ui-state-highlight'}); $( '#",inputId,"' ).disableSelection(); });"))
  )
}

updateReturnOrder <- function(session, inputId, vars){
  session$sendInputMessage(inputId, list(value = vars))
}

html_list <- function(vars, id) {
  hl <- paste0("<ul id=\'",id,"\' class='stab'>")
  for(i in vars) hl <- paste0(hl, "<li class='ui-state-default stab'><span class='label'>",i,"</span></li>")
  paste0(hl, "</ul>")
}

my_load <- function(f_name){
  con <- file(f_name, "rb")
  retval <- unserialize(con)
  close(con)
  return(retval)
}

my_save <- function(obj, f_name){
  con <- file(f_name, "wb")
  serialize(obj, con, ascii = F)
  close(con)
}

logiclTransformCIPHE <- function(flow.frame, value = NULL, markers = NULL){
  
  if(is.null(markers)){
    if(is.null(flow.frame@description[["SPILL"]])){
      markers.transform <- colnames(flow.frame)
    } else {
      markers.transform <- colnames(flow.frame@description[["SPILL"]])
    }
  } else {
    markers.transform <- markers
  }
  
  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)
  
  if(is.null(value)){
    if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]])){
      r.values <- unlist(lapply(list.index, function(x)
        as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
      )
    } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]])) {
      r.values <- unlist(lapply(list.index, function(x)
        as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
      )
    } else {
      r.values <- rep(90, length(list.index))
    }
  }
  else {
    r.values <- rep(value, length(list.index))
  }
  
  w.values <- (4.5-log10(262143/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5
  
  for(t in 1:length(markers.transform)){
    lgcl <- logicleTransform(w=w.values[t])
    flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
  }
  
  return(flow.frame)
}

inversLogiclTransformCIPHE <- function(flow.frame, value = NULL, markers = NULL){
  if(is.null(markers)){
    if(is.null(flow.frame@description[["SPILL"]])){
      markers.transform <- colnames(flow.frame)
    } else {
      markers.transform <- colnames(flow.frame@description[["SPILL"]])
    }
  } else {
    markers.transform <- markers
  }
  
  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)
  
  if(is.null(value)){
    if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]])) {
      r.values <- unlist(lapply(list.index, function(x) 
        as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
      ) 
    } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]])) {
      r.values <- unlist(lapply(list.index, function(x) 
        as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
      )   
    } else {
      r.values <- rep(90, length(list.index))
    }
  }
  else {
    r.values <- rep(value, length(list.index))
  }
  
  w.values <- (4.5-log10(262144/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5
  
  flow.frame.inv <- flow.frame
  
  for(t in 1:length(markers.transform)){
    invLgcl <- inverseLogicleTransform(trans = logicleTransform(w=w.values[t]))
    flow.frame.inv <- transform(flow.frame.inv, transformList(markers.transform[t],invLgcl))
  }
  
  return(flow.frame.inv)
}

arcsinhTransCIPHE <- function(flow.frame, marker=NULL, arg){
  raw <- flow.frame@exprs
  mat <- flow.frame@exprs
  if(is.null(marker) || length(marker)<1){
    marker <- colnames(flow.frame)
  }
  mat <- mat[,marker]
  colnames(mat) <- marker
  mat <- asinh(mat/arg)
  raw[,marker] <- mat[,marker]
  flow.frame@exprs <- raw
  return(flow.frame)
}

inversArcsinhTransCIPHE <- function(flow.frame, marker_untrans, arg){
  raw <- flow.frame@exprs
  mat <- flow.frame@exprs
  marker_untrans_index <- which(colnames(flow.frame)%in%marker_untrans)
  mat <- mat[,-marker_untrans_index]
  marker <- colnames(mat)
  mat <- sinh(mat)*arg
  raw[,marker] <- mat[,marker]
  flow.frame@exprs <- raw
  return(flow.frame)
}

deCompensateFlowFrame <- function(x, spillover) {
  if(!is.null(spillover)){
    cols <- colnames(spillover)
    sel <- cols %in% colnames(x)
    if(!all(sel)) {
      stop(keyword(x)[["FILENAME"]], "\\nThe following parameters in the spillover matrix are not present in the flowFrame:\\n",
           paste(cols[!sel], collapse=", "), call.=FALSE)
    }
    e <- exprs(x)
    e[, cols] <- e[, cols] %*% spillover
    exprs(x) = e
    return(x)
  } else {
    return(x)
  }
}

######## 1. CLUSTERING #########################################################
################################################################################

pre_process_fcs <- function(flow.frames, arg, transformation, compens, marker_untrans){
  flow.frames <- lapply(flow.frames, function(flow.frame){
    if (compens == TRUE) {
      comp <- found.spill.CIPHE(flow.frames[[1]])
      if(length(comp)>0){
        print("Found compensation matrix, applying...")
        comp <- description(flow.frame)[comp][[1]]
        if(is.character(comp)){
          comp <- strsplit(comp, ",")[[1]]
          num.channels <- as.numeric(comp[1])
          m <- matrix(nrow = num.channels, byrow = T, data = as.numeric(comp[(num.channels + 2):length(comp)]))
          colnames(m) <- comp[2:(1 + num.channels)]
          comp <- m
        }
        flow.frame <- compensate.CIPHE(flow.frame)
      }
    }
    
    #Any transformation method can be added here, but you must also modify server.R and ui.R accordingly.
    if(transformation == "Asinh"){
      marker_trans <- colnames(flow.frame)[-which(colnames(flow.frame)%in%marker_untrans)]
      flow.frame <- arcsinhTransCIPHE(flow.frame, marker_trans, arg)
    } else if(transformation == "Logicle"){
      flow.frame <- logiclTransformCIPHE(flow.frame, arg) # transform value
    }
    return(flow.frame)
  })
  return(flow.frames)
}

#Function used to transform data via Logicle transform, determining automaticly the better factor.
run_clustering <- function(flow.frames, methods, args, nb.cluster, params, 
                           outputDir, groups.clustering=NULL, ncores=1, transComp, marker_untrans,
                           method.matrix="median"){
  
  if(!is.null(groups.clustering) || length(groups.clustering)<length(flow.frames)){
    flow.frames <- lapply(groups.clustering, function(x){
      ff <- flow.frames[[x[1]]]
      file.params <- rep(1, dim(ff)[1])
      if(length(x)>1){
        lapply(c(2:length(x)), function(y){
          ff@exprs <<- rbind2(ff@exprs, cbind2(flow.frames[[x[y]]]@exprs), names(groups.clustering)[[y]])
          file.params <<- as.matrix(c(file.params, rep(y, dim(flow.frames[[x[y]]])[1])))
        })
      }
      file.params <- as.matrix(file.params)
      colnames(file.params) <- "sample"
      if("sample" %in% colnames(ff@exprs)){
        ff@exprs[, "sample"] <- as.matrix(file.params)
      } else {
        ff <- cbind2(ff, as.matrix(file.params))
      }
      return(ff)
    })
  } 
  
  marker.enrich <- paste0(methods, ".", dim(flow.frames[[1]]@exprs)[2]+1)
  
  cluster_flow_frame <- function(flow.frame, methods, outputDir, params, 
                                 nb.cluster, name, groups.clustering, method.matrix, args) {
    
    # source("ModifyFCS.R")
    library('BiocGenerics') ## for the combine function
    library("flowCore")
    library("flowAI")
    library("plyr")
    
    enrich.FCS.CIPHE <- function(original, new.column){
      
      new_p <- parameters(original)[1,]
      
      ## Now, let's change it's name from $P1 to $P26 (or whatever the next new number is)
      new_p_number <- as.integer(dim(original)[2]+1)
      rownames(new_p) <- c(paste0("$P", new_p_number))
      
      ## Now, let's combine the original parameter with the new parameter
      
      allPars <- combine(parameters(original), new_p)
      
      ## Fix the name and description of the newly added parameter, say we want to be calling it cluster_id
      new_p_name <- colnames(new.column)
      allPars@data$name[new_p_number] <- new_p_name
      allPars@data$desc[new_p_number] <- new_p_name
      
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
    
    m <- flow.frame@exprs[,params]
    
    labels <- pData(parameters(flow.frame))[,2]
    params.names <- pData(parameters(flow.frame))[,1]
    labels[which(is.na(labels))] <- colnames(flow.frame)[c(which(is.na(labels)))]
    names(params.names) <- labels
    
    # clustering
    if(methods == "CLARA"){
      groups <- cluster::clara(m,k=nb.cluster,samples=args)$clustering
    } else if(methods == "FlowSOM"){
      groups <- SOM(ReadInput(flow.frame[,params]),xgrid=args,ygrid=args)$map$mapping[,1]
    } else if(methods == "kmeans"){
      groups <- stats::kmeans(m, centers=nb.cluster,iter.max = 1000)$cluster
    }
    new_col <- as.matrix(groups)
    marker <- paste0(methods, ".", dim(flow.frame@exprs)[2]+1)
    colnames(new_col) <- marker
    
    flow.frame.e <- flow.frame
    flow.frame.e <- enrich.FCS.CIPHE(flow.frame, new_col)
    
    # tab <- flow.frame.e@exprs
    # colnames(tab) <- c(names(params.names), marker)
    # tab <- as.data.frame(tab)
    
    tab <- as.data.frame(flow.frame.e@exprs)
    names(tab) <- c(names(params.names), marker)
    
    if(method.matrix == "mean"){
      tab.me <- ddply(.data = tab, .variables = marker, .fun = colwise(mean, is.numeric))
    } else if(method.matrix == "median"){
      tab.me <- ddply(.data = tab, .variables = marker, .fun = colwise(median, is.numeric))
    }
    pop.size <- ddply(.data = tab, .variables = marker, .fun = nrow)
    names(pop.size) <- gsub("V1", "popsize", names(pop.size))
    res <- merge(tab.me, pop.size, by = marker)
    
    if(!is.null(groups.clustering)) {
      if(method.matrix == "mean"){
        tab.me.by.sample <- ddply(.data = tab, .variables = c(marker, "sample"), .fun = colwise(mean, is.numeric))
      } else if(method.matrix == "median"){
        tab.me.by.sample <- ddply(.data = tab, .variables = c(marker, "sample"), .fun = colwise(median, is.numeric))
      }
      
      pop.size.by.sample <- ddply(.data = tab, .variables = as.formula(paste0("~", marker, " * sample")), .fun = nrow)
      names(pop.size.by.sample) <- gsub("V1", "popsize", names(pop.size.by.sample))
      tab.me.by.sample <- merge(tab.me.by.sample, pop.size.by.sample, by = c(marker, "sample"), all.x = T)
      #Rotate the by.sample table
      temp <- reshape::melt(tab.me.by.sample, id = c(marker, "sample"))       
      temp$variable <- paste(temp$variable, temp$sample, sep = "@")
      temp$sample <- NULL
      formula <- paste(marker, "variable", sep = " ~ ")
      temp <- reshape::cast(temp, as.formula(formula))
      res <- merge(res, temp, by = marker, all.x = T)
    }
    
    return(list(flow.frame.e,res))
  }
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  flow.frames.e <- foreach(i = c(1:length(flow.frames))) %dopar% 
    cluster_flow_frame(flow.frames[[i]], methods, outputDir, params, nb.cluster, 
                       names(flow.frames)[[i]],groups.clustering, method.matrix, args)
  
  stopCluster(cl)
  
  # print(params)
  # print(length(flow.frames))
  # flow.frames.e <- mclapply(c(1:length(flow.frames)),mc.cores=ncores,FUN=function(i){
  #   cluster_flow_frame(flow.frames[[i]], methods, outputDir, params, nb.cluster=nb.cluster, 
  #                      names(flow.frames)[[i]], groups.clustering = groups.clustering,args=args)
  # })
  
  # i <- 1
  # flow.frames.e <- cluster_flow_frame(flow.frames[[i]], methods, outputDir, params, nb.cluster, 
  #   names(flow.frames)[[i]], groups.clustering, method.matrix, args)
  
  f1 <- lapply(flow.frames.e, function(a){return(a[[1]])})
  f2 <- lapply(flow.frames.e, function(b){return(b[[2]])})
  
  names(f1) <- names(flow.frames)
  names(f2) <- names(flow.frames)
  
  lapply(c(1:length(f2)), function(x) {
    lapply(c(1:length(groups.clustering[[x]])), function(y) {
      colnames(f2[[x]]) <<- gsub(paste0("@", y), paste0("@", groups.clustering[[x]][[y]]), colnames(f2[[x]]))
    })
  })
  
  lapply(c(1:length(f1)), function(x){
    outname <- paste0(outputDir,"/enriched_", names(f1)[x])
    ff <- f1[[x]]
    if (!is.na(transComp[4])) {
      if (transComp[2] == "Logicle") {
        print("Logicle detransformation...")
        ff <- inversLogiclTransformCIPHE(ff)
      } else if (transComp[2] == "Asinh") {
        ff <- inversArcsinhTransCIPHE(flow.frame = ff, marker_untrans = c(marker.enrich,marker_untrans,"sample"), arg = as.integer(transComp[3]))
      }
      if(transComp[1] == T)
        print("Decompensation.")
      ff <- deCompensateFlowFrame(ff, ff@description[["SPILL"]])
    }
    print("Writing FCS...")
    ff <- updateFlowFrameKeywordsCIPHE(ff)
    write.FCS(ff, outname)
  })
  
  return(list(f1,f2))
}

get_matrix_from_fcs <- function(i,
                                flow.frames, method.matrix="median", marker, groups.clustering=NULL) {
  
  flow.frame <- flow.frames[[i]]
  
  #avoid bugs when sample column missing
  if(!"sample"%in%colnames(flow.frame@exprs))
  {
    flow.frame@exprs <- cbind(flow.frame@exprs, sample=1)
    pData(flow.frame@parameters) <- rbind(pData(flow.frame@parameters), sample=c("sample", "sample", 1, 0, 2))
  }
  tab <- as.data.frame(flow.frame@exprs)
  
  labels <- pData(parameters(flow.frame))[,2]
  params.names <- pData(parameters(flow.frame))[,1]
  labels[which(is.na(labels))] <- colnames(flow.frame)[c(which(is.na(labels)))]
  labels[which(labels=="<NA>")] <- colnames(flow.frame)[c(which(labels=="<NA>"))]
  names(params.names) <- labels
  names(tab) <- names(params.names)
  
  if(method.matrix == "mean"){
    tab.me <- ddply(.data = tab, .variables = marker, .fun = colwise(mean, is.numeric))
  } else if(method.matrix == "median"){
    tab.me <- ddply(.data = tab, .variables = marker, .fun = colwise(median, is.numeric))
  }
  pop.size <- ddply(.data = tab, .variables = marker, .fun = nrow)
  names(pop.size) <- gsub("V1", "popsize", names(pop.size))
  res <- merge(tab.me, pop.size, by = marker)
  
  if (!is.null(groups.clustering)) {
    if(method.matrix == "mean"){
      tab.me.by.sample <- ddply(.data = tab, .variables = c(marker, "sample"), .fun = colwise(mean, is.numeric))
    } else if(method.matrix == "median"){
      tab.me.by.sample <- ddply(.data = tab, .variables = c(marker, "sample"), .fun = colwise(median, is.numeric))
    }
    pop.size.by.sample <- ddply(.data = tab, .variables = as.formula(paste0("~", marker, " * sample")), .fun = nrow)
    names(pop.size.by.sample) <- gsub("V1", "popsize", names(pop.size.by.sample))
    tab.me.by.sample <- merge(tab.me.by.sample, pop.size.by.sample, by = c(marker, "sample"), all.x = T)
    #Rotate the by.sample table
    temp <- reshape::melt(tab.me.by.sample, id = c(marker, "sample"))       
    temp$variable <- paste(temp$variable, temp$sample, sep = "@")
    temp$sample <- NULL
    formula <- paste(marker, "variable", sep = " ~ ")
    temp <- reshape::cast(temp, as.formula(formula))
    res <- merge(res, temp, by = marker, all.x = T)
  }
  
  lapply(c(1:length(groups.clustering[[i]])), function(y) {
    colnames(res) <<- gsub(paste0("@", y), paste0("@", groups.clustering[[i]][[y]]), colnames(res))
  })
  
  return(as.matrix(res))
}

######## 2. ANALYSIS ###########################################################
################################################################################

run_analysis_gated <- function(flow.frames, clusteredFiles, map.clustedFiles.names, outputDir, 
                               groups.clustering, col.names.gated, col.names.matrix, inter.cluster.connections, 
                               col.names.inter_cluster, inter_cluster.weight_factor, overlap_method, ew_influence){
  
  print(inter_cluster.weight_factor)
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

load_attractors_from_gated_data <- function(flow.frames, col.names, ...){
  res <- NULL
  for(i in c(1:length(flow.frames))) {

    population <- names(flow.frames)[i]
    
    fcs <- flow.frames[[i]]
    
    temp <- unlist(lapply(col.names,function(i){
      return(getChannelMarker(fcs,i)[1])
    }))
    
    tab <- exprs(fcs)[,temp]
    m <- as.matrix(tab)
    tab <- data.frame(m)
    names(tab) <- col.names
    
    tab <- as.matrix(tab)
    tab[tab < 0] <- 0
    tab <- as.data.frame(tab)
    
    popvector <- rep(population, dim(tab)[1])
    tab$population <- popvector
    res <- rbind(res, tab)
    
  }
  downsampled.data <- downsample_by(res, "population", 1000)
  names(downsampled.data) <- gsub("population", "cellType", names(downsampled.data))
  
  #Change cellType to be numbers
  k <- unique(res$population)
  k <- data.frame(population = k, cellType = seq_along(k), stringsAsFactors = F)
  res <- merge(res, k)
  res <- res[, grep("population", names(res), invert = T)]
  
  ## compute method for landmark with other like mean or mode
  res <- ddply(res, ~cellType, colwise(median))
  return(list(downsampled.data = downsampled.data, tab.attractors = res, cellType_key = k))
}

downsample_by <- function(tab, col.name, size){
  print(sprintf("Downsampling to %d events", size))
  return(ddply(tab, col.name, function(x, size) {
    if(nrow(x) <= size){
      return(x)
    } else {
      return(x[sample(1:nrow(x), size),])
    }
  }, size = size))
}

process_files <- function(clusteredFiles, map.clustedFiles.names=NULL, G.attractors, 
                          tab.attractors, att.labels, col.names.gated, col.names.matrix, scaffold.mode, 
                          ref.scaffold.markers = NULL, ew_influence = NULL,col.names.inter_cluster = NULL, ...){
  
  ret <- list(graphs = list(), clustered.data = list())
  
  names.mapping <- col.names.gated
  names(names.mapping) <- col.names.matrix
  map_names <- names_map_factory(names.mapping)
  
  if(is.null(map.clustedFiles.names)){
    map.clustedFiles.names <- names(clusteredFiles)[1]
  }
  
  print(paste("Processing ",map.clustedFiles.names, sep = " "))
  tab <- clusteredFiles[[map.clustedFiles.names]]
  
  if(""%in%colnames(tab))
  {
    index.null <- which(colnames(tab)=="")
    if (index.null[1] == 1){
      tab <- tab[,2:length(colnames(tab))]
      index.null <- which(colnames(tab)=="")
    }
    if(length(index.null) < 0) 
    {
      for (i in index.null){
        colnames(tab)[i] <- paste0("NA.", i)
      }
    }
  }
  
  index <- c()
  for(i in 1:length(col.names.matrix)){
    index <- c(index, which(col.names.matrix[i]==colnames(tab)))
  }
  
  colnames(tab)[index] <- col.names.gated   
  
  col.names.inter_cluster <- map_names(col.names.inter_cluster)
  
  if(is.null(ew_influence)) ew_influence <- ceiling(length(col.names.gated) / 3)
  
  #browser()
  tab <- tab[!apply(tab[, col.names.gated], 1, function(x) {all(x == 0)}),]
  names(tab) <- gsub("cellType", "groups", names(tab))
  names(tab) <- gsub("^X", "", names(tab))
  
  print(sprintf("Running with Edge weight: %f", ew_influence))
  res <- process_data(tab, map.clustedFiles.names, G.attractors, tab.attractors,
                      col.names.gated = col.names.gated, col.names.matrix = col.names.gated, 
                      att.labels = att.labels, already.clustered = T, ew_influence = ew_influence,
                      col.names.inter_cluster = col.names.inter_cluster)
  
  
  G.complete <- get_highest_scoring_edges(res$G.complete)
  
  ret$graphs[map.clustedFiles.names] <- list(G.complete)
  
  G.attractors <- res$G.attractors
  
  
  for(i in c(1:length(clusteredFiles))) {
    if(names(clusteredFiles)[i]==map.clustedFiles.names){
      NULL
    } else {
      print(paste("Processing", names(clusteredFiles)[[i]], sep = " "))
      tab <- clusteredFiles[[i]]
      
      if(""%in%colnames(tab))
      {
        index.null <- which(colnames(tab)=="")
        if (index.null[1] == 1){
          tab <- tab[,2:length(colnames(tab))]
          index.null <- which(colnames(tab)=="")
        }
        if(length(index.null) < 0) 
        {
          for (i in index.null){
            colnames(tab)[i] <- paste0("NA", i)
          }
        }
      }
      
      colnames(tab)[index] <- col.names.gated
      # colnames(tab) <- map_names(colnames(tab))
      # names(tab) <- colnames(tab)
      
      col.names.inter_cluster <- map_names(col.names.inter_cluster)
      
      if(scaffold.mode == "existing") {
        #Some markers in the reference scaffold file have been designated
        #for mapping, but they are missing from the sample files
        if(any(is.na(names(names.mapping))))
          tab <- add_missing_columns(tab, col.names, fill.data = 0)
        if(is.null(ew_influence))
          ew_influence <- ceiling(sum(!is.na(names(names.mapping))) / 3)
      } else {
        if(is.null(ew_influence)) ew_influence <- ceiling(length(col.names.gated) / 3)
      }
      
      tab <- tab[!apply(tab[, col.names.gated], 1, function(x) {all(x == 0)}),]
      names(tab) <- gsub("cellType", "groups", names(tab))
      names(tab) <- gsub("^X", "", names(tab))
      print(sprintf("Running with Edge weight: %f", ew_influence))
      res <- process_data(tab, map.clustedFiles.names, G.attractors, tab.attractors,
                          col.names.gated = col.names.gated, col.names.matrix = col.names.gated, 
                          att.labels = att.labels, already.clustered = T, ew_influence = ew_influence,
                          col.names.inter_cluster = col.names.inter_cluster, ...)
      G.complete <- get_highest_scoring_edges(res$G.complete)
      ret$graphs[names(clusteredFiles)[[i]]] <- list(G.complete)
      G.attractors <- res$G.attractors
    }
  }
  suppressWarnings( ret <- c(ret, list(dataset.statistics = get_dataset_statistics(ret))) )
  # ret <- c(ret, list(dataset.statistics = get_dataset_statistics(ret)))
  return(ret)
}



names_map_factory <- function(names.map){
  function(v)
  {
    sel <- v %in% names(names.map)
    if(any(sel))
      v[sel] <- names.map[v[sel]]
    return(v)
    
  }
}

get_highest_scoring_edges <- function(G) {
  
  E(G)$cluster_to_landmark <- 0
  E(G)$highest_scoring <- 0
  E(G)$inter_cluster <- 0
  
  e <- igraph::get.edges(G, E(G))
  E(G)$edge_type <- "cluster_to_landmark"
  e <- cbind(V(G)$type[e[,1]], V(G)$type[e[,2]])
  E(G)$edge_type[(e[,1] == 2) & (e[,2] == 2)] <- "inter_cluster"
  
  inter.cluster <- e[, 1] == 2 & e[, 2] == 2
  to.landmark <-(e[, 1] == 2 & e[, 2] == 1) |  (e[, 1] == 1 & e[, 2] == 2)
  
  E(G)$inter_cluster[inter.cluster] <- 1
  E(G)$cluster_to_landmark[to.landmark] <- 1
  
  # Remove inter-cluster edges for this calculation
  g.temp <- igraph::delete.edges(G, E(G)[inter.cluster])
  
  V(g.temp)$highest_scoring_edge <- 0
  
  for(i in 1:(igraph::vcount(g.temp))) {
    if(V(g.temp)$type[i] == 2) {
      sel.edges <- igraph::incident(g.temp, i)
      max.edge <- sel.edges[which.max(E(G)[sel.edges]$weight)]
      
      
      if(length(max.edge) != 0){
        V(g.temp)$highest_scoring_edge[i] <- max.edge
        E(G)$highest_scoring[max.edge] <- 1
        E(G)$edge_type[max.edge] <- "highest_scoring"
      }
      else {
        print("A vertex has no highest scoring edge...")
        V(g.temp)$highest_scoring_edge[i] <- 0
      }
    }
  }
  
  V(G)$highest_scoring_edge <- V(g.temp)$highest_scoring_edge
  return(G)
}

process_data <- function(tab, map.clustedFiles.names=NULL, G.attractors = NULL, 
                         tab.attractors = NULL, col.names.gated = NULL, col.names.matrix = NULL, att.labels = NULL, 
                         dist.thresh = 0.7,already.clustered = FALSE, inter.cluster.connections = FALSE, 
                         col.names.inter_cluster = NULL, inter_cluster.weight_factor = 0.7, ew_influence,
                         overlap_method = NULL){
  
  if(!already.clustered) {
    tab <- cluster_data(tab, col.names)
    tab.clustered <- ddply(tab, ~groups, colwise(median))
  }
  else
    tab.clustered <- tab
  
  if(is.null(col.names.inter_cluster) || col.names.inter_cluster == "")
    col.names.inter_cluster = col.names.gated
  if(is.null(G.attractors)) {
    G.attractors <- build_graph(tab.attractors, col.names.gated)
    
    G.complete <- add_vertices_to_attractors_graph(G.attractors, tab.clustered, tab.attractors, col.names.att = col.names.gated, col.names.matrix = col.names.matrix, dist.thresh)
    G.complete <- complete.forceatlas2(G.complete, first.iter = 50000,
                                       overlap.iter = 20000, ew_influence = ew_influence, overlap_method = "repel")
    if(inter.cluster.connections)
    {
      print("Adding inter-cluster connections with markers:")
      print(col.names.inter_cluster)
      print(sprintf("Weight factor:%f", inter_cluster.weight_factor))
      G.complete <- add_inter_clusters_connections(G.complete, col.names.inter_cluster, weight.factor = inter_cluster.weight_factor)
      G.complete <- complete.forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
                                         ew_influence = ew_influence, overlap_method = overlap_method)
    }
    V(G.attractors)$x <- V(G.complete)$x[1:vcount(G.attractors)]
    V(G.attractors)$y <- V(G.complete)$y[1:vcount(G.attractors)]
  } else {
    
    G.complete <- add_vertices_to_attractors_graph(G.attractors, tab.clustered, tab.attractors, col.names.att = col.names.gated, col.names.matrix = col.names.matrix, dist.thresh)
    
    fixed <- rep(FALSE, vcount(G.complete))
    fixed[1:vcount(G.attractors)] <- TRUE
    
    G.complete <- complete.forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
                                       overlap_method = "repel", ew_influence = ew_influence, fixed = fixed)
    if(inter.cluster.connections) {
      print("Adding inter-cluster connections with markers:")
      print(col.names.inter_cluster)
      print(sprintf("Weight factor:%f", inter_cluster.weight_factor))
      G.complete <- add_inter_clusters_connections(G.complete, col.names.inter_cluster, weight.factor = inter_cluster.weight_factor)
      G.complete <- complete.forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
                                         overlap_method = overlap_method, ew_influence = ew_influence, fixed = fixed)
    }
    
  }
  
  G.complete <- add_attractors_labels(G.complete, att.labels)
  V(G.complete)$name <- gsub(".fcs", "", V(G.complete)$name)
  
  # G.complete <- get_highest_scoring_edges(G.complete)
  # V(G.complete)$celltype <- att.labels
  # G.complete <- add_landmarks_labels(G.complete)
  # V(G.complete)$name <- gsub(".fcs", "", G.complete)
  # print(V(G.complete)$name)
  # return(list(G.attractors = G.attractors, G.complete = G.complete))
  
  return(list(G.attractors = G.attractors, G.complete = G.complete, tab.attractors = tab.attractors, tab = tab, col.names.gated = col.names.gated, col.names.matrix = col.names.matrix))
}

add_landmarks_labels <- function(G, v) {
  w <- V(G)$type == 1
  V(G)$name[w] <- V(G)$Label[w] <- V(G)$cellType[w]
  return(G)
}

add_vertices_to_attractors_graph <- function(G, tab.clustered, tab.median, col.names.att, col.names.matrix, dist.thresh = 0.7){
  dd <- get_distances_from_attractors(tab.clustered, tab.median, col.names.att, col.names.matrix, dist.thresh)
  n <- nrow(dd)
  num.vertices <- length(V(G))
  G <- add.vertices(G, n)
  v.seq <- (num.vertices + 1):length(V(G))
  V(G)[v.seq]$name <- as.character(v.seq)
  V(G)[v.seq]$Label <- paste("c", tab.clustered[,1], sep = "")
  row.names(dd) <- as.character(v.seq)
  for(i in 1:nrow(dd))
  {
    v <- dd[i,]
    v <- v[v > 0]
    if(length(v) > 0)
    {
      e.list <- c(rbind(as.character(num.vertices + i), names(v)))
      G <- G + edges(e.list, weight = v)
    }
  }
  
  # weight <- E(G)$weight
  # E(G)$weight <- weight ^ 10
  maxx <- maxy <- rep(Inf, vcount(G))
  minx <- miny <- rep(-Inf, vcount(G))
  
  maxx[1:num.vertices] <- minx[1:num.vertices] <- V(G)$x[1:num.vertices]
  maxy[1:num.vertices] <- miny[1:num.vertices] <- V(G)$y[1:num.vertices]
  lay <- layout.kamada.kawai(G, minx = minx, maxx = maxx, miny = miny, maxy = maxy)
  colnames(lay) <- c("x", "y")
  G <- set.vertex.attribute(G, name = "x", value = lay[, "x"])
  G <- set.vertex.attribute(G, name = "y", value = lay[, "y"])
  
  V(G)[1:num.vertices]$type <- 1 #attractor
  V(G)[(num.vertices + 1):vcount(G)]$type <- 2 #cell
  
  for(i in colnames(tab.clustered))
  {
    G <- set.vertex.attribute(G, name = i, index = (num.vertices + 1):vcount(G), value = tab.clustered[, i])
  }
  
  
  
  G <- set_visual_attributes(G)
  return(G)
}

set_visual_attributes <- function(G){
  att <- V(G)$type == 1
  V(G)$r <- 79
  V(G)$g <- 147
  V(G)$b <- 222
  V(G)$size <- 10
  
  V(G)[att]$r <- 255
  V(G)[att]$g <- 117
  V(G)[att]$b <- 128
  V(G)[att]$size <- 20
  
  E(G)$r <- 180
  E(G)$g <- 180
  E(G)$b <- 180
  
  return(G)
}

build_graph <- function(tab, col.names, filtering_T = 0.7, method = "cosine"){
  # print(tab)
  print(tab$cellTyp)
  m <- as.matrix(tab[, col.names])
  row.names(m) <- tab$cellTyp
  #Any distance metric can be added here
  if (method == "cosine") {
    dd <- cosine_similarity_matrix(m)
    diag(dd) <- 0
    dd[is.na(dd)] <- 0 #This can happen if one of the attractors has all 0's for the markers of interest
    if(filtering_T >= 1) {
      dd <- filter_similarity_matrix_by_rank(dd, filtering_T)
    }
    else {
      dd <- filter_similarity_matrix(dd, filtering_T)
    }
  }
  else if (method == "euclidean") {
    dd <- euclid_similarity_matrix(m)
    diag(dd) <- 0
    dd[is.na(dd)] <- 0
  }
  G <- graph.adjacency(dd, mode = "undirected", weighted = T)
  n.vertices <- length(V(G))
  lay <- layout.kamada.kawai(G)
  colnames(lay) <- c("x", "y")
  G <- set.vertex.attribute(G, name = "x", value = lay[, "x"])
  G <- set.vertex.attribute(G, name = "y", value = lay[, "y"])
  for(i in names(tab))
    G <- set.vertex.attribute(G, name = i, value = tab[, i])
  return(G)
}

cosine_similarity_matrix <- function(m){
  ret <- t(apply(m, 1, function(x, m) {cosine_similarity_from_matrix(x, m)}, m = m))
  return(ret)
}

cosine_similarity_from_matrix <- function(v, m){
  m <- as.matrix(m[, c(1:length(v)), drop = F])
  ret <- apply(m, 1, function(x, v) {
    return(crossprod(x, v)/sqrt(crossprod(x) * crossprod(v)))}, v)
  return(ret)
}

euclid_similarity_matrix <- function(m){
  euclid <- as.matrix(dist(m, method="euclidean"))
  # ret <- 1 - (euclid/max(euclid)) #simple similarity
  ret  <- exp(-(euclid/(0.4*4.5))^2) #more complex similarity
  return(ret)
}

filter_similarity_matrix <- function(m, T){
  ret <- t(apply(m, 1, function(x)
  {
    if(max(x) <= T)
      x[x < max(x)] <- 0
    else
      x[x < T] <- 0
    return(x)
  }))
  return(ret)
}

filter_similarity_matrix_by_rank <- function(m, T){
  ret <- t(apply(m, 1, function(x)
  {
    r <- rank(x, ties.method = "first")
    r <- max(r) - r + 1
    x[r > T] <- 0
    return(x)
  }))
  return(ret)
}

get_distances_from_attractors <- function(m, tab, col.names.att, col.names.matrix, dist.thresh){
  att <- as.matrix(tab[, col.names.att])
  row.names(att) <- as.character(1:nrow(tab))
  m <- as.matrix(m[, col.names.att])
  dd <- t(apply(m, 1, function(x, att) {cosine_similarity_from_matrix(x, att)}, att))
  dd <- distance_from_attractor_hard_filter(dd, tab, col.names.att, thresh = 1) ##marque page
  dist.thresh <- quantile(dd, probs = 0.85, na.rm = T)
  dist.thresh <- max(c(dist.thresh, 0.5))
  dd[is.na(dd)] <- 0 #This can happen if one of the attractors has all 0's for the markers of interest
  dd <- filter_similarity_matrix(dd, dist.thresh)
  return(dd)
}

distance_from_attractor_hard_filter <- function(dd, tab, col.names.att, thresh = 0.5){
  tab <- tab[, col.names.att]
  w <- apply(tab[,col.names.att], 1, function(x, thresh) {all(x < thresh)}, thresh = thresh)
  if(any(w))
    print("Hard removing some connections to unstained landmarks")
  dd[,w] <- 0
  return(dd)
}

complete.forceatlas2 <- function(G, first.iter = 1000, overlap.iter, overlap_method = NULL, ...){
  
  print("First iteration")
  ret <- layout.forceatlas2(G, prevent.overlap = FALSE, iter = first.iter, ...)
  
  lay <- ret$lay
  
  G <- set.vertex.attribute(G, name = "x", value = lay[, 1])
  G <- set.vertex.attribute(G, name = "y", value = lay[, 2])
  if(!is.null(overlap_method))
  {
    if(overlap_method == "repel")
    {
      print("Second iteration with prevent overalp")
      ret <- layout.forceatlas2(G, prevent.overlap = TRUE, iter = overlap.iter, ...)
      lay <- ret$lay
      if(any(is.na(lay)))
      {
        print("Prevent overlap iteration failed")
      }
      
      else
      {
        G <- set.vertex.attribute(G, name = "x", value = lay[, 1])
        G <- set.vertex.attribute(G, name = "y", value = lay[, 2])
      }
    }
    else if(overlap_method == "expand")
      G <- adaptive_expand(G, overlap.iter)
  }
  return(G)
}

# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

layout_forceatlas2Cpp <- function(lay, F_att_orig, mass, nodes_size, edge_list, avg_displ, kgrav, iter, prevent_overlap, fixed, max_displ, stopping_tolerance, barnes_hut) {
  layout_forceatlas2Cpp(lay, F_att_orig, mass, nodes_size, edge_list, avg_displ, kgrav, iter, prevent_overlap, fixed, max_displ, stopping_tolerance, barnes_hut)
}

layout.forceatlas2 <- function(G, ew_influence = 1, kgrav = 1, iter = 1000, prevent.overlap = FALSE, fixed = rep(FALSE, vcount(G)), stopping_tolerance = 0.001, barnes_hut = FALSE){
  if(vcount(G) >= 2000)
    barnes_hut <- TRUE
  if(vcount(G) > 2000)
    stopping_tolerance <- 0.01
  else if(vcount(G) > 800)
    stopping_tolerance <- 0.005
  else
    stopping_tolerance <- 0.001
  
  if(is.null(get.vertex.attribute(G, "x")))
  {
    lay <- cbind(x = rnorm(vcount(G)), y = rnorm(vcount(G)))
  }
  else
  {
    lay <- cbind(x = V(G)$x, y = V(G)$y)
  }
  
  #This is only used with prevent.overlap
  if(is.null(get.vertex.attribute(G, "size")))
    V(G)$size <- rep(10, vcount(G))
  mass <- 1 + degree(G)
  F_att <- (E(G)$weight ^ ew_influence)
  edge_list <- get.edgelist(G, names = F) - 1 # This is gonna be used in the C code where the indexing is 0-based
  
  avg_displ <- numeric(iter)
  max_displ <- numeric(iter)
  print(system.time(layout_forceatlas2Cpp(lay, F_att, mass, V(G)$size, edge_list, avg_displ,
                                          kgrav,  iter, prevent.overlap, fixed, max_displ, stopping_tolerance, barnes_hut)))
  return(list(lay = lay, avg_displ = avg_displ, max_displ = max_displ))
}

add_inter_clusters_connections <- function(G, col.names, weight.factor, method = "cosine"){
  tab <- get_vertex_table(G)
  tab <- tab[tab$type == 2,]
  
  m <- as.matrix(tab[, col.names])
  row.names(m) <- tab$name
  #Any distance metric can be added here
  if (method == "cosine") {
    dd <- cosine_similarity_matrix(m)
    diag(dd) <- 0
    dd[is.na(dd)] <- 0 #This can happen if one of the attractors has all 0's for the markers of interest
    dist.thresh <- quantile(dd, probs = 0.85, na.rm = T)
    dist.thresh <- max(c(dist.thresh, 0.5))
    dd <- filter_similarity_matrix(dd, dist.thresh)
    dd <- filter_similarity_matrix_by_rank(dd, 3)
  }
  else if (method == "euclidean") {
    dd <- euclid_similarity_matrix(m)
    diag(dd) <- 0
    dd[is.na(dd)] <- 0
    dist.thresh <- quantile(dd, probs = 0.7, na.rm = T)
    dist.thresh <- max(c(dist.thresh, 0.5))
    dd <- filter_similarity_matrix(dd, dist.thresh)
    dd <- filter_similarity_matrix_by_rank(dd, 3)
  }
  
  e.list <- NULL
  
  for(i in 1:nrow(dd))
  {
    v <- dd[i,]
    v <- v[v > 0]
    if(length(v) > 0)
    {
      e.list <- rbind(e.list, data.frame(a = tab[i, "name"], b = names(v), weight = v, stringsAsFactors = FALSE))
      
    }
  }
  temp <- as.matrix(e.list[, c("a", "b")])
  temp <- t(apply(temp, 1, sort))
  e.list <- data.frame(temp, weight = e.list$weight, stringsAsFactors = FALSE)
  names(e.list)[1:2] <- c("a", "b")
  e.list <- e.list[!duplicated(e.list[, c("a", "b")]),]
  e.list.igraph <- c(t(as.matrix(e.list[, c("a", "b")])))
  
  #G <- G + edges(e.list, weight = (v ^ 30))
  G <- G + edges(e.list.igraph, weight = e.list$weight * weight.factor)
  return(G)
}

adaptive_expand <- function(G, max.iter){
  print("Starting adaptive expansion")
  x <- V(G)$x
  y <- V(G)$y
  m <- cbind(x, y)
  ss <- outer(V(G)$size, V(G)$size, "+")
  
  for(i in 1:max.iter)
  {
    dd <- as.matrix(dist(m), method = "euclidean")
    dd <- dd - ss
    dd <- dd[upper.tri(dd)]
    if(all(dd >= 0))
      break
    else
      m <- m * 1.2
  }
  
  print(sprintf("Expansion stopped at iteration: %d", i))
  V(G)$x <- m[, "x"]
  V(G)$y <- m[, "y"]
  
  return(G)
}

get_vertex_table <- function(G){
  att <- list.vertex.attributes(G)
  ret <- NULL
  
  for(a in att)
  {
    d <- data.frame(get.vertex.attribute(G, a), stringsAsFactors = FALSE)
    if(is.null(ret))
      ret <- d
    else
      ret <- cbind(ret, d, stringsAsFactors = FALSE)
  }
  names(ret) <- att
  return(ret)
}

add_attractors_labels <- function(G, v){
  V(G)$name[1:length(v)] <- V(G)$Label[1:length(v)] <- v
  return(G)
}

get_dataset_statistics <- function(dataset){
  graphs <- dataset$graphs
  ret <- NULL
  for(G in graphs)
  {
    V(G)$popsize.relative <- V(G)$popsize / sum(V(G)$popsize, na.rm = T)
    tab <- get.data.frame(G, what = "vertices")
    max.vals <- sapply(tab, function(x) {if(is.numeric(x)) return(max(x, na.rm  =T))})
    max.vals <- max.vals[!sapply(max.vals, is.null)]
    for(i in 1:length(max.vals))
    {
      var.name <- names(max.vals)[i]
      if(!(var.name %in% names(ret)) || max.vals[[i]] > ret[[var.name]])
        ret[var.name] <- max.vals[i]
    }
  }
  return(list(max.marker.vals = ret))
}

######## 3. Map ################################################################
################################################################################

reactiveNetwork <- function (outputId)
{
  HTML(paste("<div id=\"", outputId, "\" class=\"shiny-network-output\"><svg /></div>", sep=""))
}

get_numeric_vertex_attributes <- function(sc.data, sel.graph){
  G <- sc.data$graphs[[sel.graph]]
  d <- get.data.frame(G, what = "vertices")
  #Don't consider attributes which are only present in the landmarks
  d <- d[d$type == 2,]
  num <- sapply(d, function(x) {is.numeric(x) && !any(is.na(x))})
  v <- list.vertex.attributes(G)[num]
  v <- v[grep("@", v, invert = T)]
  exclude <- c("x", "y", "cellType", "type", "groups", "r", "g", "b", "size", "DNA1", "DNA2", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "Time", "Cell_length", "Cisplatin", "beadDist", "highest_scoring_edge")
  return(v[!(v %in% exclude)])
}

combine_marker_sample_name <- function(sel.marker, active.sample){
  if(active.sample == "All" || active.sample == "Absolute" || sel.marker == "Default")
    return(sel.marker)
  else
    return(paste(sel.marker, active.sample, sep = "@"))
}

get_sample_names <- function(sc.data, sel.graph){
  G <- sc.data$graphs[[sel.graph]]
  s <- list.vertex.attributes(G)
  if (sel.graph != names(sc.data$graphs)[1]) {
    
    ids <- which(list.vertex.attributes(sc.data$graphs[[1]])%in%s)
    s <- lapply(1:length(s), function(i) {
      if (i%in%ids) {
        return(NULL)
      } else {return(s[[i]])}
    })
  }
  s <- grep("@", s, value = T)
  ret <- sapply(strsplit(s, "@"), function (x) {x[[2]]})
  if (length(ret) == 0) {
    ret <- list(1)
  }
  return(unique(ret))
}

get_graph <- function(sc.data, sel.graph, node.size.attr, min.node.size, max.node.size, landmark.node.size){
  G <- sc.data$graphs[[sel.graph]]
  edges <- data.frame(get.edgelist(G, names = F) - 1)
  colnames(edges) <- c("source", "target")
  svg.width <- 1200
  svg.height <- 800
  svg.center <- c(svg.width/2, svg.height/2)
  
  x <- V(G)$x
  y <- V(G)$y
  
  y <- -1 * y
  x <- x + abs(min(x))
  y <- y + abs(min(y))
  num.landmarks <- sum(V(G)$type == 1)
  trans <- get_graph_centering_transform(x[V(G)$type == 1], y[V(G)$type == 1], svg.width, svg.height)
  
  x <- (x / trans$scaling) - trans$offset.x
  y <- (y / trans$scaling) - trans$offset.y
  
  vertex.size <- get_vertex_size(sc.data, sel.graph, svg.width, node.size.attr, min.node.size, max.node.size, landmark.node.size)
  edges <- cbind(edges, x1 = x[edges[, "source"] + 1], x2 = x[edges[, "target"] + 1])
  edges <- cbind(edges, y1 = y[edges[, "source"] + 1], y2 = y[edges[, "target"] + 1])
  edges <- cbind(edges, id = 1:nrow(edges))
  edges <- cbind(edges, is_highest_scoring = 0)
  edges <- cbind(edges, edge_type = "")
  #Set as true for the highest scoring edges of type 2 vertices
  edges[, "is_highest_scoring"][V(G)$highest_scoring_edge[V(G)$type == 2]] <- 1
  if("edge_type" %in% list.edge.attributes(G)) #Old graphs did not have this
    edges[, "edge_type"] <- E(G)$edge_type
  ret <- list(names = V(G)$Label, size = vertex.size / trans$scaling, type = V(G)$type, highest_scoring_edge = V(G)$highest_scoring_edge, X = x, Y = y)
  ret <- c(ret, edges = list(edges))
  return(ret)
}

get_color_for_marker <- function(
    sc.data, sel.marker, rel.to.sample, sel.graph, active.sample, color.scaling,
    stats.type, colors.to.interpolate, color.under, color.over, color.scale.limits = NULL, 
    color.scale.mid = NULL)
{
  G <- sc.data$graphs[[sel.graph]]
  if(sel.marker == "Default") {
    ret <- rep("#4F93DE", vcount(G))
    ret[V(G)$type == 1] <- "#FF7580"
    return(list(color.vector = ret, color.scale.lim = NULL))
  } else {
    v <- get.vertex.attribute(G, combine_marker_sample_name(sel.marker, active.sample))
    f <- colorRamp(colors.to.interpolate, interpolate = "linear")
    if(rel.to.sample != "Absolute") {
      rel.to.marker <- combine_marker_sample_name(sel.marker, rel.to.sample)
      if(stats.type == "Difference")
        v <- v - (get.vertex.attribute(G, rel.to.marker))
      else if(stats.type == "Ratio")
        v <- v / (get.vertex.attribute(G, rel.to.marker))
      v[is.infinite(v)] <- NA
    }
    color.scale.lim <- NULL
    if(color.scaling == "local") color.scale.lim <- list(min = min(v, na.rm = T), max = max(v, na.rm = T))
    if(!is.null(color.scale.limits)) {
      under <- v < color.scale.limits[1]
      over <- v > color.scale.limits[2]
      v[under] <- color.scale.limits[1]
      v[over] <- color.scale.limits[2]
      if(is.null(color.scale.mid))
        v <-  scales::rescale(v)
      else
        v <- scales::rescale_mid(v, mid = color.scale.mid)
      v <- f(v)
      v <- apply(v, 1, function(x) {sprintf("rgb(%s)", paste(round(x), collapse = ","))})
      v[under] <- sprintf("rgb(%s)", paste(col2rgb(color.under), collapse = ","))
      v[over] <- sprintf("rgb(%s)", paste(col2rgb(color.over), collapse = ","))
    } else {
      v <- f(scales::rescale(v)) #colorRamp needs an argument in the range [0, 1]
      v <- apply(v, 1, function(x) {sprintf("rgb(%s)", paste(round(x), collapse = ","))})
    }
    return(list(color.vector = v, color.scale.lim = color.scale.lim))
  }
}

get_number_of_cells_per_landmark <- function(sc.data, sel.graph){
  G <- sc.data$graphs[[sel.graph]]
  land <- V(G)[V(G)$type == 1]$Label
  ee <- get.edgelist(G)
  ee <- ee[V(G)[V(G)$type == 2]$highest_scoring_edge,]
  vv <- V(G)[as.numeric(ee[,2])]
  
  popsize <- V(G)[vv]$popsize
  dd <- data.frame(Landmark = ee[,1], popsize)
  dd <- ddply(dd, ~Landmark, function(x) {sum(x["popsize"])})
  dd <- cbind(dd, Percentage = dd$V1 / sum(dd$V1))
  names(dd) <- c("Landmark", "Cells", "Percentage")
  dd$Percentage <- signif(dd$Percentage * 100, digits = 4)
  return(dd)
}

get_cells_per_landmark_all_files <- function (sc.data){
  cells <- lapply(c(1:length(sc.data$graphs)), function(i) {
    name <- names(sc.data$graphs)[i]
    cells <- get_number_of_cells_per_landmark(sc.data, name)
    return(cells)
  })
  m <- matrix(nrow = length(sc.data$graphs), ncol = nrow(cells[[1]])*2)
  rownames(m) <- names(sc.data$graphs)
  colnames(m) <- sapply(c(1:ncol(m)), function(i) {
    name <- as.vector(cells[[1]][,1])[ceiling(i/2)]
    if (gtools::odd(i) == TRUE) {name <- paste0(name, " #events")}
    else {name <- paste0(name, " %percentage")}
    return(name)
  })
  i <- 0
  j <- 0
  lapply(c(1:length(cells)), function(i) {
    lapply(c(1:nrow(cells[[i]])), function(j) {
      match <- pmatch(cells[[i]][j, ][[1]], as.vector(cells[[1]]$Landmark))
      m[i, (match*2)-1] <<- cells[[i]][j, ][[2]]
      m[i, (match*2)] <<- cells[[i]][j, ][[3]]
    })
  })
  m[is.na(m)] <- 0
  return(m)
}

get_summary_table <- function(sc.data, sel.graph, sel.nodes){
  G <- sc.data$graphs[[sel.graph]]
  col.names <- get_numeric_vertex_attributes(sc.data, sel.graph)
  tab <- get.data.frame(G, what = "vertices")
  temp <-tab[tab$Label %in% sel.nodes,]
  ret <- temp[, col.names]
  ret <- rbind(ret, apply(ret, 2, median, na.rm = T))
  popsize <- data.frame(Cells = temp$popsize, Percentage = temp$popsize / sum(tab$popsize[tab$type == 2]))
  popsize <- rbind(popsize, colSums(popsize))
  ret <- cbind(popsize, ret)
  ret <- data.frame(Label = c(temp$Label, "Summary"), ret)
  ret$Percentage <- signif(ret$Percentage * 100, digits = 4)
  return(ret)
}

plot_cluster <- function(data, clusters, graph.name, col.names, pool.cluster.data, plot.type){
  G <- data$graphs[[graph.name]]
  gated_data <- data$landmarks.data
  clustered_data <- data$clustered.data[[graph.name]]
  
  names(clustered_data) <- gsub("^X", "", names(clustered_data))
  names(gated_data) <- gsub("^X", "", names(gated_data))
  
  #This only works if the col.names are actually present in the clustered.data
  #TODO: figure out a consistent way to deal with panel mismatches
  common.names <- col.names[(col.names %in% names(clustered_data)) & (col.names %in% names(gated_data))]
  clustered_data <- clustered_data[, c(col.names, "cellType")]
  gated_data <- gated_data[, c(common.names, "cellType")]
  gated_data <- add_missing_columns(gated_data, col.names, fill.data = NA)
  #Select only the landmark nodes that are connected to these clusters
  land <- V(G)[nei(V(G)$Label %in% clusters)]$Label
  land <- V(G)[(V(G)$Label %in% land) & V(G)$type == 1]$Label
  temp <- gated_data[gated_data$cellType %in% land,]
  clus.num <- as.numeric(gsub("c", "", clusters))
  temp.clustered <- clustered_data[clustered_data$cellType %in% clus.num, ]
  if(pool.cluster.data)
    temp.clustered$cellType <- "Clusters"
  temp <- rbind(temp, temp.clustered)
  p <- NULL
  if(plot.type == "Scatterplot")
  {
    p <- density_scatterplot(temp, x_name = col.names[1], y_name = col.names[2], grouping = "cellType")
  }
  else
  {
    temp <- melt(temp, id.vars = "cellType")
    temp$variable <- as.factor(temp$variable)
    if(plot.type == "Density")
    {
      p <- ggplot(aes(x = value, color = cellType), data = temp) + geom_density() + facet_wrap(~variable, scales = "free")
    }
    else if(plot.type == "Boxplot")
    {
      p <- ggplot(aes(x = variable, fill = cellType, y = value), data = temp) + geom_boxplot()
    }
  }
  plot(p)
  return(p)
}

density_scatterplot  <- function(tab, x_name, y_name, grouping){
  m <- ddply(tab, grouping, function(m, x_name, y_name)
  {
    colramp <- grDevices::colorRampPalette(c("black", "red", "yellow"))
    dens.col <- grDevices::densCols(m[, x_name], m[, y_name], colramp = colramp)
    return(data.frame(m, dens.col = dens.col))
  }, x_name = x_name, y_name = y_name)
  
  maxx <- max(m[, x_name], na.rm = T) + 0.5
  maxy <- max(m[, y_name], na.rm = T) + 0.5
  
  (p <- ggplot(aes_string(x = x_name, y = y_name, color = "dens.col", size = 1), data = m)
    + facet_wrap(grouping)
    + geom_point()
    + scale_colour_identity()
    + scale_size_identity()
    + xlim(0, maxx)
    + ylim(0, maxy)
  )
  
  return(p)
}

export_clusters <- function(working.dir, sel.graph, sel.nodes){
  d <- gsub(".txt$", ".all_events.RData", sel.graph)
  d <- file.path(working.dir, d)
  d <- my_load(d)
  clus <- as.numeric(gsub("c", "", sel.nodes))
  d <- d[d$cellType %in% clus,]
  f <- flowFrame(as.matrix(d))
  p <- sprintf("scaffold_export_%s_", gsub(".fcs.clustered.txt", "", sel.graph))
  outname <- tempfile(pattern = p, tmpdir = working.dir, fileext = ".fcs")
  # f <- cytofCore.updateFlowFrameKeywords(f)
  write.FCS(f, outname)
}

get_graph_centering_transform <- function(x, y, svg.width, svg.height){
  padding <- 50
  G.width <- max(x) - min(x)
  G.height <- max(y) - min(y)
  scaling <- max(c(G.width / (svg.width - (padding * 2)), G.height / (svg.height - (padding * 2))))
  
  x <- x / scaling
  y <- y / scaling
  
  offset.y <- min(y) - padding
  graph.x.center <- (min(x) + max(x)) / 2
  offset.x <- graph.x.center - (svg.width / 2)
  
  return(list(offset.x = offset.x, offset.y = offset.y, scaling = scaling))
}

get_vertex_size <- function(sc.data, sel.graph, figure.width, node.size.attr, min.node.size, max.node.size, landmark.node.size){
  G <- sc.data$graphs[[sel.graph]]
  size.attr <- get.vertex.attribute(G, node.size.attr)
  ret <- size.attr / sum(size.attr, na.rm = T)
  ret <- rescale_size(max.node.size, min.node.size, sc.data$dataset.statistics$max.marker.vals[["popsize.relative"]], ret)
  ret[V(G)$type == 1] <- landmark.node.size
  return(ret)
}

rescale_size <- function(max.size, min.size, max.val, x)
{
  return(((max.size - min.size) * x) / max.val + min.size);
}

add_missing_columns <- function(m, col.names, fill.data){
  v <- col.names[!(col.names %in% colnames(m))]
  print(sprintf("Adding missing columns: %s", paste(v, collapse = ", ")))
  ret <- matrix(nrow = nrow(m), ncol = length(v), data = fill.data)
  colnames(ret) <- v
  ret <- data.frame(m, ret, check.names = F)
  return(ret)
}

######## 4 Mapping dataset #####################################################
################################################################################

run_analysis_existing <- function(scaffold, 
                                  refClusteredFiles, clusteredTables, outputDir, col.names.matrix, col.names.map,
                                  inter.cluster.connections, mode, col.names.inter_cluster,
                                  inter_cluster.weight_factor, overlap_method, ew_influence){
  
  files.list <- clusteredTables
  print(sprintf("Markers used for SCAFFoLD: %s", paste(col.names.map, collapse = ", ")))
  print(paste("Using as reference", refClusteredFiles[1], sep = " "))
  ref.scaffold.file <- refClusteredFiles[2]
  ref.scaffold.data <- scaffold
  ref.scaffold.markers <- ref.scaffold.data$scaffold.col.names
  l <- load_existing_layout(ref.scaffold.data)
  tab.attractors <- l$tab.attractors
  G.attractors <- l$G.attractors
  att.labels <- V(G.attractors)$Label
  
  
  ret <- process_files(files.list, scaffold$G[[1]], G.attractors, tab.attractors, att.labels, col.names.matrix = col.names.matrix, col.names.gated = col.names.map,
                       scaffold.mode = "existing", ref.scaffold.markers = ref.scaffold.markers)
  ret <- c(list(scaffold.col.names = col.names.map, landmarks.data = ref.scaffold.data$landmarks.data), ret)
  
  # if (mode == "Concatenation") {
  init <- length(ref.scaffold.data)
  retgraphs <- as.vector(ret$graphs)
  graphnames <- names(retgraphs)
  
  for (i in c(1:length(retgraphs))) {
    ref.scaffold.data$graphs[names(retgraphs)[i]] <- retgraphs[i]
  }
  ret <- ref.scaffold.data
  my_save(ret, paste(outputDir, sprintf("%s.scaffold", refClusteredFiles[1]), sep = ""))
  # }
  # else {
  #   my_save(ret, paste(outputDir, sprintf("%s.scaffold", refClusteredFiles[1]), sep = "/"))
  # }
  return(ret)
}

load_existing_layout <- function(scaffold.data){
  G <- scaffold.data$graphs[[1]]
  G <- induced.subgraph(G, V(G)$type == 1, impl = "copy_and_delete")
  tab <- get_vertex_table(G)
  V(G)$name <- 1:vcount(G)
  return(list(G.attractors = G, tab.attractors = tab))
}

######## 5 Export Scaffold #####################################################
################################################################################

scaffold_pop_export <- function(scaffold.data){
  map <- scaffold.data
  all.table <- lapply(c(1:length(map$graphs)), function(i){
    G <- map$graphs[[i]]
    x <- V(G)$x
    y <- V(G)$y
    popsize <- V(G)$popsize
    if(is.null(popsize)){popsize <- (rep(1,length(x)))}
    tab <- data.frame(x,y,popsize)#,label)
    return(tab)
  })
  names(all.table) <- names(map$graphs)
  return(all.table)
}

scaffold_node_export <- function(scaffold.data){
  all.table <- scaffold_pop_export(scaffold.data)
  node.index.x <- which(V(scaffold.data$graphs[[1]])$type==1)
  
  table.node <- data.frame(all.table[[1]][node.index.x,c(1,2)])
  G <- scaffold.data$graphs[[1]]
  pop.name <- unlist(lapply(c(1:dim(table.node)[1]),function(x){return(names(G[[x]]))}))
  table.node <- cbind(pop.name, c(1:nrow(table.node)),table.node)
  colnames(table.node) <- c("pop Names","pop ID","Map.x","Map.y")
  
  return(table.node)
}

scaffold_cluster_export <- function(list1,list2,list.txt,scaffold.data){
  node.index.x <- which(V(scaffold.data$graphs[[1]])$type==1)
  all.table <- scaffold_pop_export(scaffold.data)
  
  multi.tables <- lapply(c(1:length(list1)), function(x){
    # print(x)
    d <- all.table[[list1[x]]]
    txt <- list.txt[[list2[x]]]
    
    table <- data.frame(d[-node.index.x,])
    
    tot <- sum(table[,"popsize"])
    # print(tot)
    temp <- as_data_frame(scaffold.data$graphs[[list1[x]]])
    pop <- temp[which(temp[,"edge_type"]=="highest_scoring"),"from"]
    
    # print(pop)
    table <- cbind(table,table[,3],pop,txt)
    table[,3] <- (table[,3]/tot)*100
    colnames(table) <- c("Scaffold.x","Scaffold.y","PercentOfTotal","Cluster.Size","pop",colnames(txt)) #clusterID,x,y,Events,Percentile,MFIs,
    
    return(table)
  })
  return(multi.tables)
}

scaffold_events_export <- function(list1, list2, list.flow.frames, scaffold.data, marker_e){

  
  landmark <- rbind(scaffold_node_export(scaffold.data), c("null.landmark", 0, NA, NA))
  all.table <- scaffold_pop_export(scaffold.data)
  print(all.table)
  node.index.x <- which(V(scaffold.data$graphs[[1]])$type==1)
  
  flow.frames.celltype <- lapply(c(1:length(list1)),function(x){
  
    fcs <- list.flow.frames[[list2[x]]]
   
    d <- all.table[[list1[x]]]
    table <- get_matrix_from_fcs(1, list(fcs), "mean",marker_e)
    # print(table)
    tot <- sum(table[,"popsize"])
    
    temp <- igraph::as_data_frame(scaffold.data$graphs[[list1[x]]])
    pop <- temp[which(temp[,"edge_type"]=="highest_scoring"),c("from", "to")]
    m.pop <- max(as.numeric(unique(pop[,2])))-min(as.numeric(unique(pop[,2])))+1
    l.pop <- length(as.numeric(unique(pop[,2])))
    
    diff <- m.pop - l.pop
    
    if (m.pop != l.pop)
    {
      v <- c(1:m.pop+10)
      for (i in v)
      {
        if (!i%in%pop[,2] && diff>0)
        {
          temp <- pop[(i-9):nrow(pop),]
          pop <- rbind(pop[1:(i-10),], c("null.landmark", i))
          pop <- rbind(pop, temp)
          diff <- diff -1
        }
      }
    }
    
    pop <- pop[,"from"]
   
    table <- cbind(table,table[,"popsize"],pop)
    table[,dim(table)[2]-1] <- (as.numeric(table[,"popsize"])/tot)*100
    colnames(table)[dim(table)[2]-1] <- "PercentOfTotal"
    
    rdata <- fcs@exprs
    new_col.1 <- matrix(rdata[,marker_e], nrow = nrow(rdata), ncol = 1, dimnames = list(NULL, marker_e))
    
    #Used when empty clusters exist with files coming from external sources
    if(length(unique(table[,marker_e])) != max(as.numeric(unique(table[,marker_e]))))
    {
      liste.true <- unique(table[,marker_e])
      liste.max <- 1:max(as.numeric(unique(table[,marker_e])))
      matrice.miss <- !(liste.max %in% liste.true)
      id <- 1
      for (i in 1:length(matrice.miss)) {
        if(matrice.miss[i]) {
          empty.clusters <- matrix(0, ncol=ncol(table), nrow = 1)
          colnames(empty.clusters) <- colnames(table)
          empty.clusters[,1] <- i
          empty.clusters[,"sample"] <- 1
          empty.clusters[,"pop"] <- "None"
          temp <- table[1:i-1,]
          temp <- rbind(temp, empty.clusters)
          table <- rbind(temp, table[i:nrow(table),])
        }
      }
    }
    
    new_col.2 <- as.vector(table[new_col.1,"pop"])
    
    pop.index <- unlist(lapply(new_col.2,function(j){return(landmark[which(j==landmark[,"pop Names"]),"pop ID"])}))
    
    popIDscaffold <- as.matrix(as.numeric(pop.index))
    colnames(popIDscaffold) <- "popIDscaffoldBis"
    
    fcs.2 <- enrich.FCS.CIPHE(fcs, popIDscaffold)
    # fcs.2 <- flowCore::cbind2(fcs, popIDscaffold)
    return(fcs.2)
  })
  return(flow.frames.celltype)
}


scaffold_pop_mfi <- function(list1, list2, list.flow.frames, scaffold.data, marker_e, methods="median"){
  
  landmark <- scaffold_node_export(scaffold.data)
  all.table <- scaffold_pop_export(scaffold.data)
  node.index.x <- which(V(scaffold.data$graphs[[1]])$type==1)
  
  table.mfi.pop <- lapply(c(1:length(list1)),function(x){
    fcs <- list.flow.frames[[list2[x]]]
    colnames(fcs@exprs) <- pData(fcs@parameters)[,"desc"]
    colnames(fcs@exprs)[which(is.na(colnames(fcs@exprs) == "NA"))] <- pData(fcs@parameters)[,"name"][which(is.na(colnames(fcs@exprs) == "NA"))]
    d <- all.table[[list1[x]]]
    
    table <- get_matrix_from_fcs(1, list(fcs), "mean",marker_e)
    
    tot <- sum(table[,"popsize"])
    
    temp <- as_data_frame(scaffold.data$graphs[[list1[x]]])
    
    pop <- temp[which(temp[,"edge_type"]=="highest_scoring"),"from"]
    
    table <- cbind(table,table[,"popsize"],pop)
    
    table[,dim(table)[2]-1] <- (as.numeric(table[,"popsize"])/tot)*100
    
    colnames(table)[dim(table)[2]-1] <- "PercentOfTotal"
    
    rdata <- fcs@exprs
    
    new_col.1 <- matrix(rdata[,marker_e], nrow = nrow(rdata), ncol = 1, dimnames = list(NULL, marker_e))
    
    # browser()
    
    #Used when empty clusters exist with files coming from external sources
    if(length(unique(table[,marker_e])) != max(as.numeric(unique(table[,marker_e]))))
    {
      liste.true <- unique(table[,marker_e])
      liste.max <- 1:max(as.numeric(unique(table[,marker_e])))
      matrice.miss <- !(liste.max %in% liste.true)
      id <- 1
      for (i in 1:length(matrice.miss)) {
        if(matrice.miss[i]) {
          empty.clusters <- matrix(0, ncol=ncol(table), nrow = 1)
          colnames(empty.clusters) <- colnames(table)
          empty.clusters[,1] <- i
          empty.clusters[,"sample"] <- 1
          empty.clusters[,"pop"] <- "None"
          temp <- table[1:i-1,]
          temp <- rbind(temp, empty.clusters)
          table <- rbind(temp, table[i:nrow(table),])
        }
      }
    }
    
    new_col.2 <- as.vector(table[new_col.1,"pop"])
    pop.index <- unlist(lapply(new_col.2,function(j){return(landmark[which(j==landmark[,"pop Names"]),"pop ID"])}))
    
    popIDscaffold <- as.matrix(as.numeric(pop.index))
    colnames(popIDscaffold) <- "popIDscaffold"
    
    fcs <- enrich.FCS.CIPHE(fcs, popIDscaffold)
    # fcs <- flowCore::cbind2(fcs, popIDscaffold)
    
    res <- lapply(unique(fcs@exprs[,"popIDscaffold"]),function(y){
      mat <- fcs@exprs[which(fcs@exprs[,"popIDscaffold"]==y),]
      if(dim(mat)[1]< 2 || is.null(dim(mat))){return(c(0,0))}
      return(c(apply(mat,2,mean),dim(mat)[1]))
    })
    
    table.mfi <- do.call(rbind, res)
    colnames(table.mfi)[dim(table.mfi)[2]] <- "Count"
    
    
    return(table.mfi)
  })
  return(table.mfi.pop)
}

get.available.ram <- function(){
  available.ram <- 0
  if(Sys.info()[['sysname']] == "Windows")
  {
    available.ram <- system("wmic OS get FreePhysicalMemory /Value", intern=T)[3]
    available.ram <- strsplit(available.ram, "=", fixed = T)[[1]][[2]]
    available.ram <- as.numeric(strsplit(available.ram, "\r", fixed = T)[[1]][[1]])*1024
  }
  
  return(available.ram)
}

get.nmb.cores.max <- function(files.sizes, #Liste ou vecteur contenant les tailles des fichiers en bytes
                              available.cores, #Nombre de coeursdisponibles
                              x.cores = 0.5,  #Coefficient de r?duction du nombres de coeurs(ex: x.cores=0.5 ==> la moiti? des coeurs utilisables)
                              x.ram = 0.5, #Coefficient de r?duction de la RAM
                              correction.coef = 1, #Erreur relative au calcul de la m?moire par coeur.
                              separate.by.files = T) #Un fichier par coeur si T; Tous les coeurs utilisent tous les fichiers si F
  #correction.coef = 1 ==> pas de changement
  #correction.coef = 1.2 ==> pr?voir une utilisation 20% plus importante de ram par coeur
{
  available.ram <- get.available.ram()
  max.size <- 0
  Nk <- 0
  
  if(separate.by.files){
    max.size <- max(as.numeric(unlist(files.sizes)))
    for(k in 1:as.integer(available.cores*x.cores)){
      mean.size <- length(unlist(files.sizes))*(max.size)+(3.309e+07 + 1.584*max.size)*k
      if((mean.size*correction.coef) <= available.ram*x.ram){
        Nk <- Nk + 1
      }
    }
  } else {
    for(i in 1:length(files.sizes)){
      max.size <- max.size + files.sizes[[i]]
      max.size <- max.size + 3.309e+07 + 1.584*files.sizes[[i]]  #Clara ram peak: p = 3.3e7 + 1.6 * file_size
    }
    print(max.size/1024/1024)
    for(k in 1:as.integer(available.cores*x.cores)){
      if(k*max.size*correction.coef <= available.ram*x.ram){
        Nk <- Nk + 1
      }
    }
  }
  
  return(Nk)
} 