gamVineStructureSelect <- function(data, type = 0, selectioncrit = "AIC", 
                                  indeptest = FALSE, level = 0.05, trunclevel = NA, 
                                method = "NR", tol.rel = 0.001, n.iters = 10, 
                                verbose = FALSE) {
  
  chk <- valid.gamVineStructureSelect(bla)
  if (chk != TRUE) {
    return(chk)
  } 
  
  if (type == 0) {
    type <- "RVine"
  } else {
    type <- "CVine"
  }
  
  dataset <- data.frame(dataset)
  n <- dim(dataset)[1]
  d <- dim(dataset)[2] 
  
  if (is.null(colnames(dataset))) {
    nn <- paste("V",1:d,sep="") 
    colnames(dataset) <- nn
  } else {
    nn <- colnames(dataset)
  }  
  
  if (is.na(trunclevel)) {
    trunclevel <- d
  }
  
  out <- list(tree <- vector("list", d), graph <- vector("list", d))
  
  graph <- initFirstGraph(data)
  mst <- findMST(graph, type)
  tree <- fitFirstTree(mst, data, familyset, selectioncrit, indeptest, level)
  
  out$tree[[1]] <- tree
  out$graph[[1]] <- graph
  oldtree  <- tree
  
  
  for(i in 2:(d-1)){    
    if(trunclevel == i-1) {
      familyset <- 0
    }
    
    graph <- buildNextGraph(tree)
    mst <- findMST(graph, type) 
    tree <- fitTree(mst, tree, familyset, selectioncrit, 
                    indeptest, level, progress)
    
    out$Tree[[i]] <- tree
    out$Graph[[i]] <- graph
  }
  
  return(as.GVC(out))
}

initFirstGraph <- function(data) {
  
  C <- cor(data,method="kendall")
  rownames(C) <- colnames(C) <- colnames(data)
  
  g <- graph.adjacency(C, mode="lower", weighted=TRUE, diag=FALSE)
  E(g)$tau <- E(g)$weight
  E(g)$weight <- 1-abs(E(g)$weight)
  E(g)$name <- paste(get.edgelist(g)[,1], get.edgelist(g)[,2], sep=",")  
  for(i in 1:ecount(g)){
    E(g)$conditionedSet[[i]] <- get.edges(g,i)
  }
  
  return(g)
}

findMST <- function(g, mode = "RVine") {
  
  if (mode == "RVine") {
    return(minimum.spanning.tree(g, weights=E(g)$weight))
  } else {
    M <- abs(get.adjacency(g, attr="weight", sparse=0))
    root <- which.min(rowSums(M))    
    Ecken <- get.edges(g,1:ecount(g))
    pos <- Ecken[,2] == root | Ecken[,1] == root
    return(delete.edges(g, E(g)[!pos]))
  }
}

fitFirstTree <- function(mst, data, type, copulaSelectionBy, 
                         testForIndependence, testForIndependence.level) {
  
  d <- ecount(mst)
  
  parameterForACopula <- vector("list", d)
  
  for(i in 1:d) {
    parameterForACopula[[i]] <- list()
    
    a <- get.edges(mst,i)
  
    parameterForACopula[[i]]$zr1 <- data[,a[1]]
    parameterForACopula[[i]]$zr2 <- data[,a[2]]
    
    E(mst)[i]$Data1 <-  list(data[,a[1]])
    E(mst)[i]$Data2 <-  list(data[,a[2]])
    
    if (is.null(V(mst)[a[1]]$name)) {
      E(mst)[i]$CondName1 <- a[1]
    } else {
      E(mst)[i]$CondName1 <- V(mst)[a[1]]$name
    }
         
    if(is.null(V(mst)[a[2]]$name)) {
      E(mst)[i]$CondName2 <- a[2]
    } else {
      E(mst)[i]$CondName2 <- V(mst)[a[2]]$name
    }
    
    if(is.null(V(mst)[a[1]]$name) || is.null(V(mst)[a[2]]$name)) {
      E(mst)[i]$Name <- paste(a[1], a[2], sep=" , ")
    } else {
      E(mst)[i]$Name <- paste(V(mst)[a[1]]$name, 
                                     V(mst)[a[2]]$name, sep=" , ")
    } 
  }
  
  outForACopula <- lapply(X <- parameterForACopula, FUN <- wrapper.fitACopula, 
                          type, copulaSelectionBy, testForIndependence, 
                          testForIndependence.level)
  
  for(i in 1:d)
  {
    E(mst)$param[[i]] <- c(outForACopula[[i]]$par,outForACopula[[i]]$par2)
    E(mst)[i]$type <- outForACopula[[i]]$family
    E(mst)[i]$out <- list(outForACopula[[i]])
    
    E(mst)[i]$CondData1 <- list(outForACopula[[i]]$CondOn1)
    E(mst)[i]$CondData2 <- list(outForACopula[[i]]$CondOn2)
  }
  
  return(mst)
}

fitTree <- function(mst, oldVineGraph, type,copulaSelectionBy,testForIndependence,testForIndependence.level,progress,weights=NA)
{
  d <- ecount(mst)
  
  parameterForACopula <- list()
  
  for(i in 1:d)
  {
    parameterForACopula[[i]] <- list()
    
    con <- get.edge(mst,i)
    tmp <- get.edges(oldVineGraph,con)
    
    if ((tmp[1,1] == tmp[2,1])|| (tmp[1,2] == tmp[2,1])) {
      same <- tmp[2,1]
    } else {
      if ((tmp[1,1] == tmp[2,2]) || (tmp[1,2] == tmp[2,2])) {
        same <- tmp[2,2]
      }
    }
    
    other1 <- tmp[1,tmp[1,] != same]
    other2 <- tmp[2,tmp[2,] != same]
    
    if(tmp[1,1] == same){
      zr1 <- E(oldVineGraph)[con[1]]$CondData2
      n1 <- E(oldVineGraph)[con[1]]$CondName2
    }else{
      zr1 <- E(oldVineGraph)[con[1]]$CondData1
      n1 <- E(oldVineGraph)[con[1]]$CondName1
    }
    
    if(tmp[2,1] == same){
      zr2 <- E(oldVineGraph)[con[2]]$CondData2
      n2 <- E(oldVineGraph)[con[2]]$CondName2
    }else{
      zr2 <- E(oldVineGraph)[con[2]]$CondData1
      n2 <- E(oldVineGraph)[con[2]]$CondName1
    }
    
    if(is.list(zr1)){
      zr1a=as.vector(zr1[[1]])
      zr2a=as.vector(zr2[[1]])
      n1a=as.vector(n1[[1]])
      n2a=as.vector(n2[[1]])
    }
    else{
      zr1a=zr1
      zr2a=zr2
      n1a=n1
      n2a=n2
    }
    
    if(progress == TRUE) message(n1a," + ",n2a," --> ", E(mst)[i]$name)
    
    parameterForACopula[[i]]$zr1 <- zr1a
    parameterForACopula[[i]]$zr2 <- zr2a
    
    E(mst)[i]$Data1 <-  list(zr1a)
    E(mst)[i]$Data2 <-  list(zr2a)
    
    E(mst)[i]$CondName2 <- n1a
    E(mst)[i]$CondName1 <- n2a
  }
  
  outForACopula <- lapply(X <- parameterForACopula, FUN <- wrapper.fitACopula, type,copulaSelectionBy,testForIndependence,testForIndependence.level,weights)
  
  for(i in 1:d)
  {
    E(mst)$param[[i]] <- c(outForACopula[[i]]$par,outForACopula[[i]]$par2)
    E(mst)[i]$type <- outForACopula[[i]]$family
    E(mst)[i]$out <- list(outForACopula[[i]])
    
    E(mst)[i]$CondData2 <- list(outForACopula[[i]]$CondOn1)
    E(mst)[i]$CondData1 <- list(outForACopula[[i]]$CondOn2)
  }
  
  return(mst)
}	

buildNextGraph <- function(graph) {
  
  EL <- get.edgelist(graph)
  d <- ecount(graph)
   
  g <- graph.full(d)
  V(g)$name <- E(graph)$name
  V(g)$conditionedSet <- E(graph)$conditionedSet
  
  if (!is.null(E(graph)$conditioningSet)) {
    V(g)$conditioningSet <- E(graph)$conditioningSet
  }
  
  for(i in 1:ecount(g)){
    
    con <- get.edge(g, i)
    tmp <- get.edges(graph, con)
    
    ok <- FALSE    
    if ((tmp[1,1] == tmp[2,1])|| (tmp[1,2] == tmp[2,1])) {
      ok <- TRUE
      same <- tmp[2,1]
    } else {
      if ((tmp[1,1] == tmp[2,2]) || (tmp[1,2] == tmp[2,2])) {
        ok <- TRUE
        same <- tmp[2,2] 
      }
    }
    
    if(ok){
      other1 <- tmp[1,tmp[1,] != same]
      other2 <- tmp[2,tmp[2,] != same]
      
      if (tmp[1,1] == same) {
        zr1 <- E(graph)[con[1]]$CondData2
      } else {
        zr1 <- E(graph)[con[1]]$CondData1
      }
      
      if (tmp[2,1] == same) {
        zr2 <- E(graph)[con[2]]$CondData2
      } else {
        zr2 <- E(graph)[con[2]]$CondData1
      }
      
      if (is.list(zr1)) {
        zr1a <- as.vector(zr1[[1]])
        zr2a <- as.vector(zr2[[1]])
      } else{
        zr1a <- zr1
        zr2a <- zr2
      }
      noNAs <- !(is.na(zr1a) | is.na(zr2a))
      
      #### TO MODIFY HERE
      #### TO MODIFY HERE
      #### TO MODIFY HERE
      #### TO MODIFY HERE
      #### TO MODIFY HERE
      browser()
      E(g)[i]$weight <- cor(zr1[noNAs], zr2[noNAs], method="kendall")
      
      name.node1 <- strsplit(V(g)[con[1]]$name,split=" *[,|] *")[[1]]
      name.node2 <- strsplit(V(g)[con[2]]$name,split=" *[,|] *")[[1]]   
      intersection <- c()
      
      for (j in 1:length(name.node1)) {
        for (k in 1:length(name.node2)) {
          if (name.node1[j] == name.node2[k]) {
            intersection <- c(intersection,name.node1[j])
            name.node1[j] <- ""
            name.node2[k] <- ""
            break
          }
        }
      }
      
      difference <- c()
      for (j in 1:length(name.node1)) {
        if (name.node1[j] != "") {
          difference <- c(difference, name.node1[j])
        }
      }
      
      for (j in 1:length(name.node2)) {
        if(name.node2[j] != "") {
          difference <- c(difference, name.node2[j])
        }
      }
      
      E(g)[i]$name <- paste(paste(difference, collapse= ","),
                            paste(intersection, collapse= ","), sep= " | ")
      
      if (is.list(V(g)[con[1]]$conditionedSet)) {
        l1 <- c(as.vector(V(g)[con[1]]$conditionedSet[[1]]),
                as.vector(V(g)[con[1]]$conditioningSet[[1]]))
        l2 <- c(as.vector(V(g)[con[2]]$conditionedSet[[1]]),
                as.vector(V(g)[con[2]]$conditioningSet[[1]]))
      } else {
        l1 <- c(V(g)[con[1]]$conditionedSet,V(g)[con[1]]$conditioningSet)
        l2 <- c(V(g)[con[2]]$conditionedSet,V(g)[con[2]]$conditioningSet)
      }
      out <- intersectDifference(l1,l2)
      
      suppressWarnings({E(g)$conditionedSet[i] <- list(out$difference)})
      suppressWarnings({E(g)$conditioningSet[i]  <- list(out$intersection)})
    }
    
    E(g)[i]$todel <- !ok
  }
  
  E(g)$tau <- E(g)$weight
  g <- delete.edges(g, E(g)[E(g)$todel])
  
  return(g)
}

wrapper.fitACopula <- function(parameterForACopula,type,...)
{
  return(fitACopula(parameterForACopula$zr1,parameterForACopula$zr2,type,...))
}

intersectDifference <- function(liste1, liste2) {
  out <- list()	
  out$intersection <- c()
  out$difference <- c()
  
  for(j in 1:length(liste1)){
    for(k in 1:length(liste2)){
      if(!is.na(liste2[k]) && liste1[j] == liste2[k]){
        out$intersection <- c(out$intersection, liste1[j])
        liste1[j] <- NA
        liste2[k] <- NA
        break
      }
    }
  }
  
  for(j in 1:length(liste1)){
    if(!is.na(liste1[j])){
      out$difference <- c(out$difference, liste1[j])
    }
  }
  for(j in 1:length(liste2)){
    if(!is.na(liste2[j])){
      out$difference <- c(out$difference, liste2[j])
    }
  }
  
  return(out)
}

fitACopula <- function(u1, u2, familyset=NA, selectioncrit="AIC", 
                        indeptest=FALSE,level=0.05) {
  
  out <- BiCopSelect(u1, u2, familyset, selectioncrit, indeptest, level)
  if(out$family%in%c(23,24))
  {
    out$family <- out$family+10
  }
  else if(out$family%in%c(33,34))
  {
    out$family <- out$family-10
  }
  tmp <- .C("Hfunc1", as.integer(out$family), as.integer(length(u1)), 
            as.double(u1), as.double(u2), as.double(out$par),
            as.double(out$par2), as.double(rep(0,length(u1))), 
            PACKAGE='VineCopula')[[7]]
  out$CondOn1 <- tmp
  
  tmp <- .C("Hfunc2", as.integer(out$family), as.integer(length(u1)), 
            as.double(u1), as.double(u2), as.double(out$par),
            as.double(out$par2), as.double(rep(0,length(u1))), 
            PACKAGE='VineCopula')[[7]]
  out$CondOn2 <- tmp
  
  return(out)
}

as.RVM <- function(RVine){
  
  n <- length(RVine$Tree)+1
  con <- list()
  nam <- V(RVine$Tree[[1]])$name
  
  conditionedSets <- NULL
  corresppondingParams <- list()
  corresppondingTypes <- list()
  
  if(is.list(E(RVine$Tree[[n-1]])$conditionedSet))
  {
    conditionedSets[[n-1]][[1]] <- (E(RVine$Tree[[n-1]])$conditionedSet[[1]])	
    for(k in 1:(n-2)){
      #conditionedSets[[k]] <- E(RVine$Tree[[k]])$conditionedSet[[1]]
      conditionedSets[[k]] <- E(RVine$Tree[[k]])$conditionedSet
      corresppondingParams[[k]] <- as.list(E(RVine$Tree[[k]])$param)
      corresppondingTypes[[k]] <- as.list(E(RVine$Tree[[k]])$type)
    }
    
    corresppondingParams[[n-1]] <- list()
    corresppondingParams[[n-1]] <- as.list(E(RVine$Tree[[n-1]])$param)
    corresppondingTypes[[n-1]] <- as.list(E(RVine$Tree[[n-1]])$type)
    #print(corresppondingParams)
  }
  else{
    conditionedSets[[n-1]][[1]] <- (E(RVine$Tree[[n-1]])$conditionedSet)
    for(k in 1:(n-2)){
      conditionedSets[[k]] <- E(RVine$Tree[[k]])$conditionedSet
      corresppondingParams[[k]] <- as.list(E(RVine$Tree[[k]])$param)
      corresppondingTypes[[k]] <- as.list(E(RVine$Tree[[k]])$type)
    }
    #print(conditionedSets)
    corresppondingParams[[n-1]] <- list()
    corresppondingParams[[n-1]] <- as.list(E(RVine$Tree[[n-1]])$param)
    corresppondingTypes[[n-1]] <- as.list(E(RVine$Tree[[n-1]])$type)
  }
  
  Param <- array(dim=c(n,n))
  Params2 <- array(0,dim=c(n,n))
  Type <- array(dim=c(n,n))
  M <- matrix(NA,n,n)
  
  for(k in 1:(n-1)){
    w <- conditionedSets[[n-k]][[1]][1]
    
    M[k,k] <- w
    M[(k+1),k] <- conditionedSets[[n-k]][[1]][2]
    
    Param[(k+1),k] <- corresppondingParams[[n-k]][[1]][1]
    Params2[(k+1),k] <- corresppondingParams[[n-k]][[1]][2]
    
    Type[(k+1),k] <- corresppondingTypes[[n-k]][[1]]
    
    if(k == (n-1)){
      M[(k+1),(k+1)] <- conditionedSets[[n-k]][[1]][2]
    }else{
      for(i in (k+2):n){
        for(j in 1:length(conditionedSets[[n-i+1]])){
          cs <- conditionedSets[[n-i+1]][[j]]
          cty <- corresppondingTypes[[n-i+1]][[j]]
          if(cs[1] == w){
            M[i,k] <- cs[2]
            Type[i,k] <- cty #Mathias 21.3.
            break
          } else if(cs[2] == w){
            M[i,k] <- cs[1]
            if(any(cty == c(23,24,26))){ Type[i,k] <- cty+10} #Mathias 21.3.
            if(any(cty == c(33,34,36))){ Type[i,k] <- cty-10} #Mathias 21.3.
            if(!any(cty == c(23,24,26,33,34,36))){ Type[i,k] <- cty} #Mathias 21.3.
            break
          }
        }
        Param[i,k] <- corresppondingParams[[n-i+1]][[j]][1]
        Params2[i,k] <- corresppondingParams[[n-i+1]][[j]][2]
        # changed Mathias 21.3.				Type[i,k] <- corresppondingTypes[[n-i+1]][[j]]
        
        conditionedSets[[n-i+1]][[j]] <- NULL
        corresppondingParams[[n-i+1]][[j]] <- NULL
        corresppondingTypes[[n-i+1]][[j]] <- NULL
      }
    }
    
  }
  
  M <- M#+1
  M[is.na(M)]=0
  Type[is.na(Type)]=0
  
  return(RVineMatrix(M, family <- Type, par <- Param, par2 <- Params2, names <- nam))
  
}

valid.gamVineStructureSelect <- function(bla) {
#   if (!is.element(type,c(0,1)) 
#       return("Vine model not implemented.")
#   
#   if(selectioncrit != "AIC" && selectioncrit != "BIC") 
#     return("Selection criterion not implemented.")
#   if(level < 0 & level > 1) 
#     return("Significance level has to be between 0 and 1.")
  
  return(TRUE)
}
