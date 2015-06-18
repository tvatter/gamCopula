#' Structure selection and estimation of a GAM-Vine model.
#' 
#' This function select the structure and estimates the parameter(s) of a 
#' Generalized Additive model 
#' (GAM) Vine model, where GAMs for individual edges are specified either for
#' the copula parameter or Kendall's tau.
#' It solves the maximum penalized likelihood estimation for the copula families 
#' supported in this package by reformulating each Newton-Raphson iteration as 
#' a generalized ridge regression, which is solved using 
#' the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param data A matrix or data frame containing the data in [0,1]^d.
#' @param type \code{type = 0} (default) for a R-Vine and \code{type = 1} for a
#' C-Vine 
#' @param familyset An integer vector of pair-copula families to select from 
#' (the independence copula MUST NOT be specified in this vector unless one 
#' wants to fit an independence vine!). The vector has to include at least one 
#' pair-copula family that allows for positive and one that allows for negative 
#' dependence. Not listed copula families might be included to better handle 
#' limit cases. If \code{familyset = NA} (default), selection among all 
#' possible families is performed. Coding of pair-copula families:  
#' \code{1} Gaussian, \code{2} Student t, 
#' \code{3} Clayton, \code{4} Gumbel, \code{13} Survival Clayton, 
#' \code{14} Survival Gumbel,  \code{23} Rotated (90 degrees) Clayton, 
#' \code{24} Rotated (90 degrees) Gumbel, 
#' \code{33} Rotated (270 degrees) Clayton and 
#' \code{34} Rotated (270 degrees) Gumbel.
#' @param rotations If \code{TRUE}, all rotations of the families in familyset 
#' are included.
#' @param selectioncrit Character indicating the criterion for bivariate copula 
#' selection. Possible choices: \code{selectioncrit = 'AIC'} (default) or 
#' \code{'BIC'}, as in \code{\link{BiCopSelect}} from the 
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package. 
#' @param SAtestOptions TODO;TODO;TODO;TODO;TODO;TODO;TODO!!!
#' @param indeptest Logical; whether a hypothesis test for the simplifying 
#' assumption and the independence of 
#' \code{u1} and \code{u2} is performed before bivariate copula selection 
#' (default: \code{indeptest = TRUE}; see \code{\link{BiCopIndTest}} and
#' \code{\link{SAtest}}). 
#' The independence copula is chosen for a (conditional) pair if both the 
#' simplifying assumption and the null 
#' hypothesis of independence cannot be rejected.
#' @param level Numerical; significance level of the test (default: 
#' \code{level = 0.05}).
#' @param trunclevel Integer; level of truncation.
#' @param tau \code{TRUE} (default) for a calibration fonction specified for 
#' Kendall's tau or \code{FALSE} for a calibration function specified 
#' for the Copula parameter.
#' @param method \code{'NR'} for Newton-Raphson
#' and  \code{'FS'} for Fisher-scoring (default).
#' @param tol.rel Relative tolerance for \code{'FS'}/\code{'NR'} algorithm.
#' @param n.iters Maximal number of iterations for 
#' \code{'FS'}/\code{'NR'} algorithm.
#' @param verbose \code{TRUE} if informations should be printed during the 
#' estimation and \code{FALSE} (default) for a silent version.
#' from \code{\link[mgcv:mgcv-package]{mgcv}}.
#' @return \code{gamVineSeqEst} returns a \code{\link{gamVine-class}} object.
#' @examples
#' require(VineCopula)
#' set.seed(0)
#' 
#' ## A first example with a 3-dimensional GAM-Vine
#' 
#' # Define a R-vine tree structure matrix
#' d <- 3
#' Matrix <- c(2,3,1,0,3,1,0,0,1)
#' Matrix <- matrix(Matrix,d,d)
#' nnames <- paste("x", 1:d, sep = "")
#' 
#' # Copula families for each edge
#' fam <- c(3,4,1)
#' 
#' # Parameters for the first tree (two unconditional copulas)
#' par <- c(1,2)
#' 
#' # Pre-allocate the GAM-Vine model list
#' count <- 1
#' model <- vector(mode = "list", length = d*(d-1)/2)
#' 
#' # The first tree contains only the two unconditional copulas
#' for (i in 1:(d-1)) {
#'   model[[count]] <- list(family = fam[count], par = par[count], par2 = 0)
#'   count <- count + 1
#' }
#' 
#' # The second tree contains a unique conditional copula
#' # In this first example, we take a linear calibration function (10*x-5)
#' 
#' # Set-up a dummy dataset
#' tmp <- data.frame(u1 = runif(1e2), u2 = runif(1e2), x1 = runif(1e2))
#' 
#' # Set-up an arbitrary linear model for the calibration function
#' model[[count]] <- gamBiCopEst(tmp, ~ x1, fam[count])$res
#' 
#' # Update the coefficients of the model
#' attr(model[[count]],"model")$coefficients <- c(-5, 10)
#' 
#' # Define gamVine object
#' GVC <- gamVine(Matrix = Matrix, model = model, names = nnames)
#' GVC
#' 
#' # Simulate new data
#' N <- 1e3
#' simData <- data.frame(gamVineSim(N, GVC))
#' colnames(simData) <- nnames
#' 
#' # Fit data
#' summary(fitGVC <- gamVineSeqEst(simData, GVC))
#' summary(fitGVC2 <- gamVineStructureSelect(simData, GVC))
#' 
#' @seealso \code{\link{gamVine-class}}, \code{\link{gamVineSim}}, 
#' \code{\link{gamVineSeqEst}} and \code{\link{gamBiCopEst}}.
gamVineStructureSelect <- function(data, type = 0, familyset = NA, 
                                   rotations = TRUE, selectioncrit = "AIC", 
                                   SAtestOptions = "ERC",
                                   indeptest = TRUE, level = 0.05,
                                   trunclevel = NA, tau = TRUE, method = "FS",
                                   tol.rel = 0.001, n.iters = 10, 
                                   verbose = FALSE) {
  
  chk <- valid.gamVineStructureSelect(TODO)
  if (chk != TRUE) {
    return(chk)
  } 
  
  if (type == 0) {
    type <- "RVine"
  } else {
    type <- "CVine"
  }
  
  data <- data.frame(data)
  n <- dim(data)[1]
  d <- dim(data)[2] 
  
  if (is.null(colnames(data))) {
    nn <- paste("V",1:d,sep="") 
    colnames(data) <- nn
  } else {
    nn <- colnames(data)
  }  
  
  if (is.na(trunclevel)) {
    trunclevel <- d
  }
  if (trunclevel == 0) {
    familyset <- 0
  }
  if (is.na(familyset)) {
    familyset <- 1:4
  }
  if (rotations) {
    familyset <- withRotations(familyset)
  }
  
  out <- list(tree <- vector("list", d), graph <- vector("list", d))
  
  graph <- initFirstGraph(data)
  mst <- findMST(graph, type)
  tree <- fitFirstTree(mst, data, familyset, selectioncrit, indeptest, level)
  out$Tree[[1]] <- tree
  out$Graph[[1]] <- graph
  oldtree  <- tree

  for(i in 2:(d-1)){    
    if(trunclevel == i-1) {
      familyset <- 0
    }

    graph <- buildNextGraph(tree, data, SAtestOptions)
    mst <- findMST(graph, type) 
    tree <- fitTree(mst, tree, data, familyset, selectioncrit, SAtestOptions,
                    indeptest, level, tau, method, tol.rel, n.iters, verbose)
    
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

fitFirstTree <- function(mst, data, familyset, 
                         selectioncrit, indeptest, level) {
  
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
  
  outForACopula <- lapply(parameterForACopula, wrapper.fitACopula, 
                          familyset, selectioncrit, indeptest, level)
  
  for(i in 1:d) {
    E(mst)[i]$model <- list(outForACopula[[i]]$model)
    E(mst)[i]$CondData1 <- list(outForACopula[[i]]$CondOn1)
    E(mst)[i]$CondData2 <- list(outForACopula[[i]]$CondOn2)
  }
  
  return(mst)
}

fitTree <- function(mst, oldVineGraph, data, familyset, selectioncrit, 
                    SAtestOptions, indeptest, level, 
                    tau, method, tol.rel, n.iters, verbose) {

  d <- ecount(mst)
  
  parameterForACopula <- list()
  
  for(i in 1:d) {
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
      zr1a <- as.vector(zr1[[1]])
      zr2a <- as.vector(zr2[[1]])
      n1a <- as.vector(n1[[1]])
      n2a <- as.vector(n2[[1]])
    }
    else{
      zr1a <- zr1
      zr2a <- zr2
      n1a <- n1
      n2a <- n2
    }
    
    if(verbose == TRUE) message(n1a," + ",n2a," --> ", E(mst)[i]$name)
    
    E(mst)[i]$Data1 <-  list(zr1a)
    E(mst)[i]$Data2 <-  list(zr2a)
    
    E(mst)[i]$CondName2 <- n1a
    E(mst)[i]$CondName1 <- n2a
    
    tmp <- E(mst)$conditioningSet[[i]]
    #print(names(data)[tmp])
    parameterForACopula[[i]] <- data.frame(u1 = as.numeric(zr1a), 
                                           u2 = as.numeric(zr2a), data[,tmp])
    names(parameterForACopula[[i]])[-c(1,2)] <- names(data)[tmp]
  }

  outForACopula <- lapply(parameterForACopula, wrapper.fitACopula, familyset, 
                          selectioncrit, SAtestOptions, indeptest, level, 
                          tau, method, tol.rel, n.iters)

  for (i in 1:d) {
    E(mst)[i]$model <- list(outForACopula[[i]]$model)
    E(mst)[i]$CondData2 <- list(outForACopula[[i]]$CondOn1)
    E(mst)[i]$CondData1 <- list(outForACopula[[i]]$CondOn2)
  }
  
  return(mst)
}	

buildNextGraph <- function(graph, data, SAtestOptions) {
  
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
      
      E(g)[i]$weight <- SAtest(cbind(as.numeric(zr1a[noNAs]), 
                                     as.numeric(zr2a[noNAs])), 
                               data[noNAs, out$intersection], 
                               SAtestOptions)$pValue
      #E(g)[i]$weight <- cor(zr1a[noNAs], zr2a[noNAs], method="kendall")
    }
    
    E(g)[i]$todel <- !ok
  }
  
  E(g)$tau <- E(g)$weight
  g <- delete.edges(g, E(g)[E(g)$todel])
  
  return(g)
}

wrapper.fitACopula <- function(parameterForACopula,...) {
  if (length(parameterForACopula) == 2) {
    return(fitACopula(parameterForACopula$zr1,parameterForACopula$zr2,...))
  } else {
    return(fitAGAMCopula(parameterForACopula,...))
  } 
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
  

  ## select family and estimate parameter(s) for the pair copula
  out <- BiCopSelect(u1, u2,
                     familyset,
                     selectioncrit,
                     indeptest,
                     level,
                     rotations = FALSE)
  
  if (out$family %in% c(23,24)) {
    out$family <- out$family+10
  } else if(out$family %in% c(33,34)) {
    out$family <- out$family-10
  }
  
  ## store pseudo-observations for estimation in next tree
  model <- list(family = out$family, par = out$par, par2 = out$par2)
  tmp <- BiCopHfunc(u1, u2, out$family, out$par, out$par2)
  out <- list(model = model, CondOn1 = tmp$hfunc2, CondOn2 = tmp$hfunc1)
  
  return(out)
}

fitAGAMCopula <- function(data, familyset = NA, selectioncrit = "AIC", 
                          SAtestOptions = "ERC", indeptest = FALSE, 
                          level = 0.05, tau = TRUE, 
                          method = "FS", tol.rel = 0.001, n.iters = 10) {
  out <- list()
  u1 <- data[,1]
  u2 <- data[,2]

  ## perform independence test (if asked for)
  if (indeptest == TRUE && familyset != 0) {
    p1 <- SAtest(data[,1:2], data[,-c(1,2)], SAtestOptions)$pValue
    if (p1 >= level) {
      p2 <- BiCopIndTest(u1, u2)$p.value
    } else {
      p2 <- NA
    }
    out$pValue <- c(p1, p2)
  } else {
    out$pValue <- rep(NA, 2)
  }  
  if (familyset == 0 || (!any(is.na(out$pValue)) && out$pValue[2] >= level)) {
    out$model <- list(family = 0, par = 0, par2 = 0)    
    par <- rep(0, length(u1))
    fam <- 0
    par2 <- 0
  } else {
    if (!any(is.na(out$pValue)) && out$pValue[2] < level) {
      out$model <- BiCopSelect(u1, u2, familyset, selectioncrit,
                               indeptest = FALSE, rotations = FALSE)
      par <- rep(out$model$par, length(u1))
      fam <- out$model$family
      par2 <- out$model$par2
    } else {
      #browser()
      tmp <- gamBiCopSel(data, familyset, selectioncrit, tau,
                               method, tol.rel, n.iters, parallel = FALSE)
      if (tmp$conv == 0) {
        out$model <- tmp$res
        par <- gamBiCopPred(out$model, target = "par")$par
        fam <- out$model@family
        par2 <- out$model@par2
      } else {
        #browser()
        out$model <- BiCopSelect(u1, u2, familyset, selectioncrit,
                                 indeptest = FALSE, rotations = FALSE)
        par <- rep(out$model$par, length(u1))
        fam <- out$model$family
        par2 <- out$model$par2
      }
    }
  }
  
  if (fam %in% c(23,24)) {
    fam <- fam+10
  } else if (fam %in% c(33,34)) {
    fam <- fam-10
  }

  if (isS4(out$model)) {
    attr(out$model, "family") <- fam
  } else {
    out$model$family <- fam
  }

  ## store pseudo-observations for estimation in next tree
  tmp <- t(sapply(1:length(par), function(x) 
    BiCopHfunc(u1[x], u2[x], fam, par[x], par2)))
  out$CondOn1 <- tmp[,2]
  out$CondOn2 <- tmp[,1]
  
  return(out)
}

as.GVC <- function(GVC){
  n <- length(GVC$Tree)+1
  con <- list()
  nam <- V(GVC$Tree[[1]])$name
  
  conditionedSets <- NULL
  correspondingModel <- NULL
  
  if (is.list(E(GVC$Tree[[n-1]])$conditionedSet)) {
    conditionedSets[[n-1]][[1]] <- (E(GVC$Tree[[n-1]])$conditionedSet[[1]])  
  }
  else{
    conditionedSets[[n-1]][[1]] <- (E(GVC$Tree[[n-1]])$conditionedSet)
  }
  
  for (k in 1:(n-2)) {
    conditionedSets[[k]] <- E(GVC$Tree[[k]])$conditionedSet
    correspondingModel[[k]] <- as.list(E(GVC$Tree[[k]])$model)
  }
  
  correspondingModel[[n-1]] <- E(GVC$Tree[[n-1]])$model
  
  model.count <- get.modelCount(d)  
  model <- vector("list", n*(n-1)/2)
  M <- matrix(NA,n,n)
  
  for(k in 1:(n-1)){
    w <- conditionedSets[[n-k]][[1]][1]
    
    M[k,k] <- w
    M[(k+1),k] <- conditionedSets[[n-k]][[1]][2]
    
    model[[model.count[(k+1),k]]] <- correspondingModel[[n-k]][[1]]
    if(k == (n-1)){
      M[(k+1),(k+1)] <- conditionedSets[[n-k]][[1]][2]
    } else {
      for(i in (k+2):n){
        for(j in 1:length(conditionedSets[[n-i+1]])){
          cs <- conditionedSets[[n-i+1]][[j]]
          tmp <- correspondingModel[[n-i+1]][[j]]
          if(cs[1] == w){
            M[i,k] <- cs[2]
            model[[model.count[i,k]]] <- tmp
            break
          } else if(cs[2] == w){
            M[i,k] = cs[1]
            if (isS4(tmp)) {
              fam <- attr(tmp, "family")
            } else {
              fam <- tmp$family
            }
            
            if (fam %in% c(23,24)) {
              fam <- fam+10
            } else if (fam %in% c(33,34)) {
              fam <- fam-10
            }
            
            if (isS4(tmp)) {
              attr(tmp, "family") <- fam
            } else {
              tmp$family <- fam
            }
            model[[model.count[i,k]]] <- tmp
            break
          }
        }
        
        conditionedSets[[n-i+1]][[j]] <- NULL
        correspondingModel[[n-i+1]][[j]] <- NULL
      }
    }
  }
  M[is.na(M)] <- 0
  return(gamVine(M, model, nam))
  
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
