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
#' @param covariates Vector of names for the covariates.
#' @param simplified If \code{TRUE}, then a simplified PCC is fitted (which is
#' possible only if there are exogenous covariates). If \code{FALSE} (default),
#' then a non-simplified PCC is fitted.
#' @param type \code{type = 0} (default) for a R-Vine and \code{type = 1} for a
#' C-Vine.
#' @param familyset An integer vector of pair-copula families to select from 
#' (the independence copula MUST NOT be specified in this vector unless one 
#' wants to fit an independence vine!). 
#' Not listed copula families might be included to better handle 
#' limit cases. If \code{familyset = NA} (default), selection among all 
#' possible families is performed. Coding of pair-copula families:  
#' \code{1} Gaussian, 
#' \code{2} Student t, 
#' \code{5} Frank, 
#' \code{301} Double Clayton type I (standard and rotated 90 degrees), 
#' \code{302} Double Clayton type II (standard and rotated 270 degrees), 
#' \code{303} Double Clayton type III (survival and rotated 90 degrees), 
#' \code{304} Double Clayton type IV (survival and rotated 270 degrees), 
#' \code{401} Double Gumbel type I (standard and rotated 90 degrees), 
#' \code{402} Double Gumbel type II (standard and rotated 270 degrees), 
#' \code{403} Double Gumbel type III (survival and rotated 90 degrees), 
#' \code{404} Double Gumbel type IV (survival and rotated 270 degrees).
#' @param rotations If \code{TRUE}, all rotations of the families in familyset 
#' are included.
#' @param familycrit Character indicating the criterion for bivariate copula 
#' selection. Possible choices: \code{familycrit = 'AIC'} (default) or 
#' \code{'BIC'}, as in \code{\link{BiCopSelect}} from the 
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package. 
#' @param treecrit Character indicating how pairs are selected in each tree.
#' \code{treecrit = "Kendall"} uses the maxmium spanning tree of the Kendall's tau 
#' (i.e., the tree of maximal overall dependence), and 
#' \code{treecrit = "SAtest"} builds the minimum spanning tree of p-values of a
#' test of the simplifying assumption (i.e., the tree of maximal variability 
#' in conditional dependence).
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
#' @param parallel \code{TRUE} (default) for parallel selection of copula 
#' family at each edge or \code{FALSE} for the sequential version.
#' for the Copula parameter.
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
#' fam <- c(301,401,1)
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
#' # Fit data using sequential estimation assuming true model known
#' summary(fitGVC <- gamVineSeqEst(simData, GVC))
#' 
#' # Fit data using structure selection and sequential estimation
#' summary(fitGVC2 <- gamVineStructureSelect(simData))
#' 
#' @seealso  \code{\link{gamVineSeqEst}},\code{\link{gamVineCopSelect}}, 
#'  \code{\link{gamVine-class}}, \code{\link{gamVineSim}} and 
#'  \code{\link{gamBiCopEst}}.
gamVineStructureSelect <- function(data, covariates = NA, simplified = FALSE,
                                   type = 0, familyset = NA, 
                                   rotations = TRUE, familycrit = "AIC", 
                                   treecrit = "Kendall", SAtestOptions = "ERC",
                                   indeptest = TRUE, level = 0.05,
                                   trunclevel = NA, tau = TRUE, method = "FS",
                                   tol.rel = 0.001, n.iters = 10, 
                                   parallel = FALSE, verbose = FALSE) {
  
  tmp <- valid.gamVineStructureSelect(data, covariates, simplified, type, 
                                      familyset, rotations, familycrit, 
                                      treecrit, SAtestOptions, 
                                      indeptest, level, trunclevel, 
                                      tau, method, tol.rel, n.iters, 
                                      parallel, verbose)
  if (tmp != TRUE) {
   stop(tmp)
  }     
  
  ## Transform to dataframe, get dimensions, etc (see in utilsPrivate)
  tmp <- prepare.data(data, covariates, trunclevel, familyset, rotations)
  n <- tmp$n
  d <- tmp$d
  l <- tmp$l
  nn <- tmp$nn
  data <- tmp$data
  covariates <- tmp$covariates
  trunclevel <- tmp$trunclevel
  familyset <- tmp$familyset

  if (type == 0) {
    type <- "RVine"
  } else {
    type <- "CVine"
  } 
  
  out <- list(tree <- vector("list", d), graph <- vector("list", d))
  #browser()
  graph <- initFirstGraph(data[,which(!is.element(names(data),covariates))])
  mst <- findMST(graph, type)
  tree <- fitFirstTree(mst, data,  l, covariates, familyset, familycrit, 
                       treecrit, SAtestOptions, indeptest, level, 
                       tau, method, tol.rel, n.iters, parallel, verbose)
  out$Tree[[1]] <- tree
  out$Graph[[1]] <- graph
  oldtree  <- tree

  for(i in 2:(d-1)){    
    if(trunclevel == i-1) {
      familyset <- 0
    }
    graph <- buildNextGraph(tree, data, l, covariates, 
                            simplified, treecrit, SAtestOptions)
    mst <- findMST(graph, type) 
    tree <- fitTree(mst, tree, data, l, covariates, simplified,
                    familyset, familycrit, 
                    treecrit, SAtestOptions, indeptest, level, 
                    tau, method, tol.rel, n.iters, parallel, verbose)
    out$Tree[[i]] <- tree
    out$Graph[[i]] <- graph
  }
  
  return(as.GVC(out,covariates))
}

initFirstGraph <- function(data) {
  
  C <- TauMatrix(data)#cor(data,method="kendall")
  rownames(C) <- colnames(C) <- colnames(data)
  
  g <- graph.adjacency(C, mode="lower", weighted=TRUE, diag=FALSE)
  E(g)$tau <- E(g)$weight
  E(g)$weight <- 1-abs(E(g)$weight)
  E(g)$name <- paste(get.edgelist(g)[,1], get.edgelist(g)[,2], sep=",")  
  for(i in 1:ecount(g)){
    E(g)$conditionedSet[[i]] <- ends(g,i,FALSE)
  }
  
  return(g)
}

findMST <- function(g, mode = "RVine") {
  
  if (mode == "RVine") {
    return(minimum.spanning.tree(g, weights=E(g)$weight))
  } else {
    M <- abs(get.adjacency(g, attr="weight", sparse=0))
    root <- which.min(rowSums(M))    
    Ecken <- ends(g,1:ecount(g),FALSE)
    pos <- Ecken[,2] == root | Ecken[,1] == root
    return(delete.edges(g, E(g)[!pos]))
  }
}

fitFirstTree <- function(mst, data, l, covariates, familyset, familycrit, 
                         treecrit, SAtestOptions, indeptest, level, 
                         tau, method, tol.rel, n.iters, parallel, verbose) {
  
  d <- ecount(mst)
  
  parameterForACopula <- vector("list", d)
  
  for(i in 1:d) {
 
    a <- ends(mst,i,FALSE)
    #     E(mst)[i]$Data1 <-  list(data[,a[1]])
    #     E(mst)[i]$Data2 <-  list(data[,a[2]])
    
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
    
    if (l == 0) {
      parameterForACopula[[i]] <- list()
      parameterForACopula[[i]]$zr1 <- data[,a[1]]
      parameterForACopula[[i]]$zr2 <- data[,a[2]]
    } else {
      parameterForACopula[[i]] <- data.frame(as.numeric(data[,a[1]]), 
                                             as.numeric(data[,a[2]]),
                                             data[,covariates])
      names(parameterForACopula[[i]]) <- c("u1","u2",covariates)
    }
  }

  if (l == 0) {
    outForACopula <- lapply(parameterForACopula, wrapper.fitACopula, 
                            familyset, familycrit, indeptest, level)
  } else {
    outForACopula <- lapply(parameterForACopula, wrapper.fitACopula, 
                            familyset, familycrit, 
                            treecrit, SAtestOptions, indeptest, level, 
                            tau, method, tol.rel, n.iters, parallel)
  }
  
  for(i in 1:d) {
    E(mst)[i]$model <- list(outForACopula[[i]]$model)
    E(mst)[i]$CondData1 <- list(outForACopula[[i]]$CondOn1)
    E(mst)[i]$CondData2 <- list(outForACopula[[i]]$CondOn2)
  }
  
  return(mst)
}

fitTree <- function(mst, oldVineGraph, data, l, covariates, simplified,
                    familyset, familycrit, 
                    treecrit, SAtestOptions, indeptest, level, 
                    tau, method, tol.rel, n.iters, parallel, verbose) {

  d <- ecount(mst)
  
  parameterForACopula <- list()
  
  for(i in 1:d) {
    con <- ends(mst,i,FALSE)
    tmp <- ends(oldVineGraph,con,FALSE)
    
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
    
    #     E(mst)[i]$Data1 <-  list(zr1a)
    #     E(mst)[i]$Data2 <-  list(zr2a)
    
    E(mst)[i]$CondName2 <- n1a
    E(mst)[i]$CondName1 <- n2a
 
    tmp <- E(mst)$conditioningSet[[i]]
    parameterForACopula[[i]] <- data.frame(u1 = as.numeric(zr1a), 
                                           u2 = as.numeric(zr2a))
    if (simplified) {
      parameterForACopula[[i]] <- cbind(parameterForACopula[[i]], 
                                        data[,covariates])
      names(parameterForACopula[[i]])[-c(1,2)] <- covariates
    } else {
      parameterForACopula[[i]] <- cbind(parameterForACopula[[i]], data[,tmp])
      names(parameterForACopula[[i]])[-c(1,2)] <- names(data)[tmp]
      if (l != 0) {
        parameterForACopula[[i]] <- cbind(parameterForACopula[[i]], 
                                          data[,covariates])
        names(parameterForACopula[[i]])[-c(1,2,3:(2+length(tmp)))] <- covariates
      }
    }
  }
  
#   if (d < 5) {
#     browser()
#   }
  outForACopula <- lapply(parameterForACopula, wrapper.fitACopula, 
                          familyset, familycrit, 
                          treecrit, SAtestOptions, indeptest, level, 
                          tau, method, tol.rel, n.iters, parallel)

  for (i in 1:d) {
    E(mst)[i]$model <- list(outForACopula[[i]]$model)
    E(mst)[i]$CondData2 <- list(outForACopula[[i]]$CondOn1)
    E(mst)[i]$CondData1 <- list(outForACopula[[i]]$CondOn2)
  }
  
  return(mst)
}	

buildNextGraph <- function(graph, data,  l, covariates, simplified, 
                           treecrit, SAtestOptions) {
  
  EL <- get.edgelist(graph)
  d <- ecount(graph)
   
  g <- graph.full(d)
  V(g)$name <- E(graph)$name
  V(g)$conditionedSet <- E(graph)$conditionedSet
  
  if (!is.null(E(graph)$conditioningSet)) {
    V(g)$conditioningSet <- E(graph)$conditioningSet
  }
  
  for(i in 1:ecount(g)){
    
    con <- ends(g, i, FALSE)
    tmp <- ends(graph, con, FALSE)
    
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
      
      if (treecrit == "Kendall") {
        E(g)[i]$weight <- 1-abs(fasttau(zr1a[noNAs], zr2a[noNAs]))
      } else if (treecrit == "SAtest") {
        if (simplified) {
          tmp <- data[noNAs,covariates]
        } else {
          tmp <- data[noNAs, out$intersection]
          if (l != 0) {
            tmp <- cbind(tmp, data[noNAs,covariates])
          }
        }
        E(g)[i]$weight <- SAtest(cbind(as.numeric(zr1a[noNAs]), 
                                       as.numeric(zr2a[noNAs])), 
                                 tmp, SAtestOptions)$pValue
      }
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

fitACopula <- function(u1, u2, familyset=NA, familycrit="AIC", 
                        indeptest=FALSE,level=0.05,rotation=TRUE) {
  
  ## transform the familyset to codes for the VineCopula package
  fams <- famTrans(familyset, inv = FALSE, set = TRUE)
  
  ## select family and estimate parameter(s) for the pair copula
  out <- BiCopSelect(u1, u2,
                     fams,
                     familycrit,
                     indeptest,
                     level,
                     rotations = FALSE)
  fam <- famTrans(out$family, inv = TRUE, 
                  par = cor(u1,u2), familyset = familyset)
  if (rotation == TRUE) {
    if (fam %in% c(301,303,401,403)) {
      fam <- fam+1
    } else if (fam %in% c(302,304,402,404)) {
      fam <- fam-1
    } 
  }
  out$family <- fam
  
  ## store pseudo-observations for estimation in next tree
  tmp <- bicoppd1d2(cbind(u1,u2,out$par,out$par2), out$family, p=FALSE, h=TRUE)
  
  ## save the model and pseudo-observations
  model <- list(family = out$family, par = out$par, par2 = out$par2)
  out <- list(model = model, CondOn1 = tmp[2,], CondOn2 = tmp[1,])
  
  return(out)
}

fitAGAMCopula <- function(data, familyset, familycrit, 
                          treecrit, SAtestOptions, indeptest, level, 
                          tau, method, tol.rel, n.iters, parallel, 
                          rotation = TRUE) {
  out <- list()
  u1 <- data[,1]
  u2 <- data[,2]
  
  ## perform independence test (if asked for)
  if (indeptest == TRUE && familyset != 0) {
    if (treecrit == "SAtest") {
      p1 <- SAtest(data[,1:2], data[,-c(1,2)], SAtestOptions)$pValue
    } else {
      p1 <- 0
    }
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
    ## independence copula
    
    out$model <- list(family = 0, par = 0, par2 = 0)    
    par <- rep(0, length(u1))
    fam <- 0
    par2 <- 0
  } else {
    ## transform the familyset to codes for the VineCopula package
    fams <- famTrans(familyset, inv = FALSE, set = TRUE)
    
    if (!any(is.na(out$pValue)) && out$pValue[2] < level) {
      ## unconditional copula
      
      ## select family and estimate parameter(s) for the pair copula
      tmp <- BiCopSelect(u1, u2, fams, familycrit,
                         indeptest = FALSE, rotations = FALSE)
      par <- rep(tmp$par, length(u1))
      out$model$par <- tmp$par
      out$model$family <- fam <- tmp$family
      out$model$par2 <- par2 <- tmp$par2
    } else {
      ## conditional copula
      tmp <- gamBiCopSel(data, familyset, FALSE, familycrit, level, 1.5, tau,
                               method, tol.rel, n.iters, parallel)
      if (!is.character(tmp)) {
        nvar <- unique(all.vars(tmp$res@model$pred.formula))
      }
      
      if (!is.character(tmp) && length(nvar) > 0) {
        out$model <- tmp$res
        par <- gamBiCopPred(out$model, target = "par")$par
        fam <- out$model@family
        par2 <- out$model@par2
      } else {
        tmp <- BiCopSelect(u1, u2, fams, familycrit,
                                 indeptest = FALSE, rotations = FALSE)
        par <- rep(tmp$par, length(u1))
        out$model$par <- tmp$par
        out$model$family <- fam <- tmp$family
        out$model$par2 <- par2 <- tmp$par2
      }
    }
    fam <- famTrans(fam, inv = TRUE, 
                    par = cor(u1,u2), familyset = familyset)
    if (rotation == TRUE) {
      if (fam %in% c(301,303,401,403)) {
        fam <- fam+1
      } else if (fam %in% c(302,304,402,404)) {
        fam <- fam-1
      } 
    }
    if (isS4(out$model)) {
      attr(out$model, "family") <- fam
    } else {
      out$model$family <- fam
    }
  }

  ## store pseudo-observations for estimation in next tree
  tmp <- bicoppd1d2(cbind(u1,u2,par,par2), family = fam, p = FALSE, h = TRUE)
  out$CondOn1 <- tmp[2,]
  out$CondOn2 <- tmp[1,]
  
  return(out)
}

as.GVC <- function(GVC, covariates){
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
  
  model.count <- get.modelCount(n)  
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
            if (fam %in% c(301,303,401,403)) {
              fam <- fam+1
            } else if (fam %in% c(302,304,402,404)) {
              fam <- fam-1
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
  return(gamVine(M, model, nam, covariates))
  
}

valid.gamVineStructureSelect <- function(data, covariates, simplified, type,
                                          familyset, rotations, familycrit, 
                                          treecrit, SAtestOptions, 
                                          indeptest, level, trunclevel, 
                                          tau, method, tol.rel, n.iters, 
                                          parallel, verbose) {  
  
  if (!is.matrix(data) && !is.data.frame(data)) {
    return("data has to be either a matrix or a data frame")
  } 
  
  covariates <- tryCatch(as.character(covariates), error = function(e) e)
  if (!is.vector(covariates) || any(class(covariates) != "character")) {
    return("covariates should be or be coercisable to a character vector.")
  }
  if (!(length(covariates) == 1 && is.na(covariates))) {
    l <- length(covariates)
  } else {
    l <- 0
  }
  n <- dim(data)[1]
  d <- dim(data)[2] - l
  
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
    
  if (d < 2) {
    return("Number of dimensions has to be at least 2.")
  }
  if (n < 2) {
    return("Number of observations has to be at least 2.")
  }
  if (any(data[,1:d] > 1) || any(data[,1:d] < 0)) {
    return("Data has be in the interval [0,1].")
  }
  
  if (!valid.logical(simplified)) {
    return(msg.logical(var2char(simplified)))
  }
  
  if (l == 0 && simplified == TRUE) {
    return("When there are no covariates, the PCC can't be simplified.")
  }
   
  if (is.null(type) || length(type) != 1 || !is.element(type,c(0,1))) {
    return("Vine model not implemented.")
  }

  if (!valid.familyset(familyset)) {
    return(return(msg.familyset(var2char(familyset))))
  }
  
  if (!valid.logical(rotations)) {
    return(msg.logical(var2char(rotations)))
  }
  
  if(length(familycrit) != 1 || (familycrit != "AIC" && familycrit != "BIC")) {
    return("Selection criterion for copula family not implemented.")
  } 
  
  if(length(treecrit) != 1 || (treecrit != "Kendall" && treecrit != "SAtest")) {
    return("Selection criterion for the pair selection not implemented.")
  } 
    
  if(!valid.unif(level)) {
    return(msg.unif(var2char(level)))
  }   
  
  if (!is.na(trunclevel) && !valid.posint(trunclevel)) {
    return("'trunclevel' should be a positive integer or NA.")
  }
  
  names(data)[1:2] <- c("u1","u2")
  tmp <- valid.gamBiCopEst(data, n.iters, FALSE, tol.rel, method, verbose, 1)
  if (tmp != TRUE) {
    return(tmp)
  }
      
  if (!valid.logical(parallel)) {
    return(msg.logical(var2char(parallel)))
  }
  
  return(TRUE)
}