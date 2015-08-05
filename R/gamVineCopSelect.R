#' Sequential pair-copula selection and maximum penalized likelihood estimation 
#' of a GAM-Vine model.
#' 
#' This function select the copula family and estimates the parameter(s) of a 
#' Generalized Additive model 
#' (GAM) Vine model, where GAMs for individual edges are specified either for
#' the copula parameter or Kendall's tau.
#' It solves the maximum penalized likelihood estimation for the copula families 
#' supported in this package by reformulating each Newton-Raphson iteration as 
#' a generalized ridge regression, which is solved using 
#' the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param data A matrix or data frame containing the data in [0,1]^d.
#' @param Matrix Lower triangular \code{d x d} matrix that defines the R-vine 
#' tree structure.
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
#' @return \code{gamVineCopSelect} returns a \code{\link{gamVine-class}} object.
#' @examples
#' TODO
#' @seealso \code{\link{gamVine-class}}, \code{\link{gamVineSim}} and 
#' \code{\link{gamBiCopEst}}.
gamVineCopSelect <- function(data, Matrix, familyset = NA, 
                             rotations = TRUE, familycrit = "AIC", 
                             treecrit = "Kendall", SAtestOptions = "ERC",
                             indeptest = TRUE, level = 0.05,
                             trunclevel = NA, tau = TRUE, method = "FS",
                             tol.rel = 0.001, n.iters = 10, 
                             parallel = TRUE, verbose = FALSE) {
  
  chk <- valid.gamVineCopSelect(data, Matrix,
                                familyset, rotations, familycrit, 
                                treecrit, SAtestOptions, 
                                indeptest, level, trunclevel, 
                                tau, method, tol.rel, n.iters, 
                                parallel, verbose)
  if (chk != TRUE) {
    stop(chk)
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
  if (length(familyset) == 1 && is.na(familyset)) {
    familyset <- 1:4
  }
  if (rotations) {
    familyset <- withRotations(familyset)
  }
  
  oldMat <- Matrix
  o <- diag(oldMat)
  oo <- o[length(o):1]
  Mat <- reorderRVineMatrix(Matrix)
  data <- data[,oo]
  
  MaxMat <- createMaxMat(Mat)
  CondDistr <- neededCondDistr(Mat)
  
  V <- list()
  V$direct <- array(NA, dim = c(d, d, n))
  V$indirect <- array(NA, dim = c(d, d, n))
  V$direct[d, , ] <- t(data[, d:1])
  
  model.count <- get.modelCount(d)
  model <- vector("list", d*(d-1)/2)
  for (i in (d - 1):1) {
    for (k in d:(i + 1)) {
      
      mki <- model.count[k, i]
      m <- MaxMat[k, i]
      zr1 <- V$direct[k, i, ]
      
      if (m == Mat[k, i]) {
        zr2 <- V$direct[k, (d - m + 1), ]
      } else {
        zr2 <- V$indirect[k, (d - m + 1), ]
      }
      
      if (verbose == TRUE) {
        if (k == d) 
          message(oldMat[i, i], ",", oldMat[k, i]) 
        else 
          message(oldMat[i, i], ",", oldMat[k, i], "|", 
                  paste(oldMat[(k +1):d, i], collapse = ","))
      }
      
      if (d + 1 - k > trunclevel) {
        tmp <- fitACopula(zr2, zr1, 0, familycrit, indeptest, level)
      } else {
        if (k == d) {
          tmp <- fitACopula(zr2, zr1, familyset, familycrit, indeptest, level)
        } else {
          cond <- Mat[(k +1):d, i]
          tmp <- data.frame(cbind(zr2, zr1, data[,cond]))
          names(tmp) <- c("u1","u2",nn[oo[cond]])
          tmp <- fitAGAMCopula(tmp, familyset, familycrit, 
                               treecrit, SAtestOptions, indeptest, level, 
                               tau, method, tol.rel, n.iters, parallel)
        }
      }
      model[[mki]] <- tmp$model
      if (CondDistr$direct[k - 1, i]) {
        V$direct[k - 1, i, ] <- as.numeric(tmp$CondOn2)
      }
      
      if (CondDistr$indirect[k - 1, i]) {
        V$indirect[k - 1, i, ] <- as.numeric(tmp$CondOn1)
      }
    }
  }
  return(gamVine(Matrix,model,nn))
} 

valid.gamVineCopSelect <- function(data, Matrix, 
                                   familyset, rotations, familycrit, 
                                   treecrit, SAtestOptions, 
                                   indeptest, level, trunclevel, 
                                   tau, method, tol.rel, n.iters, 
                                   parallel, verbose) {
  
  chk <- valid.gamVineStructureSelect(data, 0,
                                      familyset, rotations, familycrit, 
                                      treecrit, SAtestOptions, 
                                      indeptest, level, trunclevel, 
                                      tau, method, tol.rel, n.iters, 
                                      parallel, verbose)
  if (chk != TRUE) {
    return(chk)
  }  
  d <- dim(data)[2]

  if (!is.matrix(Matrix) || dim(Matrix)[1] != dim(Matrix)[2] || 
        max(Matrix) > d || RVineMatrixCheck(Matrix) != 1) {
    return("Matrix is not a valid R-vine matrix")
  }
    
  if (dim(Matrix)[1] != d) {
    return("The dimension of Matrix is incorrect.")
  }  
  
  return(TRUE)
}
