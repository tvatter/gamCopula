#' Sequential maximum penalized likelihood estimation of a GAM-Vine model.
#' 
#' This function estimates the parameter(s) of a Generalized Additive model 
#' (GAM) Vine model, where GAMs for individual edges are specified either for
#' the copula parameter or Kendall's tau.
#' It solves the maximum penalized likelihood estimation for the copula families 
#' supported in this package by reformulating each Newton-Raphson iteration as 
#' a generalized ridge regression, which is solved using 
#' the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param data A matrix or data frame containing the data in [0,1]^d.
#' @param GVC \code{\link{gamVine-class}} object.
#' @param method \code{'NR'} for Newton-Raphson
#' and  \code{'FS'} for Fisher-scoring (default).
#' @param tol.rel Relative tolerance for \code{'FS'}/\code{'NR'} algorithm.
#' @param n.iters Maximal number of iterations for 
#' \code{'FS'}/\code{'NR'} algorithm.
#' @param verbose \code{TRUE} if informations should be printed during the 
#' estimation and \code{FALSE} (default) for a silent version.
#' @param ... Additional parameters to be passed to \code{\link{gam}} 
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
#' # Set-up a dummy data
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
#' 
#' # The second tree contains a unique conditional copula
#' # In this second example, we take a smooth calibration function
#' 
#' # Set-up an arbitrary smooth model with 10 knots for the calibration function
#' model[[count]] <- gamBiCopEst(tmp, ~ s(x1, k = 10, fx = T), 
#'                               fam[count], n.iters = 1)$res
#' 
#' # Update the coefficients of the model
#' attr(model[[count]],"model")$coefficients <- rnorm(11)
#' plot(attr(model[[count]],"model"), se = F)
#' 
#' # Update the gamVine object
#' GVC <- gamVine(Matrix = Matrix, model = model, names = nnames)
#' summary(GVC)
#' 
#' # Simulate new data
#' N <- 1e3
#' simData <- data.frame(gamVineSim(N, GVC))
#' colnames(simData) <- nnames
#' 
#' # Fit data
#' summary(fitGVC <- gamVineSeqEst(simData, GVC))
#' plot(fitGVC)
#' @seealso \code{\link{gamVine-class}}, \code{\link{gamVineSim}} and 
#' \code{\link{gamBiCopEst}}.
gamVineSeqEst <- function(data, GVC, 
                          method = "NR", tol.rel = 0.001, n.iters = 10, 
                          verbose = FALSE) {
  
  chk <- valid.gamVineSeqEst(data, GVC, method, tol.rel, n.iters, verbose)
  if (chk != TRUE) {
    return(chk)
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
  
  oldGVC <- GVC
  oldMat <- GVC@Matrix
  o <- diag(oldMat)
  oo <- o[length(o):1]
  if (any(o != length(o):1)) {
    GVC <- gamVineNormalize(GVC)
    data <- data[,oo]
  }
  
  Mat <- GVC@Matrix
  fam <- gamVineFamily(GVC)
  MaxMat <- createMaxMat(Mat)
  CondDistr <- neededCondDistr(Mat)
  
  V <- list()
  V$direct <- array(NA, dim = c(d, d, n))
  V$indirect <- array(NA, dim = c(d, d, n))
  V$direct[d, , ] <- t(data[, d:1])
  
  model.count <- get.modelCount(d)
  
  for (i in (d - 1):1) {
    for (k in d:(i + 1)) {
      #print(model.count[k, i])
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
      
      mki <- model.count[k, i]
      if (k == d || valid.gamBiCop(GVC@model[[mki]]) != TRUE) {
        temp <- BiCopEst(zr2, zr1, fam[k, i])
        GVC@model[[mki]]$par <- temp$par
        GVC@model[[mki]]$par2 <- temp$par2
        par <- rep(temp$par, n)
        par2 <- temp$par2
      } else {
        cond <- Mat[(k +1):d, i]
        tmp <- data.frame(cbind(zr2, zr1, data[,cond]))
        names(tmp) <- c("u1","u2",nn[oo[cond]])
        GVC@model[[mki]] <- gamBiCopEst(tmp, GVC@model[[mki]]@model$formula, 
                                        fam[k, i], GVC@model[[mki]]@tau, method,
                                        tol.rel, n.iters)$res
        par <- gamBiCopPred(GVC@model[[mki]], target = "par")$par
        par2 <- GVC@model[[mki]]@par2
      }
      
      if (CondDistr$direct[k - 1, i]) {
        tmp <- rep(0, n)
        tmp <- sapply(1:n, function(x) .C("Hfunc1", as.integer(fam[k, i]), 
          as.integer(1), as.double(zr1[x]), as.double(zr2[x]), 
          as.double(par[x]), as.double(par2), as.double(tmp[x]), 
          PACKAGE = "VineCopula")[[7]])
        V$direct[k - 1, i, ] <- tmp
      }
      
      if (CondDistr$indirect[k - 1, i]) {
        tmp <- rep(0, n)
        tmp <- sapply(1:n, function(x) .C("Hfunc2", as.integer(fam[k, i]), 
          as.integer(1), as.double(zr2[x]), as.double(zr1[x]), 
          as.double(par[x]), as.double(par2), as.double(tmp[x]), 
          PACKAGE = "VineCopula")[[7]])
        V$indirect[k - 1, i, ] <- tmp
      }
    }
  }
  
  oldGVC@model <- GVC@model
  return(oldGVC)
} 

valid.gamVineSeqEst <- function(data, GVC, 
                                method, tol.rel, n.iters, 
                                verbose) {
  
  if (!is.matrix(data) && !is.data.frame(data)) {
    return("data has to be either a matrix or a data frame")
  } 
  
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  if (d < 2) {
    return("Number of dimensions has to be at least 2.")
  }
  if (n < 2) {
    return("Number of observations has to be at least 2.")
  }
  if (any(data > 1) || any(data < 0)) {
    return("Data has be in the interval [0,1].")
  }
  
  if (!valid.gamVine(GVC))
    return("gamBiVineSeqEst can only be used to estimate from gamVine objects")
  
  o <- diag(GVC@Matrix)
  if (length(o) != d)
    return("The dimension of the gamVine object is incorrect.")
  
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  
  names(data)[1:2] <- c("u1","u2")
  tmp <- valid.gamBiCopEst(data, n.iters, FALSE, tol.rel, method, verbose, 1)
  if (tmp != TRUE)
    return(tmp)
  
  return(TRUE)
}
