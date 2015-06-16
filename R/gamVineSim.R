#' Simulation from a \code{\link{gamVine-class}} object
#'
#' @param N number of d-dimensional observations to simulate.
#' @param GVC \code{\link{gamVine-class}} object.
#' @param U (similar as \code{\link{RVineSim}} from the from the 
#' \code{\link{VineCopula}} package) If not NULL, an (N,d)-matrix of U[0,1] 
#' random variates to be transformed to the copula sample.
#' @return A Nxd matrix of data simulated from the given 
#' \code{\link{gamVine-class}} object.
#' @examples
#' require(VineCopula)
#' 
#' ## Example adapted from RVineSim
#' 
#' ## Define 5-dimensional R-vine tree structure matrix
#' Matrix <- c(5, 2, 3, 1, 4,
#'             0, 2, 3, 4, 1,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 1)
#' Matrix <- matrix(Matrix, 5, 5)
#' 
#' ## Define R-vine pair-copula family matrix
#' family <- c(0, 1, 3, 4, 4,
#'             0, 0, 3, 4, 1,
#'             0, 0, 0, 4, 1,
#'             0, 0, 0, 0, 3,
#'             0, 0, 0, 0, 0)
#' family <- matrix(family, 5, 5)
#' 
#' ## Define R-vine pair-copula parameter matrix
#' par <- c(0, 0.2, 0.9, 1.5, 3.9,
#'          0, 0, 1.1, 1.6, 0.9,
#'          0, 0, 0, 1.9, 0.5,
#'          0, 0, 0, 0, 4.8,
#'          0, 0, 0, 0, 0)
#' par <- matrix(par, 5, 5)
#' 
#' ## Define second R-vine pair-copula parameter matrix
#' par2 <- matrix(0, 5, 5)
#' 
#' ## Define RVineMatrix object
#' RVM <- RVineMatrix(Matrix = Matrix, family = family,
#'                    par = par, par2 = par2,
#'                    names = c("V1", "V2", "V3", "V4", "V5"))
#' 
#' ## Convert to gamVine object
#' GVC <- RVM2GVC(RVM)
#' 
#' ## U[0,1] random variates to be transformed to the copula sample
#' n <- 1e2
#' d <- 5
#' U <- matrix(runif(n*d), nrow = n)
#' 
#' ## The output of gamVineSim correspond to that of RVineSim
#' sampleRVM <- RVineSim(n,RVM,U)
#' sampleGVC <- gamVineSim(n,GVC,U)
#' all.equal(sampleRVM, sampleGVC)
#' 
#' ## Fit the two models and compare the estimated parameter
#' fitRVM <- RVM2GVC(RVineSeqEst(sampleRVM,RVM)$RVM)
#' fitGVC <- gamVineSeqEst(sampleGVC,GVC)
#' all.equal(simplify2array(attr(fitRVM, "model")),
#' simplify2array(attr(fitGVC, "model")))
#' @export
gamVineSim <- function(N, GVC, U = NULL) {
  stopifnot(N >= 1)
  if (any(!isS4(GVC), !is(GVC, "gamVine"))) {
    stop("'GVC' has to be an gamVine object.")
  }
  
  o <- diag(GVC@Matrix)
  d <- length(o)
  GVC <- gamVineNormalize(GVC)
  fam <- gamVineFamily(GVC)
  MaxMat <- createMaxMat(GVC@Matrix)
  CondDistr <- neededCondDistr(GVC@Matrix)
  
  model.count <- rep(0, d^2)
  temp <- 1:(d * (d - 1)/2)
  t1 <- 1
  sel <- seq(d, d^2 - d, by = d)
  for (i in 1:(d - 1)) {
    t2 <- t1 + d - i - 1
    model.count[sel[1:(d - i)] - i + 1] <- temp[t1:t2]
    t1 <- t2 + 1
  }
  model.count <- matrix(model.count, d, d)
  
  rotate <- function(x) t(apply(x, 2, rev))
  m <- rotate(rotate(GVC@Matrix))
  fam <- rotate(rotate(fam))
  model.count <- rotate(rotate(model.count))
  maxmat <- rotate(rotate(MaxMat))
  conindirect <- rotate(rotate(CondDistr$indirect))
  
  Vdirect <- Vindirect <- array(dim = c(d, d, N))
  takeU <- !is.null(U)
  if (takeU) {
    if (!is.matrix(U)) 
      U <- rbind(U, deparse.level = 0L)
    if ((d <- ncol(U)) < 2) 
      stop("U should be at least bivariate")  # should be an (N, n) matrix
    U <- U[, rev(o)]
  } else {
    U <- matrix(runif(N * d), ncol = d)
  }
  for (i in 1:d) {
    Vdirect[i, i, ] <- U[, i]
  }
  Vindirect[1, 1, ] <- Vdirect[1, 1, ]
  count <- 1
  for (i in 2:d) {
    for (k in (i - 1):1) {
      #print(model.count[k,i])
      model <- GVC@model[[model.count[k, i]]]
      mm <- maxmat[k, i]
      if (mm == m[k, i]) {
        zz <- Vdirect[k, mm, ]
      } else {
        zz <- Vindirect[k, mm, ]
      }
      
      if (isS4(model) && is(model, "gamBiCop")) {
        vars <- all.vars(model@model$pred.formula)
        nvars <- length(vars)
        if (nvars == 1) {
          newdata <- data.frame(Vdirect[1, 1:nvars, ])
          names(newdata) <- vars
        } else {
          newdata <- data.frame(t(Vdirect[1, 1:nvars, ]))
          names(newdata) <- vars
        }
        par <- gamBiCopPred(model, newdata, target = "par")$par
        par2 <- model@par2
        #pp <- par
      } else {
        par <- rep(model$par, N)
        par2 <- model$par2
      }
      
      tmp <- rep(0, N)
      tmp <- sapply(1:N, function(x) 
        .C("Hinv1", as.integer(fam[k, i]), as.integer(1), 
           as.double(Vdirect[k + 1, i, x]), as.double(zz[x]), as.double(par[x]), 
           as.double(par2), as.double(tmp[x]), PACKAGE = "VineCopula")[[7]])
      Vdirect[k, i, ] <- tmp
      
      if (i < d) {
        if (conindirect[k + 1, i] == TRUE) {
          tmp <- sapply(1:N, function(x) 
            VineCopula::BiCopHfunc(Vdirect[k,i, x], zz[x], 
                                   fam[k, i], par[x], par2)$hfunc1)
          Vindirect[k + 1, i, ] <- tmp
        }
      }
    }
  }
  
  out <- t(Vdirect[1, , ])
  if (!is.null(GVC@names)) {
    colnames(out) <- GVC@names
  }
  out <- out[, sort(o[length(o):1], index.return = TRUE)$ix]
  
  return(out)
} 
