#' Simulation from a \code{\link{gamVine-class}} object
#'
#' @param n number of d-dimensional observations to simulate.
#' @param GVC A \code{\link[gamCopula:gamVine-class]{gamVine}} object.
#' @param U If not \code{NULL}, \code{U} is an (N,d)-matrix of U[0,1] random 
#' variates to be transformed to the copula sample.
#' @param newdata If not \code{NULL}, which is mandatory when 
#' the attribute \code{covariates} from \code{GVC} is not \code{NA}, 
#' \code{newdata} is a data frame containing the values of 
#' the model covariates at which simulations are required.
#' @return A Nxd matrix of data simulated from the given 
#' \code{\link[gamCopula:gamVine-class]{gamVine}} object.
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
#' ## The output of gamVineSimulate correspond to that of RVineSim
#' sampleRVM <- RVineSim(n, RVM,U)
#' sampleGVC <- gamVineSimulate(n, GVC,U)
#' all.equal(sampleRVM, sampleGVC)
#' 
#' ## Fit the two models and compare the estimated parameter
#' fitRVM <- RVM2GVC(RVineSeqEst(sampleRVM,RVM))
#' fitGVC <- gamVineSeqFit(sampleGVC,GVC)
#' all.equal(simplify2array(attr(fitRVM, "model")),
#' simplify2array(attr(fitGVC, "model")))
#' @export
gamVineSimulate <- function(n, GVC, U = NULL, newdata = NULL) {
  
  tmp <- valid.gamVineSimulate(n, GVC, U, newdata)
  if (tmp != TRUE) {
    stop(tmp)
  }
  
  covariates <- GVC@covariates
  if (!(length(covariates) == 1 && is.na(covariates))) {
    l <- length(covariates)
  } else {
    l <- 0
  }
  
  if (!is.null(newdata)) {
    newdata <- as.data.frame(newdata)
    n <- dim(newdata)[1]
  }
  
  if (!is.null(U)) {
    U <- as.matrix(U)
    n <- dim(U)[1]
  }
  
  o <- diag(GVC@Matrix)
  d <- length(o)
  GVC <- gamVineNormalize(GVC)
  fam <- gamVineFamily(GVC)
  MaxMat <- createMaxMat(GVC@Matrix)
  CondDistr <- neededCondDistr(GVC@Matrix)
  nn <- GVC@names
  
  model.count <- get.modelCount(d)
  
  rotate <- function(x) t(apply(x, 2, rev))
  m <- rotate(rotate(GVC@Matrix))
  fam <- rotate(rotate(fam))
  model.count <- rotate(rotate(model.count))
  maxmat <- rotate(rotate(MaxMat))
  conindirect <- rotate(rotate(CondDistr$indirect))
  
  Vdirect <- Vindirect <- array(dim = c(d, d, n))
  takeU <- !is.null(U)
  if (takeU) {
    U <- U[, rev(o)]
  } else {
    U <- matrix(runif(n * d), ncol = d)
  }
  for (i in 1:d) {
    Vdirect[i, i, ] <- U[, i]
  }
  Vindirect[1, 1, ] <- Vdirect[1, 1, ]
  
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
        vars <- vars[!is.element(vars,covariates)]
        nvars <- length(vars)
        if (nvars == 1) {
          data <- data.frame(Vdirect[1, 1, ])
          names(data) <- vars
        } else if (nvars > 1) {
          data <- data.frame(t(Vdirect[1, 1:nvars, ]))
          names(data) <- vars[rank(sapply(vars, function(x) which(nn == x)))]
        }
        if (l != 0) {
          if (nvars == 0) {
            data <- newdata
          } else {
            data <- cbind(data, newdata)
            names(data)[(nvars+1):dim(data)[2]] <- names(newdata)
          }
        }
        par <- gamBiCopPredict(model, data, target = "par")$par
        par2 <- model@par2
        #pp <- par
      } else {
        par <- rep(model$par, n)
        par2 <- model$par2
      }
      
      
      fams <- vapply(1:length(par),
                     function(j) famTrans(fam[k, i], inv = FALSE, par = par[j]),
                     numeric(1))
      if (is.null(par2))
        par2 <- 0
      
      ## bound parameters; otherwise C code may crash
      # par[par > 100] <- 100
      # par[par < -100] <- -100
      Vdirect[k, i, ] <- BiCopHinv1(zz, Vdirect[k + 1, i, ],
                                    fams, par, par2,
                                    check.pars = FALSE)
      if (i < d) {
        if (conindirect[k + 1, i] == TRUE) {
          Vindirect[k + 1, i, ] <- BiCopHfunc2(zz, Vdirect[k, i, ],
                                               fams, par, par2,
                                               check.pars = FALSE)
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

valid.gamVineSimulate <- function(n, GVC, U, newdata) {
  
  tmp <- valid.gamVine(GVC)
  if (tmp != "TRUE") {
    return(tmp)
  }
  d <- dim(GVC)
  
  covariates <- GVC@covariates
  if (!(length(covariates) == 1 && is.na(covariates))) {
    l <- length(covariates)
  } else {
    l <- 0
  }
  
  if (l != 0) {
    if (is.null(newdata)) {
      return("With this gamVine object, the argument newdata can't be null.")
    } 
    newdata <- tryCatch(as.data.frame(newdata), error = function(e) e)
    msg <- "should be or be coercisable to a data frame."
    if (any(class(newdata) != "data.frame")) {
      return(paste("newdata", msg))
    }
    if (any(!is.element(names(newdata), covariates))) {
      return("newdata should contain all the covariates required by GVC.")
    }
    n <- dim(newdata)[1]
  } 
  
  if (!is.null(U)) {
    U <- tryCatch(as.matrix(U), error = function(e) e)
    msg <- "should be or be coercisable to a matrix."
    if (any(class(U) != "matrix")) {
      return(paste("U", msg))
    }   
    if (!is.numeric(U) || dim(U)[2] != d || any(U <= 0) || any(U >= 1)) {
      msg <- paste("If not null, U should have a number of column equal to",
                   "the dimension of GVC and with elements in (0,1).")
      return(msg)
    }
    if (!is.null(newdata) && dim(U)[1] != n) {
      msg <- paste("If not null, U should have a number of rows equal to",
                   "the number of rows in newdata.")
      return(msg)
    }
    n <- dim(U)[1]
  } 
  
  if (!valid.posint(n)) {
    return(msg.posint(var2char(n)))
  }
  
  return(TRUE)
}
