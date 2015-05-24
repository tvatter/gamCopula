#' Predict method of a Generalized Additive model for the copula parameter 
#' or Kendall's tau
#'
#' @param object \code{\link{gamBiCop-class}} object.
#' @param newdata (Same as in \code{\link{predict.gam}} from the 
#' \code{\link[mgcv:mgcv-package]{mgcv}} package) A matrix or data frame 
#' containing the values of the model covariates at which predictions are 
#' required. If this is not provided then predictions corresponding to the 
#' original data are returned. If newdata is provided then it should contain all 
#' the variables needed for prediction: a warning is generated if not.
#' @param target Either \code{'calib'}, \code{'par'} or \code{'tau'} or a 
#' combination of those. \code{'calib'} (default) corresponds to the calibration 
#' function, \code{'par'} to the copula parameter and 
#' \code{'tau'} to Kendall's tau.
#' @param alpha In (0,1) to return the corresponding confidence interval.
#' @param type (Similar as in \code{\link{predict.gam}} from the 
#' \code{\link[mgcv:mgcv-package]{mgcv}} package, only active for 
#' \code{type = 'calib'}). When this has the value \code{'link'} (default), the 
#' calibration function is returned.  When \code{type = 'terms'} each component 
#' of the linear predictor is returned seperately (possibly with standard 
#' errors): this includes parametric model components, followed by each smooth 
#' component, but excludes any offset and any intercept. When 
#' \code{type = 'lpmatrix'} then a matrix is returned which yields the values of
#' the linear predictor (minus any offset) when postmultiplied by the parameter 
#' vector (in this case alpha is ignored).
#' @return If \code{target = 'calib'}, then a list with 1 item \code{calib}. 
#' If \code{target = 'par'}, \code{target = 'tau'} or 
#' \code{target = c('par', 'tau')}, 
#' then a list with 2, 2 or 3 items, namely \code{calib} and \code{par},  
#' \code{tau} and \code{par}, or  \code{calib}, \code{tau} and \code{par}.
#' 
#' If \code{alpha} is in (0,1), then a additional items of the list are 
#' \code{calib.CI} as well as e.g. \code{par.CI} and/or \code{tau.CI} depending 
#' on the value of \code{target}.
#'  
#' Otherwhise, if \code{type = 'lpmatrix'} (only active for 
#' \code{type = 'calib'}), then a matrix is returned which will give a vector of
#' linear predictor values (minus any offest) at the supplied covariate values, 
#' when applied to the model coefficient vector (similar as 
#' \code{\link{predict.gam}} from the \code{\link[mgcv:mgcv-package]{mgcv}}). 
#' @seealso \code{\link{gamBiCop}} and \code{\link{gamBiCopEst}}.
#' @examples 
#' set.seed(0)
#' 
#' ## Simulation parameters (sample size, correlation between covariates,
#' ## Clayton copula family)
#' n <- 2e2
#' rho <- 0.5
#' fam <- 3 
#' 
#' ## A calibration surface depending on three variables
#' eta0 <- 1
#' calib.surf <- list(
#'   calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
#'     Tm <- (Tf - Ti)/2
#'     a <- -(b/3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
#'     return(a + b * (t - Tm)^2)},
#'   calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1) {
#'     a <- b * (1 - 2 * Tf * pi/(f * Tf * pi +
#'                                  cos(2 * f * pi * (Tf - Ti))
#'                                - cos(2 * f * pi * Ti)))
#'     return((a + b)/2 + (b - a) * sin(2 * f * pi * (t - Ti))/2)},
#'   calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf/8) {
#'     Tm <- (Tf - Ti)/2
#'     a <- (b * s * sqrt(2 * pi)/Tf) * (pnorm(0, Tm, s) - pnorm(Tf, Tm, s))
#'     return(a + b * exp(-(t - Tm)^2/(2 * s^2)))})
#' 
#' ## 3-dimensional matrix X of covariates
#' covariates.distr <- copula::mvdc(copula::normalCopula(rho, dim = 3),
#'                                  c("unif"), list(list(min = 0, max = 1)),
#'                                  marginsIdentical = TRUE)
#' X <- copula::rMvdc(n, covariates.distr)
#' colnames(X) <- paste("x",1:3,sep="")
#' 
#' ## U in [0,1]x[0,1] with copula parameter depending on X
#' U <- CondBiCopSim(fam, function(x1,x2,x3) {eta0+sum(mapply(function(f,x)
#'   f(x), calib.surf, c(x1,x2,x3)))}, X[,1:3], par2 = 6, return.par = TRUE)
#' 
#' ## Merge U and X
#' data <- data.frame(U$data,X)
#' names(data) <- c(paste("u",1:2,sep=""),paste("x",1:3,sep=""))
#' 
#' ## Model fit with penalized cubic splines (via min GCV)
#' basis <- c(3, 10, 10)
#' formula <- ~s(x1, k = basis[1], bs = "cr") + 
#'   s(x2, k = basis[2], bs = "cr") + 
#'   s(x3, k = basis[3], bs = "cr")
#' system.time(fit <- gamBiCopEst(data, formula, fam, method = "FS"))
#' 
#' ## Extract the gamBiCop objects and show various methods
#' (res <- fit$res)
#' EDF(res)
#' pred <- gamBiCopPred(fit$res, X, target = c("calib", "par", "tau"))
#' 
#' @export
gamBiCopPred <- function(object, newdata = NULL, 
                         target = "calib", alpha = 0, type = "link") {
  
  if (!valid.gamBiCop(object)) {
    stop("gamBiCopPred can only be used to predict from gamBiCop objects")
  }
  if (!is.character(target) || !is.null(dim(target))) {
    targerr <- TRUE
  } else if (length(target) == 1 && 
               !is.element(target, c("calib", "par", "tau"))) {
    targerr <- TRUE
  } else if (length(target) > 1 && 
               !all(is.element(target, c("calib", "par", "tau")))) {
    targerr <- TRUE
  } else if (length(target) > 3) {
    targerr <- TRUE
  } else {
    targerr <- FALSE
  }
  if (targerr) {
    warning("Unknown target, reset to calib.")
    target <- "calib"
  }
  
  if (target == "calib" && type != "link" 
      && type != "terms" && type != "lpmatrix") {
    warning("Unknown type, reset to terms.")
    type <- "terms"
  }
  
  options(warn = -1)
  if (!is.na(as.double(alpha))) {
    if ((as.double(alpha) < 0) || (as.double(alpha) > 1)) {
      stop("Alpha should be a real number in [0,1]!")
    } else {
      alpha <- as.double(alpha)
    }
  } else {
    stop("Alpha should be a real number in [0,1]!")
  }
  options(warn = 0)
  
  mm <- object@model
  
  rotated <- family <- object@family
  if (is.element(rotated, c(13, 23, 33))) {
    rotated <- 3
  } else if (is.element(family, c(14, 24, 34))) {
    rotated <- 4
  }
  if (rotated %in% c(3, 4)) {
    par2tau.fun <- function(x) BiCopPar2Tau(rotated, x)
    tau2par.fun <- function(x) BiCopTau2Par(rotated, x)
    eta2par.fun <- function(x) BiCopEta2Par(rotated, x)
  }else{
    par2tau.fun <- function(x) BiCopPar2Tau(1, x)
    tau2par.fun <- function(x) BiCopTau2Par(1, x)
    eta2par.fun <- function(x) BiCopEta2Par(1, x)
  }

  out <- list()
  sel <- length(target) == 1 && (target == "calib")
  if (sel) {
    if (is.null(newdata)) {
      out$calib <- mgcv::predict.gam(mm, type = type)
    } else {
      out$calib <- mgcv::predict.gam(mm, as.data.frame(newdata), type = type)
    }
  } else {
    if (is.null(newdata)) {
      out$calib <- mgcv::predict.gam(mm)
    } else {
      out$calib <- mgcv::predict.gam(mm, as.data.frame(newdata))
    }
  }
  
  quantile.fun <- function(x) quantile(x, c((1 - alpha)/2, 1 - (1 - alpha)/2))
  if (!(type == "lpmatrix") && (alpha != 0) && (alpha != 1)) {
    Xp <- mgcv::predict.gam(mm, as.data.frame(newdata), type = "lpmatrix")
    b <- coef(mm)
    Vp <- vcov(mm)
    br <- MASS::mvrnorm(10000, b, Vp)
    calib <- br %*% t(Xp)
    out$calib.CI <- t(apply(calib, 2, quantile.fun))
  } else {
    calib <- NULL
  }
  
  if (any(is.element(target, "par"))) {
    if (object@tau) {
      tmp <- tanh(out$calib/2)
      if (rotated %in% c(3, 4)) {
        out$par <- sapply((1 + tmp)/2, tau2par.fun)
        if (family %in% c(23, 24, 33, 34)) {
          out$par <- -out$par
        }
        if (!is.null(calib)) {
          tmp <- sapply((1 + tanh(calib/2))/2, tau2par.fun)
          out$par.CI <- t(apply(tmp, 2, quantile.fun))
          if (family %in% c(23, 24, 33, 34)) {
          out$par.CI <- -out$par.CI
          }
        }
      } else {
        out$par <- sapply(tmp, tau2par.fun)
        if (!is.null(calib)) {
          tmp <- sapply(tanh(calib/2), tau2par.fun)
          out$par.CI <- t(apply(tmp, 2, quantile.fun))
        }
      }
    } else {
      out$par <- sapply(out$calib, eta2par.fun)
      if (family %in% c(23, 24, 33, 34)) {
        out$par <- -out$par
      }
      if (!is.null(calib)) {
        tmp <- sapply(calib, eta2par.fun)
        out$par.CI <- t(apply(tmp, 2, quantile.fun))
        if (family %in% c(23, 24, 33, 34)) {
          out$par.CI <- -out$par.CI
        }
      }
    }
  }
  
  if (any(is.element(target, "tau"))) {
    if (object@tau) {
      out$tau <- tanh(out$calib/2)
      if (rotated %in% c(3, 4)) {
        out$tau <- (1 + out$tau)/2
        if (family %in% c(23, 24, 33, 34)) {
          out$tau <- -out$tau
        }
      }
      if (!is.null(calib)) {
        out$tau.CI <- t(apply(tanh(calib/2), 2, quantile.fun))
        if (rotated %in% c(3, 4)) {
          out$tau.CI <- (1 + out$tau.CI)/2
          if (family %in% c(23, 24, 33, 34)) {
          out$tau.CI <- -out$tau.CI
          }
        }
      }
    } else {
      tmp <- sapply(out$calib, eta2par.fun)
      out$tau <- sapply(tmp, par2tau.fun)
      if (family %in% c(23, 24, 33, 34)) {
        out$tau <- -out$tau
      }
      if (!is.null(calib)) {
        tmp <- sapply(calib, eta2par.fun)
        tmp <- sapply(tmp, par2tau.fun)
        out$tau.CI <- t(apply(tmp, 2, quantile.fun))
        if (family %in% c(23, 24, 33, 34)) {
          out$tau.CI <- -out$tau.CI
        }
      }
    }
  }
  
  return(out)
} 
