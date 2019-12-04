#' Predict method of a Generalized Additive model for the copula parameter
#' or Kendall's tau
#'
#' @param object \code{\link{gamBiCop-class}} object.
#' @param newdata (Same as in \code{\link{predict.gam}} from the
#' \code{\link[mgcv:mgcv-package]{mgcv}} package) A matrix or data frame
#' containing the values of the model covariates at which predictions are
#' required. If this is not provided then predictions corresponding to the
#' original data are returned. If \code{newdata} is provided then it should contain all
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
#' of the linear predictor is returned separately (possibly with standard
#' errors): this includes parametric model components, followed by each smooth
#' component, but excludes any offset and any intercept. When
#' \code{type = 'lpmatrix'} then a matrix is returned which yields the values of
#' the linear predictor (minus any offset) when post-multiplied by the parameter
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
#' Otherwise, if \code{type = 'lpmatrix'} (only active for
#' \code{type = 'calib'}), then a matrix is returned which will give a vector of
#' linear predictor values (minus any offset) at the supplied covariate values,
#' when applied to the model coefficient vector (similar as
#' \code{\link{predict.gam}} from the \code{\link[mgcv:mgcv-package]{mgcv}}).
#' @seealso \code{\link{gamBiCop}} and \code{\link{gamBiCopFit}}.
#' @examples
#' require(copula)
#' set.seed(0)
#'
#' ## Simulation parameters (sample size, correlation between covariates,
#' ## Clayton copula family)
#' n <- 5e2
#' rho <- 0.5
#' fam <- 1
#'
#' ## A calibration surface depending on three variables
#' eta0 <- 1
#' calib.surf <- list(
#'   calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
#'     Tm <- (Tf - Ti) / 2
#'     a <- -(b / 3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
#'     return(a + b * (t - Tm)^2)
#'   },
#'   calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1) {
#'     a <- b * (1 - 2 * Tf * pi / (f * Tf * pi +
#'       cos(2 * f * pi * (Tf - Ti))
#'       - cos(2 * f * pi * Ti)))
#'     return((a + b) / 2 + (b - a) * sin(2 * f * pi * (t - Ti)) / 2)
#'   },
#'   calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf / 8) {
#'     Tm <- (Tf - Ti) / 2
#'     a <- (b * s * sqrt(2 * pi) / Tf) * (pnorm(0, Tm, s) - pnorm(Tf, Tm, s))
#'     return(a + b * exp(-(t - Tm)^2 / (2 * s^2)))
#'   }
#' )
#'
#' ## 3-dimensional matrix X of covariates
#' covariates.distr <- mvdc(normalCopula(rho, dim = 3),
#'   c("unif"), list(list(min = 0, max = 1)),
#'   marginsIdentical = TRUE
#' )
#' X <- rMvdc(n, covariates.distr)
#' colnames(X) <- paste("x", 1:3, sep = "")
#'
#' ## U in [0,1]x[0,1] with copula parameter depending on X
#' U <- condBiCopSim(fam, function(x1, x2, x3) {
#'   eta0 + sum(mapply(function(f, x)
#'     f(x), calib.surf, c(x1, x2, x3)))
#' }, X[, 1:3], par2 = 6, return.par = TRUE)
#'
#' ## Merge U and X
#' data <- data.frame(U$data, X)
#' names(data) <- c(paste("u", 1:2, sep = ""), paste("x", 1:3, sep = ""))
#'
#' ## Model fit with penalized cubic splines (via min GCV)
#' basis <- c(3, 10, 10)
#' formula <- ~ s(x1, k = basis[1], bs = "cr") +
#'   s(x2, k = basis[2], bs = "cr") +
#'   s(x3, k = basis[3], bs = "cr")
#' system.time(fit <- gamBiCopFit(data, formula, fam))
#'
#' ## Extract the gamBiCop objects and show various methods
#' (res <- fit$res)
#' EDF(res)
#' pred <- gamBiCopPredict(fit$res, X, target = c("calib", "par", "tau"))
#' @export
gamBiCopPredict <- function(object, newdata = NULL,
                            target = "calib", alpha = 0, type = "link") {
  if (is.list(object)) {
    if (any(!is.element(names(object), c("family", "par", "par2"))) ||
      !is.numeric(unlist(object))) {
      msg <- paste(
        "Object should be a",
        " valid gamBiCop object or a list containing three",
        " items (family, par, par2)."
      )
      return(msg)
    } else if (target == "calib" && type == "terms") {
      ztrafo <- function(x) (exp(x) - 1) / (exp(x) + 1)
      famnew <- famTrans(object$family, inv = FALSE, par = object$par)
      val <- ztrafo(BiCopPar2Tau(famnew, object$par, object$par2))
      out <- matrix(0,
        nrow = nrow(as.matrix(newdata)),
        ncol = ncol(as.matrix(newdata))
      )

      colnames(out) <- colnames(as.matrix(newdata))
      out <- list(calib = out)
      attr(out$calib, "constant") <- val
      return(out)
    } else {
      stop("Prediction method not implemented for simplified copulas.")
    }
  }

  stopifnot(valid.gamBiCop(object))

  mm <- object@model

  family <- object@family
  # Define links between Kendall's tau, copula parameter and calibration
  # function... the cst/cstinv make sure that the boundaries are never attained
  if (family %in% c(1, 2)) {
    cstpar <- csttau <- function(x) x * (1 - 1e-8)
  } else {
    csttau <- function(x) x * (1 - 1e-8)
    cstpar <- function(x) x
  }

  par2tau.fun <- function(x) csttau(par2tau(cstpar(x), family))
  tau2par.fun <- function(x) cstpar(tau2par(csttau(x), family))
  eta2par.fun <- function(x) cstpar(eta2par(x, family))
  par2eta.fun <- function(x) par2eta(cstpar(x), family)
  eta2tau.fun <- function(x) csttau(eta2par(x, 1))
  tau2eta.fun <- function(x) par2eta(csttau(x), 1)

  out <- list()
  sel <- length(target) == 1 && (target == "calib")
  if (sel) {
    if (is.null(newdata)) {
      out$calib <- predict.gam(mm, type = type)
    } else {
      out$calib <- predict.gam(mm, as.data.frame(newdata), type = type)
    }
  } else {
    if (is.null(newdata)) {
      out$calib <- predict.gam(mm)
    } else {
      out$calib <- predict.gam(mm, as.data.frame(newdata))
    }
  }

  quantile.fun <- function(x) quantile(x, c((1 - alpha) / 2, 1 - (1 - alpha) / 2))
  myCI <- function(x) t(apply(x, 2, quantile.fun))
  if (!sel && (alpha != 0) && (alpha != 1)) {
    if (is.null(newdata)) {
      Xp <- predict.gam(mm, type = "lpmatrix")
    } else {
      Xp <- predict.gam(mm, as.data.frame(newdata), type = "lpmatrix")
    }
    b <- coef(mm)
    Vp <- vcov(mm)
    br <- mvrnorm(1e4, b, Vp)
    calib <- br %*% t(Xp)
    out$calib.CI <- myCI(calib)
  } else {
    calib <- NULL
  }

  if (any(is.element(target, "par"))) {
    if (object@tau) {
      out$par <- tau2par.fun(eta2tau.fun(out$calib))
      if (!is.null(calib)) {
        out$par.CI <- myCI(tau2par.fun(eta2tau.fun(calib)))
      }
    } else {
      out$par <- eta2par.fun(out$calib)
      if (!is.null(calib)) {
        out$par.CI <- myCI(eta2par.fun(calib))
      }
    }
  }

  if (any(is.element(target, "tau"))) {
    if (object@tau) {
      out$tau <- eta2tau.fun(out$calib)
      if (!is.null(calib)) {
        out$tau.CI <- myCI(eta2tau.fun(calib))
      }
    } else {
      out$tau <- par2tau.fun(eta2par.fun(out$calib))
      if (!is.null(calib)) {
        out$tau.CI <- myCI(par2tau.fun(eta2par.fun(calib)))
      }
    }
  }

  return(out)
}


valid.gamBiCopPredict <- function(object, newdata, target, alpha, type) {
  if (!valid.gamBiCop(object)) {
    return("'gamBiCopPredict' can only be used to predict from 'gamBiCop' objects")
  }
  if (!is.character(target) || !is.vector(target)) {
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
    return(paste(
      "'target' should be either 'calib', 'par', 'tau', or a",
      "vector containing two or three of those."
    ))
  }

  if (target == "calib" && (length(type) != 1 ||
    !is.element(type, c("link", "terms", "lpmatrix")))) {
    return(paste(
      "When 'target' is 'calib', then 'type' should be either",
      "'link', 'terms' or 'lpmatrix'."
    ))
  }

  if (!valid.unif(alpha)) {
    return(msg.unif(var2char(alpha)))
  }

  return(TRUE)
}
