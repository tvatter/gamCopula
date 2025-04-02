#' Simulate from \code{\link{gamBiCop-class}} object
#'
#' @param object \code{\link{gamBiCop-class}} object.
#' @param newdata (same as in \code{\link[mgcv]{predict.gam}} from the \code{\link[mgcv:mgcv-package]{mgcv}} package) A matrix or data frame containing the values of the model covariates at which simulations are required.
#' If this is not provided then simulations corresponding to the original data are returned.
#' @param N sample size.
#' @param return.calib should the calibration function (\code{TRUE}) be returned or not (\code{FALSE})?
#' @param return.par should the copula parameter (\code{TRUE}) be returned or not (\code{FALSE})?
#' @param return.tau should the Kendall's tau (\code{TRUE}) be returned or not (\code{FALSE})?
#' @return A list with 1 item \code{data}.  When \code{N} is smaller or larger than the \code{newdata}'s number of rows
#' (or the number of rows in the original data if \code{newdata} is not provided),
#' then \code{N} observations are sampled uniformly (with replacement) among the row of \code{newdata}
#' (or the rows of the original data if \code{newdata} is not provided).
#'
#' If \code{return.calib = TRUE}, \code{return.par = TRUE}
#' and/or \code{return.tau = TRUE}, then the list also contains respectively items
#' \code{calib}, \code{par} and/or \code{tau}.
#' @examples
#' require(copula)
#' set.seed(1)
#'
#' ## Simulation parameters (sample size, correlation between covariates,
#' ## Gaussian copula family)
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
#' sim <- gamBiCopSimulate(fit$res, X)
#' @export
gamBiCopSimulate <- function(object, newdata = NULL, N = NULL, return.calib = FALSE,
                             return.par = FALSE, return.tau = FALSE) {
  tmp <- valid.gamBiCopSimulate(
    object, newdata, N, return.calib,
    return.par, return.tau
  )
  if (tmp != TRUE) {
    stop(tmp)
  }


  if (is.null(newdata)) {
    newdata <- object@model$data
  }

  if (is.null(N)) {
    if (!is.null(newdata)) {
      N <- dim(newdata)[1]
    } else {
      N <- dim(object@model$data)[1]
    }
  }

  dd <- dim(newdata)[1]
  if ((N < dd) || (N > dd)) {
    newdata <- newdata[sample.int(dd, N, replace = TRUE), ]
  }

  tmp <- gamBiCopPredict(object, as.data.frame(newdata),
    target = c("par", "calib", "tau")
  )
  par <- tmp$par
  family <- object@family
  if (family %in% c(1, 2, 5)) {
    family <- rep(family, length(par))
  } else {
    fam <- getFams(family)
    sel <- par > 0
    family <- rep(0, length(par))
    family[sel] <- fam[1]
    family[!sel] <- fam[2]
  }

  data <- t(mapply(BiCopSim, 1, family, par, object@par2))
  out <- list(data = data)

  if (return.calib == TRUE) {
    out$calib <- tmp$calib
  }

  if (return.par == TRUE) {
    out$par <- tmp$par
  }

  if (return.tau == TRUE) {
    out$tau <- tmp$tau
  }

  return(out)
}

valid.gamBiCopSimulate <- function(object, newdata, N, return.calib,
                                   return.par, return.tau) {
  if (!valid.gamBiCop(object)) {
    return("gamBiCopPredict can only be used to predict from gamBiCop objects")
  }

  if (is.null(N)) {
    if (!is.null(newdata)) {
      N <- dim(newdata)[1]
    } else {
      N <- dim(object@model$data)[1]
    }
  }
  if (!valid.posint(N)) {
    return(msg.posint(var2char(N)))
  }

  if (!valid.logical(return.calib)) {
    return(msg.logical(var2char(return.calib)))
  }

  if (!valid.logical(return.par)) {
    return(msg.logical(var2char(return.par)))
  }

  if (!valid.logical(return.tau)) {
    return(msg.logical(var2char(return.tau)))
  }

  return(TRUE)
}
