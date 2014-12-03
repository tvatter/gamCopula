############################# gam bivariate copulas ##
valid.gamBiCop <- function(object) {
  d <- length(attributes(object))
  if (d < 2) {
    return("A gamBiCop contains at least a copula family and a mgcv model.")
  } else if (d >= 2) {
    if (!is(object@model, "gam")) {
      return("Invalid mgcv model.")
    } else if (!(object@family %in% c(1, 2, 3, 4, 13, 14, 23, 24, 33, 34))) {
      return("Copula family not yet implemented.")
    }
  }
  if (!(is.logical(object@tau) || (object@tau == 0) || (object@tau == 1))) {
    return("Tau should takes 0/1 or FALSE/TRUE to model the copula parameter/tau's tau.")
  }
  if ((object@family == 2) && (is.null(object@par2) || is.na(as.numeric(object@par2)) || 
    as.numeric(object@par2) <= 2)) {
    return("Par2 greater than 2 is needed for the t-copula.")
  }
  return(TRUE)
}

setOldClass(c("gam"))

#'  The \code{\link{gamBiCop-class}}
#'
#'  \code{\link{gamBiCop-class}} is an S4 class to store 
#'  a Generalized Additive Model for bivariate copula a parameter or Kendall's tau.
#'
#' @slot family A copula family: \code{1} Gaussian, \code{2} Student t, 
#' \code{3} Clayton, \code{4} Gumbel, \code{13} Survival Clayton, \code{14} Survival Gumbel, 
#' \code{23} Rotated (90 degrees) Clayton, \code{24} Rotated (90 degrees) Gumbel, 
#' \code{33} Rotated (270 degrees) Clayton and \code{34} Rotated (270 degrees) Gumbel.
#' @slot model A \code{\link{gamObject}} as return by the \code{\link{gam}} function 
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#' @slot par2 Second parameter for the Studen t-copula.
#' @slot tau \code{FALSE} (default) for a calibration fonction specified for the Copula parameter 
#' or \code{TRUE} for a calibration function specified for Kendall's tau.
#' @seealso \code{\link{gamBiCop}}, \code{\link{gamBiCopEst}} and \code{\link{gamBiCopPred}}.
#' @export
setClass("gamBiCop", slots = c(family = "integer", model = "gam", par2 = "numeric", 
  tau = "logical"), validity = valid.gamBiCop)

#' Constructor of the \code{\link{gamBiCop-class}}
#' 
#' A constructor for objects of the \code{\link{gamBiCop-class}}.
#'
#' @param family A copula family: \code{1} Gaussian, \code{2} Student t, 
#' \code{3} Clayton, \code{4} Gumbel, \code{13} Survival Clayton, \code{14} Survival Gumbel, 
#' \code{23} Rotated (90 degrees) Clayton, \code{24} Rotated (90 degrees) Gumbel, 
#' \code{33} Rotated (270 degrees) Clayton and \code{34} Rotated (270 degrees) Gumbel.
#' @param model A \code{\link{gamObject}} as return by the \code{\link{gam}} function 
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#' @param par2 Second parameter for the Studen t-copula.
#' @param tau \code{FALSE} (default) for a calibration fonction specified for the Copula parameter 
#' or \code{TRUE} for a calibration function specified for Kendall's tau.
#' @return A \code{\link{gamBiCop-class}} object.
#' @seealso \code{\link{gamBiCop-class}}, \code{\link{gamBiCopEst}} and \code{\link{gamBiCopPred}}.
#' @export
gamBiCop <- function(family, model, par2 = 0, tau = FALSE) {
  if (family != 2) {
    par2 <- 0
  }
  if (as.integer(family) != family) {
    return("Family should be an integer.")
  }
  if (!(is.logical(tau) || (tau == 0) || (tau == 1))) {
    return("Tau should takes 0/1 or FALSE/TRUE to model the copula parameter/Kendall's tau.")
  }
  new("gamBiCop", family = as.integer(family), model = model, par2 = par2, tau = as.logical(tau))
} 
