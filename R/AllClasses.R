setOldClass(c("gam"))

#'  The gamBiCop Class
#'
#'  gamBiCop is an S4 class to store
#'  a Generalized Additive Model for bivariate copula a parameter or
#'  Kendall's tau. Objects can be created by calls of the form
#'  \code{new("gamBiCop", ...)}, or by function \code{\link{gamBiCop}}.
#'
#' @slot family A copula family: \code{1} Gaussian,
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
#' @slot model A \code{\link[mgcv]{gamObject}} as return by the
#'  \code{\link[mgcv]{gam}} function
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#' @slot par2 Second parameter for the Studen t-copula.
#' @slot tau \code{FALSE} (default) for a calibration fonction
#' specified for the Copula parameter
#' or \code{TRUE} for a calibration function specified for Kendall's tau.
#' @seealso \code{\link{gamBiCopFit}},
#' \code{\link{gamBiCopPredict}} and \code{\link{gamBiCopSimulate}}.
#' @name gamBiCop-class
#' @rdname gamBiCop-class
#' @export
setClass("gamBiCop", slots = c(
  family = "integer", model = "gam",
  par2 = "numeric", tau = "logical"
))

#'  The gamVine Class
#'
#'  gamVine is an S4 class to store a conditional and potentially
#'  non-simplified pair-copula construction. Objects can be created by calls of
#'  the form \code{new("gamVine", ...)}, or by function \code{\link{gamVine}}.
#'
#' @slot Matrix Lower triangular d x d matrix that defines the tree structure.
#' @slot model list containing d x (d-1)/2 lists with three numeric items
#' (family, par and par2) and/or
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}} objects.
#' @slot names vector of d names.
#' @slot covariates vector of names for the exogenous covariates.
#' @seealso \code{\link{gamVine}},
#' \code{\link[VineCopula]{RVineMatrix}},
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}
#' \code{\link{gamVineSeqFit}}, \code{\link{gamVineCopSelect}},
#' \code{\link{gamVineStructureSelect}} and \code{\link{gamVineSimulate}}.
#' @name gamVine-class
#' @rdname gamVine-class
#' @export
setClass("gamVine", slots = c(
  Matrix = "matrix", model = "list",
  names = "character", covariates = "character"
))
