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
#' \code{3} Clayton, 
#' \code{4} Gumbel, 
#' \code{13} Survival Clayton, 
#' \code{14} Survival Gumbel, 
#' \code{23} Rotated (90 degrees) Clayton, 
#' \code{24} Rotated (90 degrees) Gumbel, 
#' \code{33} Rotated (270 degrees) Clayton and 
#' \code{34} Rotated (270 degrees) Gumbel.
#' @slot model A \code{\link{gamObject}} as return by the
#'  \code{\link{gam}} function 
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#' @slot par2 Second parameter for the Studen t-copula.
#' @slot tau \code{FALSE} (default) for a calibration fonction 
#' specified for the Copula parameter 
#' or \code{TRUE} for a calibration function specified for Kendall's tau.
#' @seealso \code{\link{gamBiCopEst}}, 
#' \code{\link{gamBiCopPred}} and \code{\link{gamBiCopSim}}.
#' @name gamBiCop-class
#' @rdname gamBiCop-class
#' @export
setClass("gamBiCop", slots = c(family = "integer", model = "gam", 
                               par2 = "numeric", tau = "logical"))

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
#' @seealso \code{\link{gamVine}}, 
#' \code{\link{RVineMatrix}}, \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}
#' \code{\link{gamVineSeqEst}}, \code{\link{gamVineStructureSelect}}, 
#' \code{\link{gamVinePred}} and \code{\link{gamVineSim}}.
#' @name gamVine-class
#' @rdname gamVine-class
#' @export
setClass("gamVine", slots = c(Matrix = "matrix", model = "list", 
                              names = "character"))