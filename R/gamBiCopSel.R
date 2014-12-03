#' Selection and Maximum penalized likelihood estimation of a Generalized 
#' Additive model (gam) for the copula parameter or Kendall's tau.
#'
#' This function selects an appropriate bivariate copula family for given 
#' bivariate copula data using one of a range of methods. The corresponding 
#' parameter estimates are obtained by maximum penalized  likelihood estimation,
#' where each Newton-Raphson iteration is reformulated as a generalized ridge 
#' regression solved using the \code{\link[mgcv:mgcv-package]{mgcv}} package. 
#' 
#' @param data A list or data frame containing the model responses, (u1,u2) in 
#' [0,1]x[0,1], and covariates required by the formula.
#' @param familyset (Similar to \code{\link{BiCopSelect}} from the 
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package) 
#' Vector of bivariate copula families to select from. The vector has to include
#' at least one bivariate copula family that allows for positive and one that 
#' allows for negative dependence. If \code{familyset = NA} (default), selection
#' among all possible families is performed.   Coding of bivariate copula 
#' families: \code{1} Gaussian, \code{2} Student t, \code{3} Clayton, 
#' \code{4}Gumbel, \code{13} Survival Clayton, \code{14} Survival Gumbel, 
#' \code{23} Rotated (90 degrees) Clayton, \code{24} Rotated (90 degrees) Gumbel
#' , \code{33} Rotated (270 degrees) Clayton and \code{34} Rotated (270 degrees)
#' Gumbel.  
#' @param selectioncrit Character indicating the criterion for bivariate copula 
#' selection. Possible choices: \code{selectioncrit = 'AIC'} (default) or 
#' \code{'BIC'}, as in \code{\link{BiCopSelect}} from the 
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package.   
#' @param tau \code{FALSE} (default) for a calibration fonction specified for 
#' the Copula parameter or \code{TRUE} for a calibration function specified for 
#' Kendall's tau.  
#' @param method \code{'FS'} for Fisher-scoring and 
#' \code{'NR'} for Newton-Raphson.  
#' @param tol.rel Relative tolerance for \code{'FS'}/\code{'NR'} algorithm.  
#' @param n.iters Maximal number of iterations for 
#' \code{'FS'}/\code{'NR'} algorithm.
#' @param parallel \code{TRUE} for a parallel estimation across copula families.
#' As the code is based on mclapply, this parameter has no effect on windows.
#' @param verbose \code{TRUE} prints informations during the estimation.
#' @param ... Additional parameters to be passed to \code{\link{gam}} 
#' @return \code{gamBiCopEst} returns a list consisting of 
#' \item{res}{S4 \code{\link{gamBiCop-class}} object.} 
#' \item{method}{\code{'FS'} for Fisher-scoring and 
#' \code{'NR'} for Newton-Raphson.} 
#' \item{tol.rel}{relative tolerance for \code{'FS'}/\code{'NR'} algorithm.} 
#' \item{n.iters}{maximal number of iterations for 
#' \code{'FS'}/\code{'NR'} algorithm.} 
#' \item{trace}{the estimation procedure's trace.} 
#' \item{conv}{\code{0} if the algorithm converged and \code{1} otherwise.}
#' @seealso \code{\link{gamBiCop}} and \code{\link{gamBiCopEst}}. 
#' @examples
#' ## Simulation parameters (sample size, correlation between covariates,
#' ## Student copula with 4 degrees of freedom)
#' n <- 5e2
#' rho <- 0.9
#' fam <- 2
#' par2 <- 4
#' 
#' ## A calibration surface depending on four variables
#' eta0 <- 1
#' calib.surf <- list(calib.lin <- function(t, Ti = 0, Tf = 1, b = 2) {
#'     return(-2+4*t)},
#'   calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
#'     Tm <- (Tf - Ti)/2
#'     a <- -(b/3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
#'   return(a + b * (t - Tm)^2)},
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
#' ## 6-dimensional matrix X of covariates
#' covariates.distr <- copula::mvdc(copula::normalCopula(rho, dim = 6),
#'                                  c("unif"), list(list(min = 0, max = 1)),
#'                                  marginsIdentical = TRUE)
#' X <- copula::rMvdc(n, covariates.distr)
#' 
#' ## U in [0,1]x[0,1] depending on the four first columns of X
#' U <- CondBiCopSim(fam, function(x1,x2,x3,x4) {eta0+sum(mapply(function(f,x)
#'   f(x), calib.surf, c(x1,x2,x3,x4)))}, X[,1:4], par2 = 4, return.par = TRUE)
#' 
#' ## Merge U and X
#' data <- data.frame(U$data,X)
#' names(data) <- c(paste("u",1:2,sep=""),paste("x",1:6,sep=""))
#' 
#' ## Parallel family selection using AIC
#' system.time(best <- gamBiCopSel(data, parallel = TRUE))
#' print(best$res)
#' EDF(best$res)
#' @export
gamBiCopSel <- function(data, familyset = NA, selectioncrit = "AIC",
                        tau = FALSE, method = "FS", 
                        tol.rel = 1e-3, n.iters = 10, 
                        parallel = FALSE, verbose = FALSE, ...) {
  
  if (!(is.list(data) || is.data.frame(data))) {
    stop("data has to be either a list or a data frame")
  } else {
    if(is.list(data)){
      if(!is.null(data$xt)){
        xt <- data$xt
        data <- as.data.frame(data[-which(names(data) == "xt")])
      }else{
        data <- as.data.frame(data)
      }     
    }
    n <- dim(data)[1]
    m <- dim(data)[2]
    u1 <- data$u1
    u2 <- data$u2
    u <- cbind(u1,u2)
  }
  
  if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
    stop("u1 and/or u2 are missing.")
  if (n < 2) 
    stop("Number of observations has to be at least 2.")
  if (any(u > 1) || any(u < 0)) 
    stop("(u1, u2) have to be in [0,1]x[0,1].")
 
  k.tau <- VineCopula:::fasttau(u1,u2)

  if (is.null(familyset) || (any(is.na(familyset)) && length(familyset) != 1) || 
        any(sapply(familyset, function(x) 
          !(x %in% c(NA,1,2,3,4,13,14,23,24,33,34))))) {
    stop("Copula family not implemented.")
  } else if (length(familyset) == 1 && is.na(familyset)) {
    if ( k.tau < 0) {
      familyset <- c(1,2,23,24,33,34)
    } else { 
      familyset <- c(1,2,3,4,13,14)
    }  
  } else if ( k.tau < 0 && !any(c(1,2,23,24,33,34) %in% familyset) ) {
    stop(paste("Because Kendall's tau is negative, familyset needs at least ",
               "one bivariate copula family for negative dependence."))
  } else if ( k.tau > 0 && !any(c(1,2,3,4,13,14) %in% familyset) ) { 
    stop(paste("Because Kendall's tau is positive, familyset needs at least ",
               "one bivariate copula family for positive dependence."))
  } else if ( k.tau < 0 && any(c(1,2,3,4,13,14) %in% familyset) ) {
    stop(paste("Because Kendall's tau is negative, familyset cannot contain ",
               "bivariate copula families for positive dependence."))
  } else if ( k.tau > 0 && any(c(1,2,23,24,33,34) %in% familyset) ) { 
    stop(paste("Because Kendall's tau is positive, familyset cannot contain ",
               "bivariate copula families for negative dependence."))
  } else {
    familyset <- as.integer(familyset)
  }

  options(warn = -1)
  if (is.null(n.iters) || is.na(as.integer(n.iters)) || 
        (as.integer(n.iters) < 1) || 
        (as.integer(n.iters) !=  as.numeric(n.iters))) {
    stop("N.iters should be a positive integer.")
  } else {
    n.iters <- as.integer(n.iters)
  }
  
  if (is.null(tau) ||  is.na(tau) ||  
        !(is.logical(tau) || (tau == 0) || (tau == 1))) {
    stop("Tau should takes 0/1 or FALSE/TRUE to specify a model for 
         the copula parameter/Kendall's tau.")
  }
  
  if (is.null(tol.rel) ||  is.na(as.numeric(tol.rel)) || 
        (as.numeric(tol.rel) < 0) || (as.numeric(tol.rel) > 1)) {
    stop("Tol.rel should be a real number in [0,1].")
  } else {
    tol.rel <- as.numeric(tol.rel)
  }
  
  if (is.null(method) || !is.element(method, c("FS", "NR"))) {
    stop("Method should be a string, either NR (Newton-Raphson) 
         or FS (Fisher-scoring, faster but unstable).")
  }
  
  if (is.null(verbose) ||  is.na(verbose) || 
        !(is.logical(verbose) || (verbose == 0) || (verbose == 1))) {
    stop("Verbose should takes 0/1 or FALSE/TRUE.")
  } else {
    verbose <- as.logical(verbose)
  }
  
  if (is.null(parallel) ||  is.na(parallel) || 
        !(is.logical(parallel) || (parallel == 0) || (parallel == 1))) {
    stop("parallel should takes 0/1 or FALSE/TRUE.")
  } else {
    parallel <- as.logical(parallel)
  }
  options(warn = 0)
  
  if (parallel == FALSE) {
    res <- lapply(familyset,function(x) 
      gamBiCopVarSel(data,x,tau,method,tol.rel,n.iters,verbose,...))
  } else {
    res <- parallel::mclapply(familyset,function(x) 
      gamBiCopVarSel(data,x,tau,method,tol.rel,n.iters,verbose,...))
  }
  
  if (selectioncrit == "AIC") {
    return(res[[which.min(sapply(res, function(x) AIC(x$res)))]])
  } else {
    return(res[[which.min(sapply(res, function(x) BIC(x$res)))]])
  }
} 

gamBiCopVarSel <- function(data, family,
                        tau = FALSE, method = "FS", 
                        tol.rel = 1e-3, n.iters = 10, 
                        verbose = FALSE, ...) {
  
  ## Create a list with formulas of smooth terms corresponding to the covariates
  get.formula <- function(x,k){
    paste("s(",x,", k=",k,", bs='cr')",sep = "")
  } 
  nn <- names(data[,-which( (names(data) == "u1") | (names(data) == "u2"))])
  basis <- rep(5,ncol(data)-2)
  formula.expr <- mapply(get.formula,nn,basis)

  ## Update the list by removing unsignificant predictors 
  sel <- FALSE
  while(!all(sel)){
    formula.tmp <- as.formula(paste("~",paste(formula.expr,collapse = " + ")))
    tmp <- gamBiCopEst(data, formula.tmp, family, tau, 
                       method, tol.rel, n.iters)
    sel <- summary(tmp$res@model)$s.pv < 5e-2
    nn <- nn[sel]
    basis <- rep(5,length(nn))
    formula.expr <- mapply(get.formula,nn,basis)
  }

  ## Create a separate list by setting as linear the predictors with EDF < 1.5
  formula.tmp <- as.formula(paste("~",paste(formula.expr,collapse = " + ")))
  tmp <- gamBiCopEst(data, formula.tmp, family, tau, 
                     method, tol.rel, n.iters)
  sel <- summary(tmp$res@model)$edf < 1.5
  get.linear <- function(x){
    names(unlist(sapply(nn,function(z)grep(z,x))))
  }
  formula.lin <- get.linear(formula.expr[sel])
  formula.expr <- formula.expr[!sel]
  basis <- basis[!sel]
  nn <- nn[!sel]
  formula.tmp <- as.formula(paste("~",paste(c(formula.lin,
                                        formula.expr),
                                        collapse = " + ")))
  tmp <- gamBiCopEst(data, formula.tmp, family, tau, 
                     method, tol.rel, n.iters)

  ## Increasing the basis size appropriately
  sel <- summary(tmp$res@model)$edf > (basis-1)/2 
  while(any(sel)){
    basis[sel] <- 2*basis[sel]
    formula.expr[sel] <- mapply(get.formula,nn[sel],basis[sel])
    formula.tmp <- as.formula(paste("~",paste(c(formula.lin,
                                              formula.expr),
                                              collapse = " + ")))
    tmp <- gamBiCopEst(data, formula.tmp, family, tau, 
                       method, tol.rel, n.iters)
    sel <- summary(tmp$res@model)$edf > (basis-1)/2
  }
  return(tmp)
}