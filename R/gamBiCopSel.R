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
#' families:
#' \code{1} Gaussian, 
#' \code{2} Student t, 
#' \code{3} Clayton, 
#' \code{4} Gumbel,
#' \code{5} Frank, 
#' \code{13} Survival Clayton, 
#' \code{14} Survival Gumbel,  
#' \code{23} Rotated (90 degrees) Clayton, 
#' \code{24} Rotated (90 degrees) Gumbel, 
#' \code{33} Rotated (270 degrees) Clayton and 
#' \code{34} Rotated (270 degrees) Gumbel.
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
#' U <- condBiCopSim(fam, function(x1,x2,x3,x4) {eta0+sum(mapply(function(f,x)
#'   f(x), calib.surf, c(x1,x2,x3,x4)))}, X[,1:4], par2 = 4, return.par = TRUE)
#' 
#' ## Merge U and X
#' data <- data.frame(U$data,X)
#' names(data) <- c(paste("u",1:2,sep=""),paste("x",1:6,sep=""))
#' 
#' ## Selection using AIC (take about 3mn on single core) 
#' ## Use parallel = TRUE to speed-up.... currently unavailable for windows 
#' ## users as the parallelization is handled by mclapply!
#' system.time(best <- gamBiCopSel(data))
#' print(best$res)
#' EDF(best$res)
#' @export
gamBiCopSel <- function(data, familyset = NA, selectioncrit = "AIC",
                        tau = FALSE, method = "FS", 
                        tol.rel = 1e-3, n.iters = 10, 
                        parallel = FALSE, verbose = FALSE, ...) {

  tmp <- valid.gamBiCopSel(data, familyset, selectioncrit, tau, method, tol.rel, 
                           n.iters, parallel, verbose)
  if (tmp != TRUE)
    stop(tmp)
  
  if (is.list(data)){
    if(!is.null(data$xt)){
      xt <- data$xt
      data <- data[-which(names(data) == "xt")]
    }
    data <- as.data.frame(data)
  }
  
  n <- dim(data)[1]
  m <- dim(data)[2]
  u1 <- data$u1
  u2 <- data$u2
  u <- cbind(u1,u2)

  k.tau <- fasttau(u1,u2)
  
  if (length(familyset) == 1 && is.na(familyset)) {
    if ( k.tau < 0) {
      familyset <- c(1,2,23,24,33,34)
    } else { 
      familyset <- c(1,2,3,4,13,14)
    }  
  }
  
  ## find families for which estimation is required (only families that allow for
  ## the empirical kendall's tau)
  if ( k.tau < 0 ) {
    todo <- c(1, 2, 5, 23, 24, 26:30, 33, 34, 36:40, 124, 134, 224, 234)
  } else {
    todo <- c(1:10, 13, 14, 16:20, 104, 114, 204, 214)
  }
  familyset <- as.integer(todo[which(todo %in% familyset)])


  parallel <- as.logical(parallel)
  if (parallel == FALSE) {
    res <- foreach(x=familyset) %do% 
      tryCatch(gamBiCopVarSel(data,x,tau,method,tol.rel,n.iters,verbose,...),
               error = function(e) e)
  } else {
    cl <- makeCluster(parallel::detectCores() - 1)
    registerDoParallel(cl, cores = detectCores() - 1)
    res <- foreach(x=familyset) %dopar% 
      tryCatch(gamBiCopVarSel(data,x,tau,method,tol.rel,n.iters,verbose,...),
               error = function(e) e)
    stopCluster(cl)
  }
  res <- res[sapply(res,length) == 6] 
  if (length(res) == 0 ) { #|| sum(sapply(res, function(x) x$conv) == 0) == 0) {
    return(paste("No convergence of the estimation for any copula family.",
                 "Try modifying the parameters (FS/NR, tol.rel, n.iters)..."))
  }
#   } else {
#     res <- res[sapply(res, function(x) x$conv) == 0]
#   }

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
  
  if (verbose == TRUE) {
    cat(paste("Model selection for family", family, "\n"))
  }
  
  n <- dim(data)[1]

  ## Create a list with formulas of smooth terms corresponding to the covariates
  get.formula <- function(x,k){
    paste("s(",x,", k=",k,", bs='cr')",sep = "")
  } 
  nn <- names(data)[-which( (names(data) == "u1") | (names(data) == "u2"))]
  basis <- rep(5,ncol(data)-2)
  formula.expr <- mapply(get.formula,nn,basis)

  ## Update the list by removing unsignificant predictors 
  if (verbose == TRUE) {
    cat("Remove unsignificant covariates.......\n")
  }
  sel <- FALSE
  while(!all(sel) && length(basis) > 0){
    formula.tmp <- as.formula(paste("~",paste(formula.expr,collapse = " + ")))
    if (verbose == TRUE) {
      cat("Model formula:\n")
      print(formula.tmp)
    }
    tmp <- gamBiCopEst(data, formula.tmp, family, tau, 
                       method, tol.rel, n.iters)
    sel <- summary(tmp$res@model)$s.pv < 5e-2
    nn <- nn[sel]
    basis <- rep(5,length(nn))
    formula.expr <- mapply(get.formula,nn,basis)
  }
  
  if (length(basis) == 0) {
    tmp <- gamBiCopEst(data, ~1, family, tau, 
                       method, tol.rel, n.iters)
    return(tmp)
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
  if (verbose == TRUE) {
    cat("Remove by setting as linear the predictors with EDF < 1.5....... \n")
    cat("Updated model formula:\n")
    print(formula.tmp)
  }
  tmp <- gamBiCopEst(data, formula.tmp, family, tau, 
                     method, tol.rel, n.iters)

  ## Increasing the basis size appropriately
  sel <- summary(tmp$res@model)$edf > (basis-1)/2 
  if (verbose == TRUE) {
    cat(paste("For the other predictors,", 
        "increase the basis size appropriately.......\n"))
  }
  while (any(sel) && all(basis < n/30)) {
    basis[sel] <- 2*basis[sel]
    #browser()
    formula.expr[sel] <- mapply(get.formula,nn[sel],basis[sel])
    formula.tmp <- as.formula(paste("~",paste(c(formula.lin,
                                              formula.expr),
                                              collapse = " + ")))
    if (verbose == TRUE) {
      cat("Updated model formula:\n")
      print(formula.tmp)
    }
    tmp <- gamBiCopEst(data, formula.tmp, family, tau, 
                       method, tol.rel, n.iters)
    sel <- summary(tmp$res@model)$edf > (basis-1)/2
  }
  return(tmp)
}

valid.gamBiCopSel <- function(data, familyset, selectioncrit, tau, method, 
                              tol.rel, n.iters, parallel, verbose) {
  
  tmp <- valid.gamBiCopEst(data, n.iters, tau, tol.rel, method, verbose, 1)
  if (tmp != TRUE) {
    return(tmp)
  }
  
  if (!valid.familyset) {
    return(return(msg.familyset(var2char(familyset))))
  }
  
  if (is.list(data)){
    if(!is.null(data$xt)){
      xt <- data$xt
      data <- data[-which(names(data) == "xt")]
    }
    data <- as.data.frame(data)
  }
  
  n <- dim(data)[1]
  m <- dim(data)[2]
  u1 <- data$u1
  u2 <- data$u2
  u <- cbind(u1,u2)
  
  tau <- fasttau(u1,u2) 
  if (!valid.familysetpos(familyset, tau)) {
    return(paste("Because Kendall's tau is positive", 
                 msg.familysetpos(var2char(familyset))))
  }
  
  if (!valid.familysetneg(familyset, tau)) {
    return(paste("Because Kendall's tau is negative,",
                 msg.familysetneg(var2char(familyset))))
  }
  
  if(is.null(selectioncrit) || length(selectioncrit) != 1 || 
       (selectioncrit != "AIC" && selectioncrit != "BIC")) {
    return("Selection criterion not implemented.")
  } 
  
  if (!valid.logical(parallel)) {
    return(msg.logical(var2char(parallel)))
  }
  options(warn = 0)
  
  
  return(TRUE)
}
