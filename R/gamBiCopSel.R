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
#' Vector of bivariate copula families to select from. 
#' If \code{familyset = NA} (default), selection
#' among all possible families is performed.   Coding of bivariate copula 
#' families:
#' \code{1} Gaussian, 
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
#' @param rotations If \code{TRUE}, all rotations of the families in familyset 
#' are included.
#' @param selcrit Character indicating the criterion for bivariate copula 
#' selection. Possible choices: \code{selcrit = 'AIC'} (default) or 
#' \code{'BIC'}, as in \code{\link{BiCopSelect}} from the 
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package. 
#' @param level Numerical; significance level of the test for removing individual
#' predictors (default: \code{level = 0.05}).  
#' @param edf Numerical; if the estimated EDF for individual predictors is 
#' smaller than \code{edf} but the predictor is still significant, then
#' it is set as linear (default: \code{edf = 1.5}). 
#' @param tau \code{FALSE} for a calibration fonction specified for 
#' the Copula parameter or \code{TRUE} (default) for a calibration function 
#' specified for Kendall's tau.  
#' @param method \code{'FS'} for Fisher-scoring (default) and 
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
#' require(copula)
#' set.seed(0)
#' 
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
#' covariates.distr <- mvdc(normalCopula(rho, dim = 6),
#'                                  c("unif"), list(list(min = 0, max = 1)),
#'                                  marginsIdentical = TRUE)
#' X <- rMvdc(n, covariates.distr)
#' 
#' ## U in [0,1]x[0,1] depending on the four first columns of X
#' U <- condBiCopSim(fam, function(x1,x2,x3,x4) {eta0+sum(mapply(function(f,x)
#'   f(x), calib.surf, c(x1,x2,x3,x4)))}, X[,1:4], par2 = 4, return.par = TRUE)
#' 
#' ## Merge U and X
#' data <- data.frame(U$data,X)
#' names(data) <- c(paste("u",1:2,sep=""),paste("x",1:6,sep=""))
#' 
#' ## Selection using AIC (take about 5mn on single core) 
#' ## Use parallel = TRUE to speed-up....
#' system.time(best <- gamBiCopSel(data))
#' print(best$res)
#' EDF(best$res)
#' plot(best$res)
#' @export
gamBiCopSel <- function(data, familyset = NA, rotations = TRUE, 
                        selcrit = "AIC", level = 5e-2, edf = 1.5, tau = TRUE, 
                        method = "FS", tol.rel = 1e-3, n.iters = 10,
                        parallel = FALSE, verbose = FALSE, ...) {

  tmp <- valid.gamBiCopSel(data, rotations, familyset, selcrit, level, edf, tau, 
                           method, tol.rel, n.iters, parallel, verbose)
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
  
  if (length(familyset) == 1 && is.na(familyset)) {
    familyset <- get.familyset()
  }
  if (rotations) {
    familyset <- withRotations(familyset)
  }

  parallel <- as.logical(parallel)
  x <- NULL
  if (parallel == FALSE) {
    res <- foreach(x=familyset) %do% 
      tryCatch(gamBiCopVarSel(data, x, tau, method, tol.rel, n.iters, 
                              level, edf, verbose,...),
               error = function(e) e)
  } else {
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl, cores = detectCores() - 1)
    res <- foreach(x=familyset) %dopar% 
      tryCatch(gamBiCopVarSel(data, x, tau, method, tol.rel, n.iters, 
                              level, edf, verbose,...),
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

  if (selcrit == "AIC") {
    return(res[[which.min(sapply(res, function(x) AIC(x$res)))]])
  } else {
    return(res[[which.min(sapply(res, function(x) BIC(x$res)))]])
  }
} 

gamBiCopVarSel <- function(data, family,
                        tau = TRUE, method = "FS",
                        tol.rel = 1e-3, n.iters = 10, 
                        level = 5e-2, edf = 1.5, 
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
    sel <- summary(tmp$res@model)$s.pv < level
    nn <- nn[sel]
    basis <- rep(5,length(nn))
    formula.expr <- mapply(get.formula,nn,basis)
  }
  
  if (length(basis) == 0) {
    tmp <- gamBiCopEst(data, ~1, family, tau, 
                       method, tol.rel, n.iters)
    return(tmp)
  }

  ## Create a separate list by setting as linear the predictors with EDF < edf
  formula.tmp <- as.formula(paste("~",paste(formula.expr,collapse = " + ")))
  tmp <- gamBiCopEst(data, formula.tmp, family, tau, 
                     method, tol.rel, n.iters)
  sel <- summary(tmp$res@model)$edf < edf
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
    cat("Remove by setting as linear the predictors with EDF < edf....... \n")
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
    if (any(sel)) {
      if (sum(sel) == 1) {
        sel[sel] <- 2*basis[sel] < length(unique(data[,nn[sel]]))
      } else {
        sel[sel] <- 2*basis[sel] < apply(data[,nn[sel]], 2, function(x) 
          length(unique(x)))
      }
      
    }
  }
  return(tmp)
}

valid.gamBiCopSel <- function(data, rotations, familyset, selcrit, level, edf, 
                              tau, method, tol.rel, n.iters, parallel, 
                              verbose) {
  
  tmp <- valid.gamBiCopEst(data, n.iters, tau, tol.rel, method, verbose, 1)
  if (tmp != TRUE) {
    return(tmp)
  }
  
  if (!valid.familyset(familyset)) {
    return(return(msg.familyset(var2char(familyset))))
  }
  
  if (!valid.logical(rotations)) {
    return(msg.logical(var2char(rotations)))
  }
  
  if(is.null(selcrit) || length(selcrit) != 1 || 
       (selcrit != "AIC" && selcrit != "BIC")) {
    return("Selection criterion not implemented.")
  } 
  
  if (!valid.logical(parallel)) {
    return(msg.logical(var2char(parallel)))
  }
  
  if (!valid.unif(level)) {
    return(msg.unif(var2char(level)))
  }
  
  if (!valid.real(edf)) {
    return(msg.real(var2char(edf)))
  }
  
  return(TRUE)
}
