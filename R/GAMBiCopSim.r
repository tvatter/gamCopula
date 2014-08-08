#' Simulate from \code{\link{gamBiCop-class}} object
#'
#' @param object \code{\link{gamBiCop-class}} object.
#' @param newdata (same as in \code{\link{predict.gam}} from the \code{\link[mgcv:mgcv-package]{mgcv}} package) A matrix or data frame containing the values of the model covariates at which simulations are required. 
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
#' ##     												
#' ##	1) Simulate some data	
#' ##  a) Kendall's tau as a sum of three smooth components
#' ##  b) Distribution of the covariates as a Gaussian copula 
#' ##	2) Estimate the model and simulate new data							
#' 
#' require(gamVineCopula)
#' 
#' set.seed(0)
#' 
#' ##' Simulation parameters
#' # Sample size
#' n <- 2e2  
#' # Correlation between the covariates
#' rho <- 0.5
#' # Copula family (Clayton here)
#' fam <- 3  
#' # Degrees of freedom (for the t-copula when fam <- 2)
#' par2 <- 4 
#' # Should the model be specified in terms of Kendall's tau (TRUE) or copula parameter
#' tau <- FALSE
#' # Newton-Raphse ("NR) or Fisher-Scoring ("FS") algorithm
#' met <- "FS"
#' # Relative tolerance for "NR"/"FS"
#' tol <- 1e-6
#' # Max number of iterations for "NR"/"FS"
#' itermax <- 25
#' 
#' ##
#' ## 1) Simulate some data
#' ##
#' 
#' ## Quadratic calibration
#' Tf <- 1
#' b.quad <- 8*Tf
#' t0.quad <- Tf/2
#' a.quad <- -(b.quad/3)*(Tf^2-3*Tf*t0.quad+3*t0.quad^2)
#' calib.quadratic <- function(t, t0, a, b){return(a + b*(t-t0)^2)}
#' 
#' ## Sinusoidal calibration
#' b.sin <- 1
#' f.sin <- 1
#' t0.sin <- 0
#' a.sin <- b.sin*(1-2*f.sin*Tf*pi/(f.sin*Tf*pi+cos(2*f.sin*pi*(Tf-t0.sin))-cos(2*f.sin*pi*t0.sin)))
#' calib.sinusoidal <- function(t, f, t0, a, b){return((a+b)/2 + (b-a)*sin(2*pi*f*(t-t0))/2)}
#' 
#' ## Exponential calibration 
#' t0.exp  <- Tf/2
#' s.exp <- Tf/8
#' b.exp <- 2
#' a.exp <- (b.exp*s.exp*sqrt(2*pi)/Tf)*(pnorm(0,t0.exp,s.exp)-pnorm(Tf,t0.exp,s.exp))
#' calib.exponential <- function(t, a, b, t0, s){return(a+b*exp(-(t-t0)^2/(2*s^2)))}
#' 
#' ## Base level
#' eta0 <- 1
#' 
#' ## Additive calibration function
#' calib.fct <- function(x1,x2,x3){
#'   return(eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#'          calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin)+
#'          calib.exponential(x3, a.exp, b.exp, t0.exp, s.exp))}
#' 
#' 
#' ## A single dataset
#' covariates.distr <- copula::mvdc(copula::normalCopula(rho, dim = 3), 
#'    c("unif"), list(list(min = 0, max = Tf)), marginsIdentical = TRUE) 
#' X <- copula::rMvdc(n,covariates.distr)
#' temp <- CondBiCopSim(fam, calib.fct, X, par2=par2, return.par=TRUE)
#' calib <- temp$eta
#' param <- temp$par
#' pseudo <- temp$data
#' dataset <- data.frame("u1" = pseudo[,1], "u2" = pseudo[,2], 
#'    "x1" = X[,1], "x2" = X[,2], "x3" = X[,3])
# 
#' ##
#' ## 2) Estimate and simulate model
#' ##
#' 
#' ## Model fit with penalized cubic splines
#' pen <- TRUE
#' basis <- c(3,10,10)
#' formula <- ~s(x1, k=basis[1], bs = "cr", fx= !pen)+
#'   s(x2, k=basis[2], bs = "cr", fx= !pen)+
#'   s(x3, k=basis[3], bs = "cr", fx= !pen)
#' system.time(fit <- gamBiCopEst(dataset, family = fam, parFrhs = formula))
#' 
#' X <- as.data.frame(copula::rMvdc(n,covariates.distr))
#' names(X) <- c("x1", "x2", "x3")
#' sim <- gamBiCopSim(fit$res, X)
#' @export
gamBiCopSim <- 
  function(object, newdata = NULL, N = NULL, return.calib = FALSE, return.par = FALSE, return.tau = FALSE){
	
    if(!valid.gamBiCop(object)){
      stop("gamBiCopPred can only be used to predict from gamBiCop objects")
    }
    
  if(is.null(N)){
    if(!is.null(newdata)){
      N <- dim(newdata)[1]
    }else{
      N <- dim(object@model$data)[1]
    }
  }
	if(!is.na(as.integer(N))){
	  if( (as.integer(N) < 1) || (as.integer(N) != as.numeric(N))){
	    stop("N should be a positive integer!")
	  }else{
	    N <- as.integer(N)
	  }
	}else{
	  stop("N should be a positive integer!")
	}
	if(!(is.logical(return.calib) || (return.calib == 0) || (return.calib == 1))){
 		warning("Return.calib should takes 0/1 or FALSE/TRUE.")
 		return.calib = FALSE
 	}
 	if(!(is.logical(return.par) || (return.par == 0) || (return.par == 1))){
 		warning("Return.par should takes 0/1 or FALSE/TRUE.")
 		return.par = FALSE
 	}
 	if(!(is.logical(return.tau) || (return.tau == 0) || (return.tau == 1))){
 		warning("Return.tau should takes 0/1 or FALSE/TRUE.")
 		return.tau = FALSE
 	}
 	if(!is.na(as.integer(N))){
		if( (as.integer(N) < 1) || (as.integer(N) != as.numeric(N))){
			stop("N should be a positive integer!")
		}else{
			N <- as.integer(N)
		}
	}else{
		stop("N should be a positive integer!")
	}
  
  if(is.null(newdata)){
    newdata <- object@model$data
  }
  dd <- dim(newdata)[1]
  if((N < dd) || (N > dd)){
    newdata <- newdata[sample.int(dd,N,replace = TRUE),]
  }

	temp <- gamBiCopPred(object, newdata, target = c("par", "calib", "tau"))
  
  data = t(mapply(BiCopSim, 1, object@family, temp$par, object@par2)) 
  out <- list(data = data)
  	
  if(return.calib==TRUE){
		out$calib <-  temp$calib
	}
  	
	if(return.par==TRUE){
		out$par <-  temp$par
	}
	
	if(return.tau==TRUE){
		out$tau <-  temp$tau
	}
	
	return(out)
}
