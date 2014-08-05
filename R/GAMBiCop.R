#############################
##  gam bivariate copulas  ##
#############################
validgamBiCop = function(object) {	
 d <- length(attributes(object))
 if(d < 2){
  	return("A gamBiCop contains at least a copula family and a mgcv model.")
 }else if(d >= 2){
  	if(!is(object@model, "gam")){
  		return("Invalid mgcv model.")	
  	}else if(!(object@family %in% c(1,2,3,4,13,14,23,24,33,34))){
		return("Copula family not yet implemented.")
	}
 } 
 if(!(is.logical(object@tau) || (object@tau == 0) || (object@tau == 1))){
 	return("Tau should takes 0/1 or FALSE/TRUE to model the copula parameter/tau's tau.")
 } 
 if((object@family == 2) && (is.null(object@par2) || is.na(as.numeric(object@par2)) || as.numeric(object@par2) <= 2)){
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
setClass("gamBiCop",
         slots = c(family="integer", model = "gam", par2 = "numeric", tau = "logical"),
         validity = validgamBiCop
)

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
gamBiCop <- function (family, model, par2 = 0, tau = FALSE) {
 	if(family != 2){
 		par2 = 0
 	}
 	if(as.integer(family) != family){
 		return("Family should be an integer.")
 	}
 	if(!(is.logical(tau) || (tau == 0) || (tau == 1))){
 		return("Tau should takes 0/1 or FALSE/TRUE to model the copula parameter/Kendall's tau.")
 	} 
  	new("gamBiCop", family = as.integer(family), model = model, par2 = par2, tau = as.logical(tau))
}

show.gamBiCop <- function(object) {
  cat("Family: ", object@family, "\n")
  if(object@tau == TRUE){
  	cat("Model: ")
  	cat("tau(z) = (exp(z)-1)/(exp(z)+1) where \n")
  }else{
  	 cat("Model for the Copula parameter:\n")
  	 if(object@family %in% c(1,2)){
  		cat("par(z) = (exp(z)-1)/(exp(z)+1) where \n")
  	 }else if(object@family %in% c(3,13)){
  	   	cat("par(z) = exp(z) where \n")	 	
  	 }else if(object@family %in% c(4,14)){
  	   	cat("par(z) = 1+exp(z) where \n")	  	 	
  	 }else if(object@family %in% c(23,33)){
  	   	cat("par(z) = -exp(z) where \n")	 	
  	 }else if(object@family %in% c(24,34)){
  	   	cat("par(z) = -1-exp(z) where \n")	
  	 }
  }
  show(object@model$formula)
}
setMethod("show", signature("gamBiCop"), show.gamBiCop)

#' Extract the Number of Obserations from a fitted \code{\link{gamBiCop-class}} object
#' 
#' Extract the number of 'observations' from a model fit. 
#' This is principally intended to be used in computing BIC (see \code{\link{AIC}}).
#'
#' @S3method nobs gamBiCop
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @return A single number, normally an integer.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @name nobs.gamBiCop
#' @rdname nobs-methods
#' @export
nobs.gamBiCop <- function(object, ...) {
  n <- dim(object@model$data)[1]
  return(n)
}
#' @docType methods
#' @rdname nobs-methods
setMethod("nobs", signature("gamBiCop"), nobs.gamBiCop)


#' Log-likelihood for a fitted \code{gamBiCop} object
#' 
#' Function to extract the log-likelihood for a fitted \code{gamBiCop-class}
#' object (note that the models are usually fitted by penalized likelihood maximization). 
#' Used by \code{\link{AIC}} and \code{\link{BIC}}.
#'
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @return Standard \code{logLik} object: see \code{\link{logLik}}.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @rdname logLik-methods
#' @export
logLik.gamBiCop <- function(object, ...) {
  family <- object@family  
  par <- gamBiCopPred(object, target = "par")$par
  data <- cbind(object@model$data[,c(3,4)],par)
  
  if(family == 2){
    par2 <- rep(object@par2, length(par))
    data <- cbind(data, par2)
    Li<- apply(data, 1, function(x) BiCopPDF(x[1], x[2],family=2,x[3],x[4]))
  }else{
    Li<- apply(data, 1, function(x) BiCopPDF(x[1], x[2],family=family,x[3]))
  }
  val <- sum(log(Li))
  df <- sum(object@model$edf)
  
  if(family == 2){
    df <- df+1    
  }
  
  attr(val, "df") <- df
  attr(val, "nobs") <- dim(data)[1]
  class(val) <- "logLik"
  return(val)
}
#' @docType methods
#' @rdname logLik-methods
setMethod("logLik", signature("gamBiCop"), logLik.gamBiCop)

#' Akaike's 'An Information Criterion' for a fitted \code{\link{gamBiCop-class}}
#' 
#' Function calculating Akaike's 'An Information Criterion' (AIC) for a fitted \code{\link{gamBiCop-class}}
#' object (note that the models are usually fitted by penalized likelihood maximization). 
#'
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @param k numeric, the penalty per parameter to be used; the default \code{k = 2} is the classical AIC.
#' @return A numeric value with the corresponding AIC.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @rdname AIC-methods
#' @export
AIC.gamBiCop <- function(object, ..., k = 2) {
  l <- logLik.gamBiCop(object)
  d <- attributes(l)$df
  return(k*d-2*l[1])
}
#' @docType methods
#' @rdname AIC-methods
setMethod("AIC", signature("gamBiCop"), AIC.gamBiCop)

#' Schwarz's Bayesian Information Criterion for a fitted \code{\link{gamBiCop-class}}
#' 
#' Function calculating the Schwarz's Bayesian Information Criterion (BIC) 
#' for a fitted \code{\link{gamBiCop-class}} object 
#' (note that the models are usually fitted by penalized likelihood maximization). 
#'
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @return A numeric value with the corresponding BIC.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @rdname BIC-methods
BIC.gamBiCop <- function(object, ...){
  return(AIC.gamBiCop(object, ..., k = log(nobs(object))))
}
#' @docType methods
#' @rdname BIC-methods
#' @export
setMethod("BIC", signature("gamBiCop"), BIC.gamBiCop)


#' \code{\link{gamBiCop-class}} formula
#' 
#' Description of the \code{\link{gam}} formula for  fitted \code{\link{gamBiCop-class}} object.
#' This function is a wrapper to \code{\link{formula.gam}}
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param x fitted \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @seealso \code{\link{formula.gam}} function 
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#' @docType methods
#' @rdname formula-methods
formula.gamBiCop <- function(x, ...){
  return(x@model$formula)
}
#' @docType methods
#' @rdname formula-methods
#' @export
setMethod("formula", signature("gamBiCop"), formula.gamBiCop)

#' Equivalent Degrees of Freedom for a fitted \code{\link{gamBiCop-class}}
#' 
#' Function calculating the Equivalent Degrees of Freedom (EDF) 
#' for a fitted \code{\link{gamBiCop-class}} object. 
#' It basically sums the edf of the \code{\link{gamObject}} 
#' for each smooth component.
#'
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @return Estimated degrees of freedom for each smooth component.
#' @docType methods
#' @rdname EDF-methods
#' @export
EDF.gamBiCop <- function(object){
  
  edf <- object@model$edf[-1]
  
  param.terms <- object@model$pterms
  ll.param <- dim(attributes(param.terms)$factors)[2]
  if(is.null(ll.param)){
    ll.param <- 0
  }
  
  smooth.terms <- object@model$smooth
  ll.smooth <- length(smooth.terms)
  bs.dim <- unlist(lapply(smooth.terms, function(x) x$bs.dim))-1
  
  out <- rep(NA, ll.param+ll.smooth+1)
  out[1:(ll.param+1)] <- 1
  sel <- c(ll.param,ll.param+cumsum(bs.dim))
  if(length(sel) > 1){for(i in 1:(length(sel)-1)){
    out[ll.param+1+i] <- sum(edf[(sel[i]+1):sel[i+1]])
  }}

  return(out)
}
  