#' Simulate from fitted \code{\link{GAMBiCop-class}} object
#'
#' @param object fitted \code{\link{GAMBiCop-class}} object.
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
#' @export
GAMBiCopSim <- 
  function(object, newdata = NULL, N = NULL, return.calib = FALSE, return.par = FALSE, return.tau = FALSE){
	
	if(!validGAMBiCop(object)){
		stop("GAMBiCopPred can only be used to predict from GAMBiCop objects")
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
  
	temp <- GAMBiCopPred(object@model, newdata, target = c("par", "eta", "tau"))
  
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
