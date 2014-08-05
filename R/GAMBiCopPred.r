#' Predict from fitted \code{\link{gamBiCop-class}} object
#'
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @param newdata (Same as in \code{\link{predict.gam}} from the \code{\link[mgcv:mgcv-package]{mgcv}} package) A matrix or data frame containing the values of the model covariates at which predictions are required. 
#' If this is not provided then predictions corresponding to the original data are returned. 
#' If newdata is provided then it should contain all the variables needed for prediction: 
#' a warning is generated if not.
#' @param target Either \code{"calib"}, \code{"par"} or \code{"tau"} or a combination of those. 
#' \code{"calib"} (default) corresponds to the calibration function, 
#' \code{"par"} to the copula parameter and \code{"tau"} to Kendall's tau.
#' @param alpha In (0,1) to return the corresponding confidence interval.
#' @param type (Similar as in \code{\link{predict.gam}} from the \code{\link[mgcv:mgcv-package]{mgcv}} package, 
#' only active for \code{type = "calib"}). When this has the value \code{"link"} (default), the calibration function is returned. 
#'  When \code{type = "terms"} each component of the linear predictor is returned seperately (possibly with standard errors): 
#'	this includes parametric model components, followed by each smooth component, but excludes any offset and any intercept. 
#'	When \code{type = "lpmatrix"} then a matrix is returned which yields the values of the linear predictor (minus any offset) 
#'	when postmultiplied by the parameter vector (in this case alpha is ignored).
#' @return If \code{target = "calib"}, then a list with 1 item \code{calib}. 
#' If \code{target = "par"}, \code{target = "tau"} or \code{target = c("par", "tau")}, 
#' then a list with 2, 2 or 3 items, namely \code{calib} and \code{par},  \code{tau} and \code{par}, 
#' or  \code{calib}, \code{tau} and \code{par}.
#' 
#' If \code{alpha} is in (0,1), then a additional items of the list are \code{calib.CI}
#'  as well as e.g. \code{par.CI} and/or \code{tau.CI} depending on the value of \code{target}.
#'  
#' Otherwhise, if \code{type = "lpmatrix"} (only active for \code{type = "calib"}),
#' then a matrix is returned which will give a vector of linear predictor values (minus any offest) at the supplied covariate
#' values, when applied to the model coefficient vector (similar as \code{\link{predict.gam}} from the \code{\link[mgcv:mgcv-package]{mgcv}}). 
#' @seealso \code{\link{gamBiCop}} and \code{\link{gamBiCopPred}}.
#' @export
gamBiCopPred <- 
  function(object, newdata = NULL, target = "calib", alpha = 0, type = "link"){
	
	if(!validgamBiCop(object)){
		stop("gamBiCopPred can only be used to predict from gamBiCop objects")
	}
    if(!is.character(target) || !is.null(dim(target))){
    	targerr <- TRUE
    }else if(length(target) == 1 && !is.element(target, c("calib", "par", "tau"))) {
        targerr <- TRUE
    }else if(length(target) > 1 && !all(is.element(target, c("calib", "par", "tau")))){
    	targerr <- TRUE
    }else if(length(target) > 3){
    	targerr <- TRUE
    }else{
    	targerr <- FALSE
    }
    if(targerr){
    	warning("Unknown target, reset to calib.")
        target <- "calib"
    }
    
    if (target == "calib" && type != "link" && type != "terms" && type != "lpmatrix") {
        warning("Unknown type, reset to terms.")
        type <- "terms"
    }
	
	options(warn = -1)	
	if(!is.na(as.double(alpha))){
		if( (as.double(alpha) < 0) || (as.double(alpha) > 1)){
			stop("Alpha should be a real number in [0,1]!")
		}else{
			alpha <- as.double(alpha)
		}
	}else{
		stop("Alpha should be a real number in [0,1]!")
	}
	options(warn = 0)

	mm <- object@model
	
	rotated <- family <- object@family
	if(is.element(rotated,c(13,23,33))){
		rotated <- 3
	}else if(is.element(family,c(14,24,34))){
		rotated <- 4
	}
	
	out <- list()
	sel <- length(target) == 1 && (target == "calib")
	if(sel){
    if(is.null(newdata)){
      out$calib <- predict(mm, type = type)      
    }else{
      out$calib <- predict(mm, newdata, type = type)
    }
	}else{
	  if(is.null(newdata)){
	    out$calib <- predict(mm)      
	  }else{
	    out$calib <- predict(mm, newdata)
	  }
	}

	if(!(type == "lpmatrix") && (alpha != 0) && (alpha != 1)){
		Xp <- predict(mm,newdata, type = "lpmatrix")
		b <- coef(mm)
		Vp <- vcov(mm)
		br <- MASS::mvrnorm(1e4,b,Vp)
		calib <- br %*% t(Xp)
		out$calib.CI <- t(apply(calib,2,function(x) quantile(x, c((1-alpha)/2,1-(1-alpha)/2))))
	}else{
		calib = NULL	
	}
	
	if(any(is.element(target, "par"))){
		if(object@tau){
			tmp <- tanh(out$calib/2)
			if(rotated %in% c(3,4)){
				out$par <- sapply((1+tmp)/2, function(x) BiCopTau2Par(rotated,x))
				if(family %in% c(23,24,33,34)){
					out$par <- -out$par
				}
				if(!is.null(calib)){
					tmp <- sapply((1+tanh(calib/2))/2, function(x) BiCopTau2Par(rotated,x))
					out$par.CI <- t(apply(tmp,2,function(x) quantile(x, c((1-alpha)/2,1-(1-alpha)/2))))
					if(family %in% c(23,24,33,34)){
						out$par.CI <- -out$par.CI
					}
				}
			}else{
				out$par <- sapply(tmp, function(x) BiCopTau2Par(1,x))
				if(!is.null(calib)){
					tmp <- sapply(tanh(calib/2), function(x) BiCopTau2Par(1,x))
					out$par.CI <- t(apply(tmp,2,function(x) quantile(x, c((1-alpha)/2,1-(1-alpha)/2))))
				}
			}	
		}else{
			out$par <- sapply(out$calib, function(x) BiCopEta2Par(rotated,x))
			if(family %in% c(23,24,33,34)){
					out$par <- -out$par
			}
			if(!is.null(calib)){
				tmp <- sapply(calib, function(x) BiCopEta2Par(rotated,x))
				out$par.CI <- t(apply(tmp,2,function(x) quantile(x, c((1-alpha)/2,1-(1-alpha)/2))))
				if(family %in% c(23,24,33,34)){
					out$par.CI <- -out$par.CI
				}
			}
		}
	}
	
	if(any(is.element(target, "tau"))){
		if(object@tau){
			out$tau <- tanh(out$calib/2)
			if(rotated %in% c(3,4)){
				out$tau <- (1+out$tau)/2
				if(family %in% c(23,24,33,34)){
					out$tau <- -out$tau
				}
			}
			if(!is.null(calib)){
				out$tau.CI <- t(apply(tanh(calib/2),2,function(x) quantile(x, c((1-alpha)/2,1-(1-alpha)/2))))
				if(rotated %in% c(3,4)){
					out$tau.CI <- (1+out$tau.CI)/2
					if(family %in% c(23,24,33,34)){
						out$tau.CI <- -out$tau.CI
					}
				}
			}
		}else{
			tmp <- sapply(out$calib, function(x) BiCopEta2Par(rotated,x))
			out$tau <- sapply(tmp, function(x) BiCopPar2Tau(rotated,x))
			if(family %in% c(23,24,33,34)){
					out$tau <- -out$tau
			}
			if(!is.null(calib)){
					tmp <- sapply(calib, function(x) BiCopEta2Par(rotated,x))
					tmp <- sapply(tmp, function(x) BiCopPar2Tau(rotated,x))
					out$tau.CI <- t(apply(tmp,2,function(x) quantile(x, c((1-alpha)/2,1-(1-alpha)/2))))
					if(family %in% c(23,24,33,34)){
						out$tau.CI <- -out$tau.CI
					}
			}
		}
	}
	
	return(out)
}
