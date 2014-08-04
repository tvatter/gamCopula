#' Predict from fitted \code{\link{GAMBiCop-class}} object
#'
#' @param object fitted \code{\link{GAMBiCop-class}} object.
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
#' @seealso \code{\link{GAMBiCop}} and \code{\link{GAMBiCopPred}}.
#' @export
#' @examples
#' ##    													
#' ##	1) Simulate and plot data	
#' ##	2) Estimate the model and display results
#' ##															
#' 
#' library(GAMVineCopula)
#' set.seed(1)
#' 
#' ##  Simulation parameters
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
#' ## 1) Simulate and plot data
#' ##
#' 
#' ## Integration grid
#' step <- 1e-3
#' ngrid <- 1/step
#' xx <- seq(0,1,length.out = ngrid)
#' Xx <- data.frame(cbind(xx,xx,xx))
#' names(Xx) <- c("x1", "x2", "x3")
#' true <- true.approx <- true.approx2 <- matrix(NA, ngrid, dim(Xx)[2])
#' 
#' ## Quadratic calibration
#' Tf <- 1
#' b.quad <- 8*Tf
#' t0.quad <- Tf/2
#' a.quad <- -(b.quad/3)*(Tf^2-3*Tf*t0.quad+3*t0.quad^2)
#' calib.quadratic <- function(t, t0, a, b){return(a + b*(t-t0)^2)}
#' true[,1] <- calib.quadratic(xx, t0.quad, a.quad, b.quad)
#' 
#' ## Sinusoidal calibration
#' b.sin <- 1
#' f.sin <- 1
#' t0.sin <- 0
#' a.sin <- b.sin*(1-2*f.sin*Tf*pi/(f.sin*Tf*pi+cos(2*f.sin*pi*(Tf-t0.sin))-cos(2*f.sin*pi*t0.sin)))
#' calib.sinusoidal <- function(t, f, t0, a, b){return((a+b)/2 + (b-a)*sin(2*pi*f*(t-t0))/2)}
#' true[,2] <- calib.sinusoidal(xx, f.sin, t0.sin, a.sin, b.sin)
#' 
#' ## Exponential calibration 
#' t0.exp  <- Tf/2
#' s.exp <- Tf/8
#' b.exp <- 2
#' a.exp <- (b.exp*s.exp*sqrt(2*pi)/Tf)*(pnorm(0,t0.exp,s.exp)-pnorm(Tf,t0.exp,s.exp))
#' calib.exponential <- function(t, a, b, t0, s){return(a+b*exp(-(t-t0)^2/(2*s^2)))}
#' true[,3] <- calib.exponential(xx, a.exp, b.exp, t0.exp, s.exp)
#' 
#' ## Base level
#' eta0 <- 1
#' 
#' ## Additive calibration function
#' calib.fct <- function(x1,x2,x3){
#'   return(eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#'          calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin)+
#'          calib.exponential(x3, a.exp, b.exp, t0.exp, s.exp))}
#' calib.fct12 <- function(x1,x2){eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#'                                  calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin)}
#' calib.fct13 <- function(x1,x2){eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#'                                  calib.exponential(x2, a.exp, b.exp, t0.exp, s.exp)}
#' calib.fct23 <- function(x1,x2){eta0+calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin) +
#'                                  calib.exponential(x1, a.exp, b.exp, t0.exp, s.exp)}
#' 
#' ## Display calibration surface
#' require(plot3D)
#' xx2 <- seq(0,Tf,length.out = 1e2)
#' Z12 <- outer(xx2,xx2,calib.fct12)
#' Z13 <- outer(xx2,xx2,calib.fct13)
#' Z23 <- outer(xx2,xx2,calib.fct23)
#' 
#' par(mfrow= c(1,3),  pty = "s", mar=c(1,1,4,1))
#' persp3D(xx2,xx2,Z12, theta = 60, phi = 30, expanded = 2, ticktype = "simple", nticks = 4, colkey = list(plot = FALSE), xlab = "X1", ylab = "X2", zlab = "", main = "Calib. fct(X1,X2)", cex.main = 1)
#' persp3D(xx2,xx2,Z13, theta = 60, phi = 30, expanded = 2, ticktype = "simple", nticks = 4, colkey = list(plot = FALSE), xlab = "X1", ylab = "X3", zlab = "", main = "Calib. fct(X1,X3)", cex.main = 1)
#' persp3D(xx2,xx2,Z23, theta = 60, phi = 30, expanded = 2, ticktype = "simple", nticks = 4, colkey = list(plot = FALSE), xlab = "X2", ylab = "X3", zlab = "", main = "Calib. fct(X2,X3)", cex.main = 1)
#' 
#' ## A single dataset
#' covariates.distr <- copula:::mvdc(copula:::normalCopula(rho, dim = 3), c("unif"), list(list(min = 0, max = Tf)), marginsIdentical = TRUE) 
#' X <- copula:::rMvdc(n,covariates.distr)
#' temp <- CondBiCopSim(fam, calib.fct, X, par2=par2, return.par=TRUE)
#' calib <- temp$eta
#' param <- temp$par
#' pseudo <- temp$data
#' dataset <- data.frame("u1" = pseudo[,1], "u2" = pseudo[,2], "x1" = X[,1], "x2" = X[,2], "x3" = X[,3])
#' 
#' ## Display the data
#' dev.off()
#' plot(pseudo[,1], pseudo[,2], xlab = "U1", ylab = "U2")
#' 
#' ## Display the copula parameter and calibration function
#' par(mfrow=c(2,3))
#' plot(X[,1],param, xlab = "X1", ylab = "Copula parameter")
#' plot(X[,2],param, xlab = "X2", ylab = "Copula parameter")
#' plot(X[,3],param, xlab = "X3", ylab = "Copula parameter")
#' plot(X[,1],calib, xlab = "X1", ylab = "Calibration function")
#' plot(X[,2],calib, xlab = "X2", ylab = "Calibration function")
#' plot(X[,3],calib, xlab = "X3", ylab = "Calibration function")
#' 
#' ##
#' ## 2) Estimate two models and display results
#' ##
#' 
#' ## Model fit with a basis size (arguably) too small and unpenalized cubic spines
#' pen <- FALSE
#' basis <- c(3,3,3)
#' formula <- ~s(x1, k=basis[1], bs = "cr", fx= !pen)+
#'         s(x2, k=basis[2], bs = "cr", fx= !pen)+
#'         s(x3, k=basis[3], bs = "cr", fx= !pen)
#' system.time(fit <- GAMBiCopEst(dataset, family = fam, parFrhs = formula))
#' 
#' ## Model fit with a better basis size and penalized cubic splines
#' pen <- TRUE
#' basis2 <- c(3,10,10)
#' formula <- ~s(x1, k=basis2[1], bs = "cr", fx= !pen)+
#'   s(x2, k=basis2[2], bs = "cr", fx= !pen)+
#'   s(x3, k=basis2[3], bs = "cr", fx= !pen)
#' system.time(fit2 <- GAMBiCopEst(dataset, family = fam, parFrhs = formula))
#' 
#' ## Extract the GAMBiCop object and use various methods
#' res <- fit$res
#' res2 <- fit2$res
#' show(res)
#' show(res2)
#' logLik(res)
#' logLik(res2)
#' AIC(res)
#' AIC(res2)
#' BIC(res)
#' BIC(res2)
#' fitted <- GAMBiCopPred(res, newdata = Xx, type="terms")$calib
#' fitted2 <- GAMBiCopPred(res2, newdata = Xx, type="terms")$calib
#' 
#' ## Spline approximation of each true smooth function for the two basis sizes
#' for(i in 1:length(basis)){
#'   temp <- rep(0,3)
#'   temp[i] <- 1
#'   temp <- data.frame("y" = calib.fct(xx*temp[1],xx*temp[2],xx*temp[3]), "x1" = xx, "x2" = 0, "x3" = 0)  
#'   form <- ~s(x1, k=basis[i], bs = "cr", fx= !pen)
#'   temp <- gam(update(form,y~.), data = temp)
#'   true.approx[,i] <- predict.gam(temp, type = "terms")
#'   
#'   temp <- rep(0,3)
#'   temp[i] <- 1
#'   temp <- data.frame("y" = calib.fct(xx*temp[1],xx*temp[2],xx*temp[3]), "x1" = xx, "x2" = 0, "x3" = 0)  
#'   form <- ~s(x1, k=basis2[i], bs = "cr", fx= !pen)
#'   temp <- gam(update(form,y~.), data = temp)
#'   true.approx2[,i] <- predict.gam(temp, type = "terms")
#' }
#' 
#' ## Compute approximation and fitted biases
#' bias.approx <- true.approx - true
#' bias.approx2 <- true.approx2 - true
#' bias.fitted2 <- fitted2 - true
#' bias.fitted <- fitted - true
#' 
#' ## Display results
#' yy <- range(true, true.approx, fitted, true.approx2, fitted2)
#' yy[1] <- yy[1]*1.5
#' yy.bias <- range(bias.approx, bias.fitted, bias.approx2, bias.fitted2)
#' par(mfrow = c(2,3))
#' plot(xx, true[,1], type = "l", ylim = yy, xlab = "Covariate 1", ylab = "Smooth 1")
#' lines(xx, fitted[,1], col = "red")
#' lines(xx, true.approx[,1], col = "red", lty = 2)
#' lines(xx, fitted2[,1], col = "green")
#' lines(xx, true.approx2[,1], col = "green", lty = 2)
#' legend("bottomleft", cex = 0.6, c("True", "Fitted", "Appox 1", "Fitted 2", "Approx 2"), col = c("black", "red", "red", "green", "green"), lty = c(1,1,2,1,2))
#' 
#' plot(xx, true[ ,2], type = "l", ylim = yy, xlab = "Covariate 1", ylab = "Smooth 1")
#' lines(xx, fitted[ ,2], col = "red")
#' lines(xx, true.approx[ ,2], col = "red", lty = 2)
#' lines(xx, fitted2[ ,2], col = "green")
#' lines(xx, true.approx2[ ,2], col = "green", lty = 2)
#' 
#' plot(xx, true[ ,3], type = "l", ylim = yy, xlab = "Covariate 1", ylab = "Smooth 1")
#' lines(xx, fitted[ ,3], col = "red")
#' lines(xx, true.approx[ ,3], col = "red", lty = 2)
#' lines(xx, fitted2[ ,3], col = "green")
#' lines(xx, true.approx2[ ,3], col = "green", lty = 2)
#' 
#' plot(xx, bias.fitted[,1], col = "red", type = "l", ylim = yy.bias, xlab = "Covariate 1", ylab = "Bias 1")
#' lines(xx, bias.approx[,1], col = "red", lty = 2)
#' lines(xx, bias.fitted2[,1], col = "green", lty = 1)
#' lines(xx, bias.approx2[,1], col = "green", lty = 2)
#' legend("bottomleft", cex = 0.6, c("Fitted", "Appox 1", "Fitted 2", "Approx 2"), col = c("red", "red", "green", "green"), lty = c(1,2,1,2))
#' 
#' plot(xx, bias.fitted[,2], col = "red", type = "l", ylim = yy.bias, xlab = "Covariate 1", ylab = "Bias 1")
#' lines(xx, bias.approx[,2], col = "red", lty = 2)
#' lines(xx, bias.fitted2[,2], col = "green", lty = 1)
#' lines(xx, bias.approx2[,2], col = "green", lty = 2)
#' 
#' plot(xx, bias.fitted[,3], col = "red", type = "l", ylim = yy.bias, xlab = "Covariate 1", ylab = "Bias 1")
#' lines(xx, bias.approx[,3], col = "red", lty = 2)
#' lines(xx, bias.fitted2[,3], col = "green", lty = 1)
#' lines(xx, bias.approx2[,3], col = "green", lty = 2)
GAMBiCopPred <- 
  function(object, newdata = NULL, target = "calib", alpha = 0, type = "link"){
	
	if(!validGAMBiCop(object)){
		stop("GAMBiCopPred can only be used to predict from GAMBiCop objects")
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