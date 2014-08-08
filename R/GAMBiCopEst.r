#' Estimate a Generalized Additive model (gam) for the copula parameter or Kendall's tau
#'
#' @param family A copula family: \code{1} Gaussian, \code{2} Student t, 
#' \code{3} Clayton, \code{4} Gumbel, \code{13} Survival Clayton, \code{14} Survival Gumbel, 
#' \code{23} Rotated (90 degrees) Clayton, \code{24} Rotated (90 degrees) Gumbel, 
#' \code{33} Rotated (270 degrees) Clayton and \code{34} Rotated (270 degrees) Gumbel.
#' @param parFrhs A gam formula (see \code{\link{gam}}, \code{\link{formula.gam}} and \code{\link{gam.models}} from the \code{\link[mgcv:mgcv-package]{mgcv}} package).
#' @param dataset A matrix or data frame containing the model responses (in [0,1]x[0,1]) and covariates required by the formula.
#' @param tau \code{FALSE} (default) for a calibration fonction specified for the Copula parameter 
#' or \code{TRUE} for a calibration function specified for Kendall's tau.
#' @param method \code{"FS"} for Fisher-scoring and \code{"NR"} for Newton-Raphson.
#' @param tol.rel Relative tolerance for \code{"FS"}/\code{"NR"} algorithm.
#' @param n.iters Maximal number of iterations for \code{"FS"}/\code{"NR"} algorithm.
#' @param verbose \code{TRUE} prints informations during the estimation.
#' @return \code{gamBiCopEst} returns a list consisting of
#' \item{res}{S4 \code{\link{gamBiCop-class}} object.}
#' \item{method}{\code{"FS"} for Fisher-scoring and \code{"NR"} for Newton-Raphson.}
#' \item{tol.rel}{relative tolerance for \code{"FS"}/\code{"NR"} algorithm.}
#' \item{n.iters}{maximal number of iterations for \code{"FS"}/\code{"NR"} algorithm.}
#' \item{trace}{the estimation procedure's trace.}
#' @seealso \code{\link{gamBiCop}} and \code{\link{gamBiCopEst}}.
#' @examples
#' ##  														
#' ##	1) Simulate and plot data	
#' ##  a) Kendall's tau as a sum of three smooth components
#' ##  b) Distribution of the covariates as a Gaussian copula 
#' ##	2) Estimate the model with two different basis size and display results 							
#' 
#' require(gamVineCopula)
#' 
#' set.seed(0)
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
#' xx2 <- seq(0,Tf,length.out = 1e2)
#' Z12 <- outer(xx2,xx2,calib.fct12)
#' Z13 <- outer(xx2,xx2,calib.fct13)
#' Z23 <- outer(xx2,xx2,calib.fct23)
#' 
#' par(mfrow= c(1,3),  pty = "s", mar=c(1,1,4,1))
#' plot3D::persp3D(xx2,xx2,Z12, theta = 60, phi = 30, expanded = 2, ticktype = "simple", 
#'    nticks = 4, colkey = list(plot = FALSE), xlab = "X1", ylab = "X2", zlab = "", 
#'    main = "Calib. fct(X1,X2)", cex.main = 1)
#' plot3D::persp3D(xx2,xx2,Z13, theta = 60, phi = 30, expanded = 2, ticktype = "simple", 
#'    nticks = 4, colkey = list(plot = FALSE), xlab = "X1", ylab = "X3", zlab = "", 
#'    main = "Calib. fct(X1,X3)", cex.main = 1)
#' plot3D::persp3D(xx2,xx2,Z23, theta = 60, phi = 30, expanded = 2, ticktype = "simple", 
#'    nticks = 4, colkey = list(plot = FALSE), xlab = "X2", ylab = "X3", zlab = "", 
#'    main = "Calib. fct(X2,X3)", cex.main = 1)
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
#' system.time(fit <- gamBiCopEst(dataset, family = fam, parFrhs = formula))
#' 
#' ## Model fit with a better basis size and penalized cubic splines
#' pen <- TRUE
#' basis2 <- c(3,10,10)
#' formula <- ~s(x1, k=basis2[1], bs = "cr", fx= !pen)+
#'   s(x2, k=basis2[2], bs = "cr", fx= !pen)+
#'   s(x3, k=basis2[3], bs = "cr", fx= !pen)
#' system.time(fit2 <- gamBiCopEst(dataset, family = fam, parFrhs = formula))
#' 
#' ## Extract the gamBiCop object and use various methods
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
#' EDF.gamBiCop(res)
#' EDF.gamBiCop(res2)
#' fitted <- gamBiCopPred(res, newdata = Xx, type="terms")$calib
#' fitted2 <- gamBiCopPred(res2, newdata = Xx, type="terms")$calib
#' 
#' ## Spline approximation of each true smooth function for the two basis sizes
#' for(i in 1:length(basis)){
#'   temp <- rep(0,3)
#'   temp[i] <- 1
#'   pen <- TRUE
#'   temp <- data.frame("y" = calib.fct(xx*temp[1],xx*temp[2],xx*temp[3]), 
#'      "x1" = xx, "x2" = 0, "x3" = 0)  
#'   form <- ~s(x1, k=basis[i], bs = "cr", fx= !pen)
#'   temp <- mgcv::gam(update(form,y~.), data = temp)
#'   true.approx[,i] <- mgcv::predict.gam(temp, type = "terms")
#'   
#'   temp <- rep(0,3)
#'   temp[i] <- 1
#'   pen <- FALSE
#'   temp <- data.frame("y" = calib.fct(xx*temp[1],xx*temp[2],xx*temp[3]), 
#'      "x1" = xx, "x2" = 0, "x3" = 0)  
#'   form <- ~s(x1, k=basis2[i], bs = "cr", fx= !pen)
#'   temp <- mgcv::gam(update(form,y~.), data = temp)
#'   true.approx2[,i] <- mgcv::predict.gam(temp, type = "terms")
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
#' legend("bottomleft", cex = 0.6, c("True", "Fitted", "Appox 1", "Fitted 2", "Approx 2"), 
#'    col = c("black", "red", "red", "green", "green"), lty = c(1,1,2,1,2))
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
#' plot(xx, bias.fitted[,1], col = "red", type = "l", ylim = yy.bias, 
#'    xlab = "Covariate 1", ylab = "Bias 1")
#' lines(xx, bias.approx[,1], col = "red", lty = 2)
#' lines(xx, bias.fitted2[,1], col = "green", lty = 1)
#' lines(xx, bias.approx2[,1], col = "green", lty = 2)
#' legend("bottomleft", cex = 0.6, c("Fitted", "Appox 1", "Fitted 2", "Approx 2"), 
#'    col = c("red", "red", "green", "green"), lty = c(1,2,1,2))
#' 
#' plot(xx, bias.fitted[,2], col = "red", type = "l", ylim = yy.bias, 
#'    xlab = "Covariate 1", ylab = "Bias 1")
#' lines(xx, bias.approx[,2], col = "red", lty = 2)
#' lines(xx, bias.fitted2[,2], col = "green", lty = 1)
#' lines(xx, bias.approx2[,2], col = "green", lty = 2)
#' 
#' plot(xx, bias.fitted[,3], col = "red", type = "l", ylim = yy.bias, 
#'    xlab = "Covariate 1", ylab = "Bias 1")
#' lines(xx, bias.approx[,3], col = "red", lty = 2)
#' lines(xx, bias.fitted2[,3], col = "green", lty = 1)
#' lines(xx, bias.approx2[,3], col = "green", lty = 2)
#' @export
gamBiCopEst <- 
function(family = 1, parFrhs = ~1, dataset, tau = FALSE, method = "FS", tol.rel = 1e-3, n.iters = 10, verbose = FALSE){
	
	if(!is.matrix(dataset) && !is.data.frame(dataset)){
		stop("Dataset has to be either a matrix or a data frame")
	}else{
		n <- dim(dataset)[1]
		m <- dim(dataset)[2]
		u1 <- dataset[,1]
		u2 <- dataset[,2]
	}

  if(is.null(u1)==TRUE || is.null(u2)==TRUE) stop("u1 and/or u2 are not set or have length zero.")
  if(length(u1) != length(u2)) stop("Lengths of 'u1' and 'u2' do not match.")
  if(length(u1)<2) stop("Number of observations has to be at least 2.")
  if(any(u1>1) || any(u1<0)) stop("Data has be in the interval [0,1].")
  if(any(u2>1) || any(u2<0)) stop("Data has be in the interval [0,1].")
  

	temp <- VineCopula:::fasttau(u1,u2)
  
  	if(!is.element(family, c(1,2,3,4,13,14,23,24,33,34))){
		stop("Copula family not yet implemented!")
	}else if(is.element(family, c(3,4,13,14)) && (temp < 0)){
		stop("This copula family cannot be used for negatively dependent data.")
	}else if(is.element(family, c(23,24,33,34)) && (temp > 0)){
		stop("This copula family cannot be used for positively dependent data.")
	}else if(is.element(family, c(13,14))){
		dataset[,c(1,2)] <- rotate.data(dataset[,c(1,2)], 180)
		rotated <- family
		if(family == 13){family <- 3}else{family <- 4}
	}else if(is.element(family, c(23,24))){
		dataset[,c(1,2)] <- rotate.data(dataset[,c(1,2)], 90)
		rotated <- family
		if(family == 23){family <- 3}else{family <- 4}
	}else if(is.element(family, c(33,34))){
		dataset[,c(1,2)] <- rotate.data(dataset[,c(1,2)], 270)
		rotated <- family
		if(family == 33){family <- 3}else{family <- 4}
	}else{
		rotated <- family
	}

	rotated <- family
	if(is.element(family,c(13,23,33))){
		family <- 3
	}else if(is.element(family,c(14,24,34))){
		family <- 4
	}

  options(warn = -1)
	if(!is.na(as.integer(n.iters))){
		if( (as.integer(n.iters) < 1) || (as.integer(n.iters) != as.numeric(n.iters))){
			stop("N.iters should be a positive integer!")
		}else{
			n.iters <- as.integer(n.iters)
		}
	}else{
		stop("N.iters should be a positive integer!")
	}
	
	if(!(is.logical(tau)||(tau == 0)||(tau == 1))){
 		stop("Tau should takes 0/1 or FALSE/TRUE to model the copula parameter/Kendall's tau.")
 	} 

	if(!is.na(as.numeric(tol.rel))){
		if((as.numeric(tol.rel) < 0) || (as.numeric(tol.rel) > 1)){
			stop("Tol.rel should be a real number in [0,1]!")
		}else{
			tol.rel <- as.numeric(tol.rel)
		}
	}else{
		stop("Tol.rel should be a real number in [0,1]!")
	}

	if(!is.element(method, c("FS","NR"))){
		stop("Method should be a string, either FS (Fisher-scoring) or NR (Newton-Raphson, unstable)!")
	}
	
	if(!(is.logical(verbose)||(verbose == 0)||(verbose == 1))){
 		stop("Verbose should takes 0/1 or FALSE/TRUE.")
 	}else{
			verbose <- as.logical(verbose)
	}
	options(warn = 0)

    init <- BiCopEst(u1,u2,family=family,method="mle")
	data <- cbind(u1,u2)

	new.pars <- list()
	new.pars$par <- rep(init$par,n)
	data <- cbind(data, new.pars$par)

	if(family == 2){
		new.pars$par2 <- rep(init$par2,n)
		data <- cbind(data, new.pars$par2)
	}

	if(tau == FALSE){
		new.pars$partrans <- sapply(new.pars$par, function(x) BiCopPar2Eta(family, x))
	}else{
		if(family == 2){
			new.pars$tau <- sapply(new.pars$par, function(x) BiCopPar2Tau(family, x, init$par2))
		}else{
			new.pars$tau <- sapply(new.pars$par, function(x) BiCopPar2Tau(family, x))
		}
		new.pars$partrans <- 2*atanh(new.pars$tau)
	}

	old.pars <- new.pars

	if(verbose == 1){
			t <- Sys.time()
			print(paste("gam iteration", 1))				
	}

	temp <- derivatives.par(data, new.pars, family, method, tau)
	temp <- as.data.frame(wz.update(temp, new.pars, family, method, tau))	
	temp <- cbind(temp, dataset)
   	par.formula <- update(parFrhs,z~.)

    w <- NULL
   	res <- tryCatch({mm <- gam(par.formula, data=temp,weights = w)},
      					error = function(err){
      						print(paste("A problem occured at the first iteration of the ", 
      						method, "algorithm. The ERROR comming from mgcv's gam function is:"))
      						stop(err)	
      						print("......The results should not be trusted!")
      						return(NULL)
      					})
    if(verbose == 1){
		print(Sys.time()-t)
	}
    stopifnot(!is.null(res))
    
   	temp <- pars.update(mm, family, temp, tau)

	new.pars$par <- data[,3] <- temp$par
	new.pars$partrans <- temp$partrans
	if(tau){
		new.pars$tau <- temp$tau
	}
		
	if(family == 2){
		if(verbose == 1){
			print(paste("DF iteration", 1))	
			t <- Sys.time()		
		}		
			
		LL <- function(nu){
			if(2+1e-8+exp(nu) == Inf){
				nu <- log(30)
			}
			-sum(log(apply(data, 1, function(x) BiCopPDF(x[1], x[2],family=2,x[3],2+1e-8+exp(nu)))), na.rm = TRUE)
		}

		nu <- optimize(LL, c(log(2), log(30)))
		data[,4] <- new.pars$par2 <- rep(2+1e-8+exp(nu$minimum),n)
			
		if(verbose == 1){
			print(data[1,4])
		 	print(Sys.time()-t)	
		}
	}
	
	trace <- numeric(n.iters)
	tt <- trace.update(old.pars$partrans, new.pars$partrans)
	trace[1] <- tt$trace

	eps <- tt$eps
	k <- 1
	while((k < n.iters) & (eps > tol.rel)){
		
		k <- k + 1
		old.pars <- new.pars

		if(verbose == 1){
			t <- Sys.time()
			print(paste("gam iteration", k))				
		}

		temp <- derivatives.par(data, new.pars, family, method, tau)
		temp <- as.data.frame(wz.update(temp, new.pars, family, method, tau))
		temp <- cbind(temp, dataset)		

      	res <- tryCatch({mm <- gam(par.formula, data=temp,weights = w, control = gam.control(keepData=TRUE))},
      					error = function(err){
      						print(paste("A problem occured at the ", k,"th iteration of the ", 
      						method, "algorithm. The ERROR comming from mgcv's gam function is:"))
      						print(err)	
      						print("......The results should not be trusted!")
      						return(NULL)
      					})
     	if(is.null(res)){
      		break
      	}
      	if(verbose == 1){
				print(Sys.time()-t)
		}
    
    	temp2 <- pars.update(mm, family, temp, tau)
		new.pars$par <- data[,3] <- temp2$par
		new.pars$partrans <- temp2$partrans
		if(tau){
			new.pars$tau <- temp2$tau
		}
		
		if((family == 2)){# && (k %% 2 == 0)){
			if(verbose == 1){
				print(paste("DF iteration", k))	
				t <- Sys.time()	
			}		
				
			LL <- function(nu){
				if(2+1e-8+exp(nu) == Inf){
					nu <- log(30)
				}
				-sum(log(apply(data, 1, function(x) BiCopPDF(x[1], x[2],family=2,x[3],2+1e-8+exp(nu)))), na.rm = TRUE)
			}
						
			nu <- optimize(LL, c(log(2), log(30)))
			data[,4] <- new.pars$par2 <- rep(2+1e-8+exp(nu$minimum),n)
			
			 if(verbose == 1){
				print(data[1,4])
			 	print(Sys.time()-t)	
			}
		}
		
		tt <- trace.update(old.pars$partrans, new.pars$partrans)
		trace[k] <- abs(tt$trace)
		eps <- tt$eps
	}

	if(family == 2){
		res <- gamBiCop(rotated, mm, data[1,4], tau)
	}else{
		res <- gamBiCop(rotated, mm, 0, tau)
	}
	out <- list(res = res, method = method, tol.rel = tol.rel, n.iters = n.iters, trace = trace[1:k])
	return(out)
}
