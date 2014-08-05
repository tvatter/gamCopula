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
