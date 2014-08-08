#' Simulation from a \code{\link{gamVineCop-class}} object
#'
#' @param N number of d-dimensional observations to simulate.
#' @param GVC \code{\link{gamVineCop-class}} object.
#' @param U (similar as \code{\link{RVineSim}} from the from the \code{\link{VineCopula}} package) 
#' If not NULL, an (N,d)-matrix of U[0,1] random variates to be transformed to the copula sample.
#' @return A Nxd matrix of data simulated from the given \code{\link{gamVineCop-class}} object.
#' @examples
#' #########
#' ## Comparison between RVine package's simulation (in C) with the same algorithm in R
#' #########
#' 
#' require(VineCopula)
#' 
#' # define 5-dimensional R-vine tree structure matrix
#' Matrix = c(5,2,3,1,4,0,2,3,4,1,0,0,3,4,1,0,0,0,4,1,0,0,0,0,1)
#' Matrix = matrix(Matrix,5,5)
#' 
#' # define R-vine pair-copula family matrix
#' family = c(0,1,3,4,4,0,0,3,4,1,0,0,0,4,1,0,0,0,0,3,0,0,0,0,0)
#' family = matrix(family,5,5)
#' 
#' # define R-vine pair-copula parameter matrix
#' par = c(0,0.2,0.9,1.5,3.9,0,0,1.1,1.6,0.9,0,0,0,1.9,0.5,
#'         0,0,0,0,4.8,0,0,0,0,0)
#' par = matrix(par,5,5)
#' 
#' # define second R-vine pair-copula parameter matrix
#' par2 = matrix(0,5,5)
#' 
#' # define RVineMatrix object
#' RVM = RVineMatrix(Matrix=Matrix,family=family,par=par,par2=par2,
#'                   names=c("V1","V2","V3","V4","V5"))
#' 
#' # convert as a gamVineCop object
#' GVC <- RVM2GVC(RVM)
#' 
#' # vector of sample sizes
#' Nn <- rep(c(1,2,5), 3)*10^c(rep(1,3),rep(2,3), rep(3,3))
#' 
#' # to store the times
#' t1 <- t2 <- rep(0, length(Nn))
#' 
#' # simulate a sample for each size / method, 
#' # store the times
#' for(i in 1:length(Nn)){
#'   N <- Nn[i]
#'   U <- matrix(runif(N*5), ncol = 5)
#'   t1[i] <- system.time(Sim1 <- RVineSim(N, RVM, U))[3]
#'   t2[i] <- system.time(Sim2 <- gamVineCopSim(N, GVC, U))[3]
#'   print(all.equal(Sim1,Sim2))
#' }
#' 
#' # display results
#' # indeed, the R code is much slower
#' plot(Nn, t1, ylim = range(t1,t2), log = "xy", 
#'  xlab = "Number of observations", ylab = "Time elapsed")
#' points(Nn, t2, pch = 2)
#' legend("topleft", c("VineCopula", "gamVineCopula"), pch = 1:2)
#' 
#' @export
gamVineCopSim <- 
  function(N, GVC, U = NULL){
    
  stopifnot(N >= 1)
  if(any(!isS4(GVC),!is(GVC, "gamVineCop"))){
    stop("'GVC' has to be an gamVineCop object.")
  } 

	o <- diag(GVC@Matrix)
	d <- length(o)
	GVC <- gamVineCopNormalize(GVC)
  fam <- gamVineCopFamily(GVC)
  MaxMat=VineCopula:::createMaxMat(GVC@Matrix)
  CondDistr=VineCopula:::neededCondDistr(GVC@Matrix)

  model.count <- rep(0,d^2)
  temp <- 1:(d*(d-1)/2)
	t1 <- 1
  sel <- seq(d,d^2-d, by = d)
	for(i in 1:(d-1)){
	  t2 <- t1+d-i-1
	  model.count[sel[1:(d-i)]-i+1] <- temp[t1:t2]
	  t1 <- t2+1
	}
	model.count <- matrix(model.count,d,d)
  
	rotate <- function(x) t(apply(x, 2, rev))
	m <- rotate(rotate(GVC@Matrix))
	fam <- rotate(rotate(fam))
	model.count <- rotate(rotate(model.count))
	maxmat <- rotate(rotate(MaxMat))
	conindirect <- rotate(rotate(CondDistr$indirect))
		
	Vdirect <- Vindirect <- array(dim = c(d,d,N))
	if(is.null(U)){
		U <- matrix(runif(N*d), ncol = d)
	}	
	for(i in 1:d){
		Vdirect[i,i,] <- U[,i]
	}	
	Vindirect[1,1,] <- Vdirect[1,1,]	
	count <- 1
	for(i in 2:d){
		for(k in (i-1):1){
      #print(model.count[k,i])
      model <- GVC@model[[model.count[k,i]]]
			mm <- maxmat[k,i]
			if(mm == m[k,i]){
				zz = Vdirect[k,mm,]
			}else{
				zz = Vindirect[k,mm,]
			}

      if(isS4(model) && is(model, "gamBiCop")){
        vars <- all.vars(model@model$pred.formula)
        nvars <- length(vars)
        if(nvars == 1){
          newdata <- data.frame(Vdirect[1,1:nvars,])
          names(newdata) <- vars
        }else{
          newdata <- data.frame(t(Vdirect[1,1:nvars,]))
            names(newdata) <- vars
        } 
        par <- gamBiCopPred(model, newdata, target = "par")$par
        par2 <- model@par2
      }else{
        par <- rep(model$par, N)
        par2 <- rep(model$par2, N)
      }

			tmp <- rep(0, N)
    	tmp <- sapply(1:N, function(x) .C("Hinv1", 
    			as.integer(fam[k,i]), 
    			as.integer(1), 
       			as.double(Vdirect[k+1,i,x]),
       			as.double(zz[x]),
       			as.double(par[x]),
       			as.double(par2), 
       			as.double(tmp[x]),
       			PACKAGE = "VineCopula")[[7]])
       Vdirect[k,i,] <- tmp

       		if(i < d){
       			if(conindirect[k+1,i] == TRUE){
					tmp <- sapply(1:N, function(x) VineCopula::BiCopHfunc(Vdirect[k,i,x], 
									zz[x], fam[k,i], par[x], par2)$hfunc1)
					Vindirect[k+1,i,] <- tmp
       			}
       		}
		}
	}
 
 	out <- t(Vdirect[1,,])
	if(!is.null(GVC@names)){
		colnames(out) = GVC@names
	}
	out = out[,sort(o[length(o):1],index.return=TRUE)$ix]

	return(out)
}