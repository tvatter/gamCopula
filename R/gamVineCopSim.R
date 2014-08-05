#' Simulation from a \code{\link{gamVineCop-class}} object
#'
#' @param N number of d-dimensional observations to simulate.
#' @param GVC fitted \code{\link{gamVineCop-class}} object.
#' @param U (similar as \code{\link{RVineSim}} from the from the \code{\link{VineCopula}} package) 
#' If not NULL, an (N,d)-matrix of U[0,1] random variates to be transformed to the copula sample.
#' @return A Nxd matrix of data simulated from the given \code{\link{gamVineCop-class}} object.
#' @export
gamVineCopSim <- 
  function(N, GVC, U = NULL)
{
  stopifnot(N >= 1)
  if(!is(GVC, "gamVineCop")) stop("'GVC' has to be an gamVineCop object.")

	o <- diag(GVC$Matrix)
	d <- length(o)
	GVC <- gamVineCopNormalize(GVC)
  fam <- gamVineCopFamily(GVC)
  
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
	m = rotate(rotate(GVC$Matrix))
	fam = rotate(rotate(fam))
	model.count = rotate(rotate(model.count))
	maxmat = rotate(rotate(GVC$MaxMat))
	conindirect = rotate(rotate(GVC$CondDistr$indirect))
		
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
      model <- GVC$model[[model.count[k,i]]]
			mm <- maxmat[k,i]
			if(mm == m[k,i]){
				zz = Vdirect[k,mm,]
			}else{
				zz = Vindirect[k,mm,]
			}
      if(k == 1){
        par <- rep(gamBiCopPred(model, target = "par")$par[1],N)
        par2 <- model@par2
      }else{
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
					tmp <- sapply(1:N, function(x) BiCopHfunc(Vdirect[k,i,x], 
									zz[x], fam[k,i], par[x], par2)$hfunc1)
					Vindirect[k+1,i,] <- tmp
       			}
       		}
		}
	}

 	out <- t(Vdirect[1,,])
	if(!is.null(GVC$names)){
		colnames(out) = GVC$names
	}
	out = out[,sort(o[length(o):1],index.return=TRUE)$ix]

	return(out)
}