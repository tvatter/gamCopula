
#' Normalize a \code{\link{gamVineCop-class}} object
#' 
#' Change the R-vine matrix in the natural order, 
#' i.e. with d:1 on the diagonal
#'
#' @param GVC fitted \code{\link{gamVineCop-class}} object.
#' @return Normalized \code{\link{gamVineCop-class}} object.
#' @seealso \code{\link{gamVineCop-class}} and \code{\link{gamVineCop}}.
#' @export
gamVineCopNormalize <- function(GVC) {
  if(any(!isS4(GVC),!is(GVC, "gamVineCop"))) stop("'GVC' has to be an gamVineCop object.")
  
  oldOrder = diag(GVC@Matrix)
  Matrix = VineCopula:::reorderRVineMatrix(GVC@Matrix)
  
  return(gamVineCop(Matrix, GVC@model, names = rev(GVC@names[oldOrder])))
}

#' Dimension of a \code{\link{gamVineCop-class}} object
#' 
#' Retrieve the dimension of a \code{\link{gamVineCop-class}} object.
#'
#' @param x fitted \code{\link{gamVineCop-class}} object.
#' @return Dimension of the \code{\link{gamVineCop-class}} object.
#' @seealso \code{\link{gamVineCop-class}} and \code{\link{gamVineCop}}.
#' @docType methods
#' @rdname dim-methods
#' @export
dim.gamVineCop = function(x){
  GVC=x
  return(dim(GVC@Matrix)[1])
}
#' @docType methods
#' @rdname dim-methods
setMethod("dim", signature("gamVineCop"), dim.gamVineCop)

#' Print a \code{\link{gamVineCop-class}} object
#'
#' @param x \code{\link{gamVineCop-class}} object.
#' @param ... un-used for this class.
## @param detail should additional details be printed (\code{detail=TRUE}) or not?
#' @seealso \code{\link{gamVineCop-class}} and \code{\link{gamVineCop}}.
#' @docType methods
#' @rdname print-methods
#' @export
print.gamVineCop = function(x, ...){
  ## TODO , detail=TRUE
  detail = FALSE
  GVC=x
  message("gam-vine matrix:")
  print(GVC@Matrix, ...)
  message("")
  message("Where")
  for(i in 1:length(GVC@names)){
    message(i," <-> ",GVC@names[[i]])
  }
  
  d=dim(GVC)
  count <- 1
  if(detail==TRUE || detail==T)
  {
    message("")
    message("Tree 1:")
    for(i in 1:(d-1))
    {
      a <- paste(GVC@names[[GVC@Matrix[i,i]]],",",GVC@names[[GVC@Matrix[d,i]]],sep="")
      if(is.numeric(GVC@model[[count]])){
        a <- paste(a,": ", BiCopName(GVC@model[[count]], short=FALSE),sep="")
      }else{
        a <- paste(a,": ", BiCopName(GVC@model[[count]]@family, short=FALSE),sep="")
      }  	
      message(a)
      # 			if(GVC@model[d,i]!=0)
      # 			{
      # 				print(GVC@model[d,i])
      # 			}
      count <- count+1
    }
    for(j in 2:(d-1))
    {
      message("")
      a <- paste("Tree ",j,":",sep="")
      message(a)
      for(i in 1:(d-j))
      {
        a <- paste(GVC@names[[GVC@Matrix[i,i]]],",",GVC@names[[GVC@Matrix[d-j+1,i]]],sep="")
        a <- paste(a,"|",sep="")
        conditioningSet=(d-j+2):d
        for(k in 1:length(conditioningSet))
        {
          if(k>1){ a <- paste(a,",",sep="")}
          a <- paste(a,GVC@names[[GVC@Matrix[conditioningSet[k],i]]],sep="")
        }       
        if(is.numeric(GVC@model[[count]])){
          a <- paste(a,": ", BiCopName(GVC@model[[count]], short=FALSE),sep="")
        }else{
          EDF <- EDF.gamBiCop(GVC@model[[count]])
          a <- paste(a,": ", BiCopName(GVC@model[[count]]@family, short=FALSE),sep="")
          a <- paste(a, paste("EDF:", paste(round(EDF[-1],3), collapse = ", ")), sep = "\n")
        }
        message(a)
        # 				if(GVC@model[d-j+1,i]!=0)
        # 				{
        # 				  if(GVC@model[d-j+1,i]!=0)
        # 				  {
        # 				    print(GVC@model[d-j+1,i])
        # 				  }
        # 				}
        count <- count+1
      }
    }
    
  }
}
#' @docType methods
#' @rdname print-methods
setMethod("print", signature("gamVineCop"), print.gamVineCop)

#' Family matrix of \code{\link{gamVineCop-class}} object
#' 
#' Return the matrix of copula family (see \code{\link{gamBiCop-class}}) corresponding 
#' to the model list in the \code{\link{gamVineCop-class}} object.
#'
#' @param GVC \code{\link{gamVineCop-class}} object.
#' @return Matrix of copula family (see \code{\link{gamBiCop-class}}) corresponding 
#' to the model list in the \code{\link{gamVineCop-class}} object.
#' @seealso \code{\link{gamVineCop-class}} and \code{\link{gamVineCop}}.
#' @export
gamVineCopFamily = function(GVC){
  
  if(any(!isS4(GVC),!is(GVC, "gamVineCop"))) stop("'GVC' has to be an gamVineCop object.")
  
  d <- dim(GVC)
  fam <- rep(0,d^2)
  temp <- sapply(GVC@model, function(x) if(isS4(GVC) && is(GVC, "gamBiCop")){x@family}else{x$family})
  sel <- seq(d,d^2-d, by = d)
  t1 <- 1
  for(i in 1:(d-1)){
    t2 <- t1+d-i-1
    fam[sel[1:(d-i)]-i+1] <- temp[t1:t2]
    t1 <- t2+1
  }
  
  return(matrix(fam,d,d))
}

#' \code{\link{RVineMatrix}} to \code{\link{gamVineCop-class}}
#' 
#' Transform a \code{\link{RVineMatrix}} object 
#' into a \code{\link{gamVineCop-class}} object.
#'
#' @param RVM \code{\link{RVineMatrix}} object.
#' @return A \code{\link{gamVineCop-class}} object.
#' @seealso \code{\link{RVineMatrix}}, \code{\link{gamVineCop-class}} and \code{\link{gamVineCop}}.
#' @export
RVM2GVC <- function(RVM){
  
  if(class(RVM) != "RVineMatrix"){
    stop("RVM has to be a RVineMatrix class object.")
  }
  
  d <- dim(RVM$Matrix)[1]
  
  sel <- seq(d,d^2-d, by = d)
  family <- RVM$family[sel]
  par <-  RVM$par[sel]
  par2 <-  RVM$par2[sel]
  for(j in 2:(d-1)){
    family <- c(family, RVM$family[sel[1:(d-j)]-j+1])
    par <- c(par, RVM$par[sel[1:(d-j)]-j+1])
    par2 <- c(par2, RVM$par2[sel[1:(d-j)]-j+1])
  }
  
  out <- cbind(family, par, par2)
  colnames(out) <- c()
  out <- apply(out, 1, function(x) list(family = x[1], par = x[2], par2 = x[3]))
  if(!is.null(RVM$names)){
    out <- gamVineCop(RVM$Matrix,out,RVM$names)
  }else{
    nnames <- paste("X", 1:d, sep = "")
    out <- gamVineCop(RVM$Matrix,out,nnames)
  }
  return(out)
}

# as.gamVineCop = function(gamVine){
# 
# 	n = length(GVC@Tree)+1
# 	con = list()
# 	nam = V(GVC@Tree[[1]])$name
# 	
# 	conditionedSets = NULL
# 	corresppondingModels = list()
# 	corresppondingTypes = list()
# 	
# 	print(is.list(E(GVC@Tree[[n-1]])$conditionedSet))
# 	
# 	conditionedSets[[n-1]][[1]] = (E(GVC@Tree[[n-1]])$conditionedSet)
# 	for(k in 1:(n-2)){
# 		conditionedSets[[k]] = E(GVC@Tree[[k]])$conditionedSet
# 		corresppondingModels[[k]] = as.list(E(GVC@Tree[[k]])$model)
# 	}
# 	corresppondingModels[[n-1]] = list()
# 	corresppondingModels[[n-1]][[1]] = (E(GVC@Tree[[n-1]])$model)
# 	
# 	model = array(dim=c(n,n))
# 	M = matrix(NA,n,n)
# 	
# 	for(k in 1:(n-1)){
# 		w = conditionedSets[[n-k]][[1]][1]
# 		
# 		M[k,k] = w
# 		M[(k+1),k] = conditionedSets[[n-k]][[1]][2]
# 		
# 		model[(k+1),k] = corresppondingModels[[n-k]][[1]][1]
# 		
# 		if(k == (n-1)){
# 			M[(k+1),(k+1)] = conditionedSets[[n-k]][[1]][2]
# 		}else{
# 			for(i in (k+2):n){
# 				for(j in 1:length(conditionedSets[[n-i+1]])){
# 					cs = conditionedSets[[n-i+1]][[j]]
# 					if(cs[1] == w){
# 						M[i,k] = cs[2]
# 						break
# 					} else if(cs[2] == w){
# 						M[i,k] = cs[1]
# 						break
# 					}
# 				}
# 				model[i,k] = corresppondingModels[[n-i+1]][[j]]
# # 				
# # 				conditionedSets[[n-i+1]][[j]] = NULL
# # 				corresppondingParams[[n-i+1]][[j]] = NULL
# # 				corresppondingTypes[[n-i+1]][[j]] = NULL
# 			}
# 		}
# 		
# 	}
# 	
# 	M = M+1 
# 	M[is.na(M)]=0
# 	Type[is.na(Type)]=0
# 	
# 	return(gamVineCop(M, family = model, names = nam))
# 	
# }
