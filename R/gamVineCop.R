#############################
##  gam vine copulas  ##
#############################
valid.gamVineCop = function(object) { 
  
  d <- length(attributes(object))
  if((d < 3) || names(attributes(object))[1:3] != c("Matrix", "model",  "names")){
    return("A gamVineCop contains at least a R-Vine matrix, a list of gamBiCop objects and 
           a vector of names.")
  }
  Matrix <- object@Matrix  
  names <- object@names
  Matrix[is.na(Matrix)]=0  
  d <- dim(Matrix)[1]
  if(dim(Matrix)[2] != d) return("Structure matrix has to be quadratic.")
  if(max(Matrix)>d) return("Error in the structure matrix.")
  if(RVineMatrixCheck(Matrix)!=1) return("'Matrix' is not a valid R-vine matrix")
  if(length(names)>1 & length(names)!= d){ 
    return("Length of the vector 'names' is not correct.")
  }else if(length(names) == 0){
    names <- paste("x", 1:d, sep = "")
  } 
  
  model <- object@model
  if(length(model)!= d*(d-1)/2){ 
    return("Length of the list 'model' is not correct.") 
  }
  count <- 1
  # First tree
  for(i in 1:(d-1)){
    mm <- model[[count]]
    if(!is.numeric(mm)){
      if(!validgamBiCop(mm)){
        return(paste("Element", count, "of the model list should 
                     be a valid gamBiCop object or 0."))
      }
      cond <- sort(all.vars(model$pred.formula))
      if(length(cond) != 0){
        print(cond)
        return(paste("The formula of element ", count, " of the model list should
                     not contain conditioning variables.", sep = ""))
      }
      }else if(mm != 0){
        return(paste("Element", count, "of the model list should 
                     be a valid gamBiCop object or 0."))
      } 
    count <- count + 1
      }
  
  # Trees 2 to (d-1)
  for(j in 2:(d-1)){
    for(i in 1:(d-j)){ 
      mm <- model[[count]]
      if(!is.numeric(mm)){
        if(!validgamBiCop(mm)){
          return(paste("Element ", count," of the model list should 
                     be a valid gamBiCop object or 0.", sep = ""))
        }
        cond <- sort(all.vars(model$pred.formula))
        cond2 <- names[sort(Matrix[(d-j+2):d,i])]
        if(!all(cond == cond2)){
          return(paste("The formula of element ", count, " of the model list should
                     not contain conditioning variables.", sep = ""))}
      }else if(mm != 0){
        return(paste("Element", count, "of the model list should 
                     be a valid gamBiCop object or 0."))
      } 
      count <- count+1  
    } 
  }
  
  return(TRUE)       	
}

#'  The \code{\link{gamVineCop-class}}
#'
#'  \code{\link{gamVineCop-class}} is an S4 class to store 
#'  a generalized additive model on a vine copula.
#'
#' @slot Matrix lower triangular d x d matrix that defines the R-vine tree structure.
#' @slot model list containing d x (d-1)/2 \code{\link{gamBiCop-class}} objects.
#' @slot names vector of d names.
#' @seealso \code{\link{gamVineCop}}, \code{\link{RVineMatrix}} and \code{\link{gamBiCop-class}}.
#' @export
setClass("gamVineCop",
         slots = c(Matrix="matrix", model = "list", names = "character"),
         validity = valid.gamVineCop
)

#' Constructor of the \code{\link{gamVineCop-class}}
#'
#'  A constructor for objects of the \code{\link{gamVineCop-class}}.
#'
#' @param Matrix lower triangular d x d matrix that defines the R-vine tree structure.
#' @param model list containing d x (d-1)/2 \code{\link{gamBiCop-class}} objects.
#' @param names vector of d names.
#' @return A \code{\link{gamVineCop-class}} object.
#' @seealso \code{\link{gamVineCop-class}}, \code{\link{RVineMatrix}} and \code{\link{gamBiCop-class}}.
gamVineCop  <- 
  function(Matrix,model,names=NA){
	#MaxMat=VineCopula:::createMaxMat(Matrix)
	#CondDistr=VineCopula:::neededCondDistr(Matrix)
  
  new("gamVineCop", Matrix = Matrix, model = model, names = as.character(names))
}

#' Normalize a \code{gamVineCop-class} object
#' 
#' Change the R-vine matrix in the natural order, 
#' i.e. with d:1 on the diagonal
#'
#' @param GVC fitted \code{\link{gamVineCop-class}} object.
#' @return Normalized \code{gamVineCop-class} object.
#' @seealso \code{gamVineCop-class} and \code{gamVineCop}.
#' @export
gamVineCopNormalize <- function(GVC) {
  stopifnot(is(GVC, "gamVineCop"))
      
  oldOrder = diag(GVC@Matrix)
  Matrix = VineCopula:::reorderRVineMatrix(GVC@Matrix)
  
  return(gamVineCop(Matrix, GVC@model, names = rev(GVC@names[oldOrder])))
}

#' Dimension of a \code{gamVineCop-class} object
#' 
#' Retrieve the dimension of a \code{gamVineCop-class} object.
#'
#' @param x fitted \code{\link{gamVineCop-class}} object.
#' @return Dimension of the \code{gamVineCop-class} object.
#' @seealso \code{gamVineCop-class} and \code{gamVineCop}.
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

#' Print a \code{gamVineCop-class} object
#' 
#' Print a \code{gamVineCop-class} object.
#'
#' @param x fitted \code{\link{gamVineCop-class}} object.
#' @param detail should additional details be printed (\code{detail=TRUE}) or not?
#' @seealso \code{gamVineCop-class} and \code{gamVineCop}.
#' @docType methods
#' @rdname print-methods
#' @export
print.gamVineCop = function(x, detail=FALSE, ...){
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
# 				show(GVC@model[d,i])
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
# 				    show(GVC@model[d-j+1,i])
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

#' Family matrix of \code{gamVineCop-class} object
#' 
#' Return the matrix of copula family (see \code{\link{gamBiCop-class}}) corresponding 
#' to the model list in the \code{gamVineCop-class} object.
#'
#' @param x fitted \code{\link{gamVineCop-class}} object.
#' @return Matrix of copula family (see \code{\link{gamBiCop-class}}) corresponding 
#' to the model list in the \code{gamVineCop-class} object.
#' @seealso \code{gamVineCop-class} and \code{gamVineCop}.
#' @export
gamVineCopFamily = function(x){
  GVC <- x
  
  d <- dim(GVC)
  fam <- rep(0,d^2)
  temp <- sapply(GVC@model, function(x) x@family)
  sel <- seq(d,d^2-d, by = d)
  t1 <- 1
  for(i in 1:(d-1)){
    t2 <- t1+d-i-1
    fam[sel[1:(d-i)]-i+1] <- temp[t1:t2]
    t1 <- t2+1
  }
  
  return(matrix(fam,d,d))
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
