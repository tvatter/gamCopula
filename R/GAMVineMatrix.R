gamVineMatrix  <- 
  function(Matrix,model=vector(mode = "list",length=dim(Matrix)[1]*(dim(Matrix)[1]-1)/2),names=NULL){
	Matrix[is.na(Matrix)]=0  
  d <- dim(Matrix)[1]
	if(dim(Matrix)[2] != d) stop("Structure matrix has to be quadratic.")
	if(max(Matrix)>d) stop("Error in the structure matrix.")
	if(RVineMatrixCheck(Matrix)!=1) stop("'Matrix' is not a valid R-vine matrix")
	if(length(names)>0 & length(names)!= d){ 
    stop("Length of the vector 'names' is not correct.")
  }else if(length(names) == 0){
        names <- paste("x", 1:d, sep = "")
  }
	
	
  if(length(model)!= d*(d-1)/2){ 
    stop("Length of the list 'model' is not correct.") 
  }
  count <- 1
	# First tree
	for(i in 1:(d-1)){
	  mm <- model[[count]]
	  if(!is.numeric(mm)){
	    if(!validgamBiCop(mm)){
	      stop(paste("Element", count, "of the model list should 
                       be a valid gamBiCop object or 0."))
	    }
	    cond <- sort(all.vars(model$pred.formula))
	    if(length(cond) != 0){
        print(cond)
	        stop(paste("The formula of element ", count, " of the model list should
                         not contain conditioning variables.", sep = ""))
	    }
	  }else if(mm != 0){
	    stop(paste("Element", count, "of the model list should 
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
	        stop(paste("Element ", count," of the model list should 
	                   be a valid gamBiCop object or 0.", sep = ""))
	      }
	      cond <- sort(all.vars(model$pred.formula))
	      cond2 <- names[sort(Matrix[(d-j+2):d,i])]
	      if(!all(cond == cond2)){
	         stop(paste("The formula of element ", count, " of the model list should
	                    not contain conditioning variables.", sep = ""))}
	    }else if(mm != 0){
	      stop(paste("Element", count, "of the model list should 
                       be a valid gamBiCop object or 0."))
	    } 
	    count <- count+1  
	  } 
	}

	MaxMat=VineCopula:::createMaxMat(Matrix)
	CondDistr=VineCopula:::neededCondDistr(Matrix)
					
	GVM = list(Matrix = Matrix, model = model, names = names, MaxMat = MaxMat, CondDistr = CondDistr)
			   
	class(GVM) = "gamVineMatrix"
	return(GVM)	
}

normalizegamVineMatrix = function(GVM){
	
	oldOrder = diag(GVM$Matrix)
	Matrix = VineCopula:::reorderRVineMatrix(GVM$Matrix)
	
  return(gamVineMatrix(Matrix, GVM$model, names = rev(GVM$names[oldOrder])))
}

# exported version of normalizegamVineMatrix
gamVineMatrixNormalize <- function(GVM) {
  stopifnot(is(GVM, "gamVineMatrix"))
      
  if(is.null(GVM$names))
    GVM$names <- paste("V",1:nrow(GVM$Matrix),sep="")
      oldOrder = diag(GVM$Matrix)
  
  return(normalizegamVineMatrix(GVM))
}

dim.gamVineMatrix = function(x){
	gamVine=x
	return(dim(gamVine$Matrix)[1])
	NextMethod("dim")
}

print.gamVineMatrix = function(x, detail=FALSE, ...){
	gamVine=x
	NextMethod("print")
	message("gam-vine matrix:")
	print(gamVine$Matrix, ...)
	
	if(!is.null(gamVine$names)){
		message("")
		message("Where")
		for(i in 1:length(gamVine$names)){
			message(i," <-> ",gamVine$names[[i]])
		}
	}

	d=dim(gamVine)
  count <- 1
	if(detail==TRUE || detail==T)
	{
		message("")
		message("Tree 1:")
		for(i in 1:(d-1))
		{
			  a <- paste(gamVine$names[[gamVine$Matrix[i,i]]],",",gamVine$names[[gamVine$Matrix[d,i]]],sep="")
        if(is.numeric(gamVine$model[[count]])){
          a <- paste(a,": ", BiCopName(gamVine$model[[count]], short=FALSE),sep="")
        }else{
          a <- paste(a,": ", BiCopName(gamVine$model[[count]]@family, short=FALSE),sep="")
        }		
   			message(a)
# 			if(gamVine$model[d,i]!=0)
# 			{
# 				show(gamVine$model[d,i])
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
				a <- paste(gamVine$names[[gamVine$Matrix[i,i]]],",",gamVine$names[[gamVine$Matrix[d-j+1,i]]],sep="")
				a <- paste(a,"|",sep="")
				conditioningSet=(d-j+2):d
				for(k in 1:length(conditioningSet))
				{
					if(k>1){ a <- paste(a,",",sep="")}
					a <- paste(a,gamVine$names[[gamVine$Matrix[conditioningSet[k],i]]],sep="")
				}       
				if(is.numeric(gamVine$model[[count]])){
				  a <- paste(a,": ", BiCopName(gamVine$model[[count]], short=FALSE),sep="")
				}else{
          EDF <- EDF.gamBiCop(gamVine$model[[count]])
				  a <- paste(a,": ", BiCopName(gamVine$model[[count]]@family, short=FALSE),sep="")
          a <- paste(a, paste("EDF:", paste(round(EDF[-1],3), collapse = ", ")), sep = "\n")
				}
				message(a)
# 				if(gamVine$model[d-j+1,i]!=0)
# 				{
# 				  if(gamVine$model[d-j+1,i]!=0)
# 				  {
# 				    show(gamVine$model[d-j+1,i])
# 				  }
# 				}
        count <- count+1
			}
		}

	}
}

family.gamVineMatrix = function(x){
  GVM <- x
  
  d <- dim(GVM)
  fam <- rep(0,d^2)
  temp <- sapply(GVM$model, function(x) x@family)
  sel <- seq(d,d^2-d, by = d)
  t1 <- 1
  for(i in 1:(d-1)){
    t2 <- t1+d-i-1
    fam[sel[1:(d-i)]-i+1] <- temp[t1:t2]
    t1 <- t2+1
  }
  
  return(matrix(fam,d,d))
  NextMethod("family")
}

# as.gamVineMatrix = function(gamVine){
# 
# 	n = length(gamVine$Tree)+1
# 	con = list()
# 	nam = V(gamVine$Tree[[1]])$name
# 	
# 	conditionedSets = NULL
# 	corresppondingModels = list()
# 	corresppondingTypes = list()
# 	
# 	print(is.list(E(gamVine$Tree[[n-1]])$conditionedSet))
# 	
# 	conditionedSets[[n-1]][[1]] = (E(gamVine$Tree[[n-1]])$conditionedSet)
# 	for(k in 1:(n-2)){
# 		conditionedSets[[k]] = E(gamVine$Tree[[k]])$conditionedSet
# 		corresppondingModels[[k]] = as.list(E(gamVine$Tree[[k]])$model)
# 	}
# 	corresppondingModels[[n-1]] = list()
# 	corresppondingModels[[n-1]][[1]] = (E(gamVine$Tree[[n-1]])$model)
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
# 	return(gamVineMatrix(M, family = model, names = nam))
# 	
# }
