# gamVineSeqEst<-function(dataset,GVC, method = "FS", tol.rel = 1e-3, n.iters = 10, verbose = FALSE)
# {
# 	data=as.matrix(data)
# 	if(any(data>1) || any(data<0)) stop("Data has be in the interval [0,1].")
# 	n = dim(GVC)
# 	N = nrow(data)
# 	if(dim(data)[2] != dim(GVC)) stop("Dimensions of 'data' and 'GVC' do not match.")
#   if(N < 2) stop("Number of observations has to be at least 2.")
#   if(!is(GVC,"gamVineCopula")) stop("'GVC' has to be an gamVineCopula object.")
# 
#   if(method!="FS" && method!="NR") stop("Estimation method has to be either 'FS' or 'NR'.")
#   if(is.logical(se)==FALSE) stop("'se' has to be a logical variable (TRUE or FALSE).")
#                                               
#   if(max.df<=1) stop("The upper bound for the degrees of freedom parameter has to be larger than 1.")
#     
#   o = diag(GVC$Matrix)
# 
#   oldGVC = GVC
# 
# 	if(any(o != length(o):1)){	
# 	 GVC = normalizegamVineCopula(GVC)
# 	 data = data[,o[length(o):1]]
#   }
# 
# 	Params=GVC$par
# 	Params2=GVC$par2
# 
# 	# if(se==TRUE)
# 	  # {
# 		# seMat1 = matrix(0,nrow=n,ncol=n)
# 		# seMat2 = matrix(0,nrow=n,ncol=n)
# 	  # }
# 
# 	V = list()
# 	V$direct = array(NA,dim=c(n,n,N))
# 	V$indirect = array(NA,dim=c(n,n,N))
# 		
# 	V$direct[n,,] = t(data[,n:1])
# 
# 	for(i in (n-1):1)
# 	{
# 		
# 		for(k in n:(i+1))
# 		{
# 			
# 			m = GVC$MaxMat[k,i]	
# 			zr1 = V$direct[k,i,]
# 				
# 			if(m == GVC$Matrix[k,i]){
# 				zr2 = V$direct[k,(n-m+1),]
# 			}else{
# 				zr2 = V$indirect[k,(n-m+1),]
# 			}
# 		
# 
# 			if(GVC$family[k,i] %in% c(2,7,8,9,10,17,18,19,20,27,28,29,30,37,38,39,40))
# 			{
# 				if(progress == TRUE){
#           if(k == n) message(oldGVC$Matrix[i,i],",",oldGVC$Matrix[k,i])
#           else message(oldGVC$Matrix[i,i],",",oldGVC$Matrix[k,i],"|",paste(oldGVC$Matrix[(k+1):n,i],collapse=","))
#         }
#         par.out <- gamBiCopEst(zr2,zr1,GVC$family[k,i], method, se, max.df, weights)
# 				#par1 <- out.par$par
# 				Params[k,i] <- par.out$model
# 				Params2[k,i] <- par.out$par2
# 				# if(se==TRUE)
# 				# {
# 				# #se1 <- par.out$se
# 				# seMat1[k,i] <- par.out$se
# 				# seMat2[k,i] <- par.out$se2
# 				# }
# 			}else{
# 				if(progress == TRUE){
#           if(k == n) message(oldGVC$Matrix[i,i],",",oldGVC$Matrix[k,i])
#           else message(oldGVC$Matrix[i,i],",",oldGVC$Matrix[k,i],"|",paste(oldGVC$Matrix[(k+1):n,i],collapse=","))
#         }
#         par.out <- gamBiCopEst(zr2,zr1,GVC$family[k,i], method, se, max.df, max.BB,weights)
# 				#par1 <- out.par$par
# 				Params[k,i] <- par.out$model
# 				Params2[k,i] <- par.out$par2
# 				# if(se==TRUE)
# 				# {
# 				# #se1 <- par.out$se
# 				# seMat1[k,i] <- par.out$se
# 				# seMat2[k,i] <- par.out$se2
# 				# }
# 			}
# 			
# 				
# 				if(GVC$CondDistr$direct[k-1,i])
# 				{
# 					
# 					V$direct[k-1,i,] = sapply(par.out$par, function(par) .C("Hfunc1",
# 					   as.integer(GVC$family[k,i]),
# 					   as.integer(length(zr1)),
# 					   as.double(zr1),
# 					   as.double(zr2),
# 					   as.double(par),
# 					   as.double(Params2[k,i]),
# 					   as.double(rep(0,length(zr1))),
# 					   PACKAGE='VineCopula')[[7]])
# 				}
# 				if(GVC$CondDistr$indirect[k-1,i])
# 				{
# 					V$indirect[k-1,i,] = sapply(par.out$par, function(par) .C("Hfunc2",
# 					   as.integer(GVC$family[k,i]),
# 					   as.integer(length(zr2)),
# 					   as.double(zr2),
# 					   as.double(zr1),
# 					   as.double(par),
# 					   as.double(Params2[k,i]),
# 					   as.double(rep(0,length(zr1))),
# 					   PACKAGE='VineCopula')[[7]])
# 				}
# 			
# 		}
# 	}
# 
#   oldGVC$par = Params
#   oldGVC$par2 = Params2
# 
# if(se==FALSE)
# 	return(list(GVC=oldGVC))
# # else
# 	# return(list(GVC=oldGVC, se=seMat1, se2=seMat2))
# }
