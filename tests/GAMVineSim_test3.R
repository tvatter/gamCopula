# set.seed(0)
# 
# ##  Simulation parameters
# # Sample size
# n <- 2e2  
# # Copula family
# fam <- 1:4
# # Degrees of freedom (for the t-copula when fam <- 2)
# par2 <- 4 
# # Should the model be specified in terms of Kendall's tau (TRUE) or copula parameter
# tau <- FALSE
# # Newton-Raphse ("NR) or Fisher-Scoring ("FS") algorithm
# met <- "FS"
# # Relative tolerance for "NR"/"FS"
# tol.rel <- 1e-6
# # Max number of iterations for "NR"/"FS"
# n.iters <- 25
# # Define a 5-dimensional R-vine tree structure matrix
# d <- 5
# Matrix <- c(5,2,3,1,4,0,2,3,4,1,0,0,3,4,1,0,0,0,4,1,0,0,0,0,1)
# Matrix <- matrix(Matrix,d,d)
# nnames <- paste("X", 1:5, sep = "")
# 
# ## Integration grid
# step <- 1e-3
# ngrid <- 1/step
# xx <- seq(0,1,length.out = ngrid)
# Xx <- data.frame(cbind(xx,xx,xx,xx,xx))
# names(Xx) <- nnames
# true <- true.approx <- matrix(NA, ngrid, dim(Xx)[2])
# 
# ## Base level
# eta0 <- 1
# 
# # define gam-vine model list
# count <- 1
# model <- vector(mode = "list", length = d*(d-1)/2)
# par <- matrix(0,d,d)
# sel <- seq(d,d^2-d, by = d)
# # First tree
# for(i in 1:(d-1))
# {
#   print(count)
#   # Select copula family
#   family <- sample(fam, 1)
#   
#   # Simulate and fit some data
#   temp <- BiCopSim(n, family, par = BiCopEta2Par(family, eta0), par2 = par2)
#   dataset <- data.frame("u1" = temp[,1], "u2" = temp[,2])      
#   model[[count]] <- gamBiCopEst(family = family, parFrhs = ~1, dataset, tau = tau,
#                                 method = met, tol.rel = tol.rel, n.iters = n.iters)$res
#   par[sel[i]] <- BiCopEta2Par(family,model[[count]]@model$coefficient[1])
#   count <- count + 1
# }
# 
# # Trees 2 to (d-1)
# for(j in 2:(d-1)){
#   for(i in 1:(d-j)){ 
#     # Select a copula family
#     family <- sample(fam, 1)  
#     
#     # Select a true underlying smooth component for each conditioning variable  
#     cond <- nnames[sort(Matrix[(d-j+2):d,i])]
#     l <- length(cond)
#     temp <- rep(0,3)
#     temp[sample(3, l)] <- 1
#     
#     # Simulate a dataset
#     if(l != 1){
#       covariates.distr <- copula:::mvdc(copula:::normalCopula(rho, dim = l), c("unif"), list(list(min = 0, max = Tf)), marginsIdentical = TRUE) 
#       X <- copula:::rMvdc(n,covariates.distr)
#       formula.expr <- c()
#       for(o in 1:l){
#         formula.expr <- c(formula.expr,paste(cond[o], sep = ""))
#       }
#       formula.expr <- paste(formula.expr, collapse = " + ")
#       formula.expr <-  paste("~", formula.expr, collapse = "", sep = "")  
#       #print(formula.expr)
#     }else{
#       X <- runif(n)
#       formula.expr <- paste("~",cond, sep = "")
#       #print(formula.expr)
#     }
#     formula.expr <- as.formula(formula.expr)
#     
#     if(l == 1){
#             calib.temp <- function(x){eta0}
#     }else if(l == 2){
#       calib.temp <- function(x1,x2){eta0}
#     }else{
#       calib.temp <- function(x1,x2,x3){eta0}
#     }
#     temp <- CondBiCopSim(family, calib.temp, X, par2=par2, return.par=FALSE, tau = tau) 
#     dataset <- data.frame("u1" = temp[,1], "u2" = temp[,2])
#     dataset <- cbind(dataset, X)
#     names(dataset)[3:(2+l)] <- cond   
#     
#     # Estimate the gam model
#     model[[count]] <- gamBiCopEst(family = family, parFrhs = formula.expr, dataset, tau = tau,
#                                   method = met, tol.rel = tol.rel, n.iters = n.iters)$res 
#     par[sel[i]-j+1] <- BiCopEta2Par(family,model[[count]]@model$coefficient[1])
#     count <- count+1  
#   } 
# }
# # define gamVineMatrix object
# GVM <- gamVineMatrix(Matrix=Matrix,model = model,names=nnames)
# print(GVM)
# fam <- family(GVM)
# par22 <- fam
# par22[which(par22 != 2)] <- 0
# par22[which(par22 == 2)] <- par2
# 
# RVM = RVineMatrix(Matrix=Matrix,family=fam,par=par,par2=par22,
#                   names=c("X1","X2","X3","X4","X5"))
# Nn <- rep(c(1,2,5), 4)*10^c(rep(1,3),rep(2,3), rep(3,3), rep(4,3))
# t1 <- t2 <- rep(0, length(Nn))
# for(i in 1:length(Nn)){
#   N <- Nn[i]
#   U <- matrix(rep(0.5,N*5), ncol = 5)
#   t1[i] <- system.time(Sim1 <- RVineSim(N, RVM, U))[3]
#   t2[i] <- system.time(Sim2 <- gamVineSim(N, GVM, U))[3]
#   print(all.equal(Sim1,Sim2,1e-3))
# }
# 
# par(mfrow = c(1,2))
# plot(Nn, t1, ylim = range(t1,t2), log = "xy")
# points(Nn, t2, pch = 2)
# plot(Nn, t2/t1, log = "xy")