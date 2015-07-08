# set.seed(1)
# 
# ##  Simulation parameters
# # Sample size
# n <- 2e2  
# # Copula family
# familyset <- c(1:4,13,14,23,24,33,34)
# # Define a 5-dimensional R-vine tree structure matrix
# d <- 4
# Matrix <- c(2,3,4,1,0,3,4,1,0,0,4,1,0,0,0,1)
# Matrix <- matrix(Matrix,d,d)
# nnames <- paste("X", 1:d, sep = "")
# 
# # Define gam-vine model list
# count <- 1
# model <- vector(mode = "list", length = d*(d-1)/2)
# par <- matrix(0,d,d)
# sel <- seq(d,d^2-d, by = d)
# 
# # First tree
# for (i in 1:(d-1)) {
#   # Select a copula family
#   family <- sample(familyset, 1)
#   model[[count]]$family <- family
#   
#   # Use the canonical link and a randomly generated parameter 
#   model[[count]]$par <- links(family)$par(rnorm(1))
#   if (family == 2) {
#     model[[count]]$par2 <- links(family)$par2(rnorm(1))
#   } else {
#     model[[count]]$par2 <- 0
#   }
#   count <- count + 1
# }
# 
# # Create a dummy dataset
# data <- data.frame(u1 = runif(1e2), u2 = runif(1e2), matrix(runif(1e2*d),1e2,d))
# 
# # Trees 2 to (d-1)
# for(j in 2:(d-1)){
#   for(i in 1:(d-j)){ 
#     # Select a copula family
#     family <- sample(familyset, 1)  
#     
#     # Select the conditiong set and create a model formula
#     cond <- nnames[sort(Matrix[(d-j+2):d,i])]
#     form <- as.formula(paste("~",paste(paste("s(", cond, ", k=10, bs='cr')",
#                                              sep = ""), collapse=" + ")))
#     
#     # Create a dummy gamBiCop object
#     tmp <- gamBiCopEst(data = data, formula = form, family = 1, n.iters = 1)$res
#     
#     # Update the copula family and the model coefficients to make them random
#     ll <- length(attr(tmp, "model")$coefficients)
#     attr(tmp, "model")$coefficients <- rnorm(ll)
#     attr(tmp, "family") <- family
#     if (family == 2) {
#       attr(tmp, "par2") <- links(family)$par2(rnorm(1))
#     }
#     model[[count]] <- tmp
#     count <- count+1  
#   } 
# }
# 
# # Create gamVineCopula object
# GVC <- gamVine(Matrix=Matrix,model = model,names=nnames)
# print(GVC)
# 
# # simulate the gamVineMatrix
# N <- 1e3
# sim <- gamVineSim(N, GVC)
# fitGVC <- gamVineSeqEst(sim, GVC, verbose = TRUE)
# # 
# par(mfrow=c(2,4))
# plot(GVC, ylim = c(-2.5,2.5))
# 
# plot(fitGVC, ylim = c(-2.5,2.5))
# 
# test <- gamVineStructureSelect(sim, verbose = TRUE)