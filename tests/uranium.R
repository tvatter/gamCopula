# library(gamCopula)
# library(copula)
# data(uranium)
# U <- pobs(uranium[,c(3,6,7)])
# fit <- gamVineStructureSelect(U)
# library(VineCopula)
# fitsimpl <- RVineStructureSelect(U)
# 
# ## The gamBiCop object you're looking for
# (C13_2 <- attr(fit,"model")[[3]])
# summary(C13_2)
# 
# ## par2
# attr(C13_2,"par2")
# 
# ## An example of function to extract the tau:
# ## a) Use the gamBiCopPred
# ## b) Take care of the name of the conditioning variable
# ## See ?gamBiCopPred
# tau13_2 <- function(u, alpha = NA) {
#   uu <- data.frame(Sc = u)  
#   
#   if (is.na(alpha)) {
#     out <- gamBiCopPred(C13_2,newdata = uu, target = "tau")$tau
#   } else {
#     tmp <- gamBiCopPred(C13_2,newdata = uu, target = "tau", alpha = alpha)
#     out <- list(tau = tmp$tau, CI = tmp$tau.CI)
#   }
# }
# 
# ## Plot only the estimated tau
# u <- seq(0,1,1e-3)
# plot(u,tau13_2(u), type = "l")
# 
# ## Plot with asymptotically valide confidence intervals (a few seconds)
# tmp <- tau13_2(u, 0.95)
# plot(u, tmp$tau, type = "l", ylim = range(tmp$CI), ylab = "tau")
# lines(u, tmp$CI[,1], lty = 2)
# lines(u, tmp$CI[,2], lty = 2)
# 
# ## Using the original data with the x-scale for the x-axis
# ## Similar as Fig. 9 from "Beyond simplified pair-copula constructions"
# tmp <- sort(uranium[,6], index.return = TRUE)
# x <- tmp$x
# tmp <- tau13_2(U[tmp$ix,2], 0.95)
# par(pty = "s")
# plot(x, tmp$tau, type = "l", ylim = c(-0.5,1), ylab = "tau")
# lines(range(x), rep(BiCopPar2Tau(1,fitsimpl$par[2,1]),2), col = "red")
# lines(x, tmp$CI[,1], lty = 2)
# lines(x, tmp$CI[,2], lty = 2)
# 
# ## With the whole dataset and truncation after level 3 
# ## (automatic plot could be improved)
# U <- pobs(uranium)
# model <- gamVineStructureSelect(U, trunclevel = 3, verbose = TRUE)
# plot(model)
# 
# ## test comment