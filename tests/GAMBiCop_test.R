# ##															
# ##	1) Simulate and plot data	
# ##	2) Estimate the model and display results
# ##															
# 
# library(gamVineCopula)
# set.seed(1)
# 
# ##  Simulation parameters
# # Sample size
# n <- 2e2  
# # Correlation between the covariates
# rho <- 0.5
# # Copula family (Clayton here)
# fam <- 3  
# # Degrees of freedom (for the t-copula when fam <- 2)
# par2 <- 4 
# # Should the model be specified in terms of Kendall's tau (TRUE) or copula parameter
# tau <- FALSE
# # Newton-Raphse ("NR) or Fisher-Scoring ("FS") algorithm
# met <- "FS"
# # Relative tolerance for "NR"/"FS"
# tol <- 1e-6
# # Max number of iterations for "NR"/"FS"
# itermax <- 25
# 
# ##
# ## 1) Simulate and plot data
# ##
# 
# ## Integration grid
# step <- 1e-3
# ngrid <- 1/step
# xx <- seq(0,1,length.out = ngrid)
# Xx <- data.frame(cbind(xx,xx,xx))
# names(Xx) <- c("x1", "x2", "x3")
# true <- true.approx <- true.approx2 <- matrix(NA, ngrid, dim(Xx)[2])
# 
# ## Quadratic calibration
# Tf <- 1
# b.quad <- 8*Tf
# t0.quad <- Tf/2
# a.quad <- -(b.quad/3)*(Tf^2-3*Tf*t0.quad+3*t0.quad^2)
# calib.quadratic <- function(t, t0, a, b){return(a + b*(t-t0)^2)}
# true[,1] <- calib.quadratic(xx, t0.quad, a.quad, b.quad)
# 
# ## Sinusoidal calibration
# b.sin <- 1
# f.sin <- 1
# t0.sin <- 0
# a.sin <- b.sin*(1-2*f.sin*Tf*pi/(f.sin*Tf*pi+cos(2*f.sin*pi*(Tf-t0.sin))-cos(2*f.sin*pi*t0.sin)))
# calib.sinusoidal <- function(t, f, t0, a, b){return((a+b)/2 + (b-a)*sin(2*pi*f*(t-t0))/2)}
# true[,2] <- calib.sinusoidal(xx, f.sin, t0.sin, a.sin, b.sin)
# 
# ## Exponential calibration 
# t0.exp  <- Tf/2
# s.exp <- Tf/8
# b.exp <- 2
# a.exp <- (b.exp*s.exp*sqrt(2*pi)/Tf)*(pnorm(0,t0.exp,s.exp)-pnorm(Tf,t0.exp,s.exp))
# calib.exponential <- function(t, a, b, t0, s){return(a+b*exp(-(t-t0)^2/(2*s^2)))}
# true[,3] <- calib.exponential(xx, a.exp, b.exp, t0.exp, s.exp)
# 
# ## Base level
# eta0 <- 1
# 
# ## Additive calibration function
# calib.fct <- function(x1,x2,x3){
#   return(eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#          calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin)+
#          calib.exponential(x3, a.exp, b.exp, t0.exp, s.exp))}
# calib.fct12 <- function(x1,x2){eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#                                  calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin)}
# calib.fct13 <- function(x1,x2){eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#                                  calib.exponential(x2, a.exp, b.exp, t0.exp, s.exp)}
# calib.fct23 <- function(x1,x2){eta0+calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin) +
#                                  calib.exponential(x1, a.exp, b.exp, t0.exp, s.exp)}
# 
# ## Display calibration surface
# require(plot3D)
# xx2 <- seq(0,Tf,length.out = 1e2)
# Z12 <- outer(xx2,xx2,calib.fct12)
# Z13 <- outer(xx2,xx2,calib.fct13)
# Z23 <- outer(xx2,xx2,calib.fct23)
# 
# par(mfrow= c(1,3),  pty = "s", mar=c(1,1,4,1))
# persp3D(xx2,xx2,Z12, theta = 60, phi = 30, expanded = 2, ticktype = "simple", nticks = 4, colkey = list(plot = FALSE), xlab = "X1", ylab = "X2", zlab = "", main = "Calib. fct(X1,X2)", cex.main = 1)
# persp3D(xx2,xx2,Z13, theta = 60, phi = 30, expanded = 2, ticktype = "simple", nticks = 4, colkey = list(plot = FALSE), xlab = "X1", ylab = "X3", zlab = "", main = "Calib. fct(X1,X3)", cex.main = 1)
# persp3D(xx2,xx2,Z23, theta = 60, phi = 30, expanded = 2, ticktype = "simple", nticks = 4, colkey = list(plot = FALSE), xlab = "X2", ylab = "X3", zlab = "", main = "Calib. fct(X2,X3)", cex.main = 1)
# 
# ## A single dataset
# covariates.distr <- copula:::mvdc(copula:::normalCopula(rho, dim = 3), c("unif"), list(list(min = 0, max = Tf)), marginsIdentical = TRUE) 
# X <- copula:::rMvdc(n,covariates.distr)
# temp <- CondBiCopSim(fam, calib.fct, X, par2=par2, return.par=TRUE)
# calib <- temp$eta
# param <- temp$par
# pseudo <- temp$data
# dataset <- data.frame("u1" = pseudo[,1], "u2" = pseudo[,2], "x1" = X[,1], "x2" = X[,2], "x3" = X[,3])
# 
# ## Display the data
# dev.off()
# plot(pseudo[,1], pseudo[,2], xlab = "U1", ylab = "U2")
# 
# ## Display the copula parameter and calibration function
# par(mfrow=c(2,3))
# plot(X[,1],param, xlab = "X1", ylab = "Copula parameter")
# plot(X[,2],param, xlab = "X2", ylab = "Copula parameter")
# plot(X[,3],param, xlab = "X3", ylab = "Copula parameter")
# plot(X[,1],calib, xlab = "X1", ylab = "Calibration function")
# plot(X[,2],calib, xlab = "X2", ylab = "Calibration function")
# plot(X[,3],calib, xlab = "X3", ylab = "Calibration function")
# 
# ##
# ## 2) Estimate two models and display results
# ##
# 
# ## Model fit with a basis size (arguably) too small and unpenalized cubic spines
# pen <- FALSE
# basis <- c(3,3,3)
# formula <- ~s(x1, k=basis[1], bs = "cr", fx= !pen)+
#         s(x2, k=basis[2], bs = "cr", fx= !pen)+
#         s(x3, k=basis[3], bs = "cr", fx= !pen)
# system.time(fit <- gamBiCopEst(dataset, family = fam, parFrhs = formula))
# 
# ## Model fit with a better basis size and penalized cubic splines
# pen <- TRUE
# basis2 <- c(3,10,10)
# formula <- ~s(x1, k=basis2[1], bs = "cr", fx= !pen)+
#   s(x2, k=basis2[2], bs = "cr", fx= !pen)+
#   s(x3, k=basis2[3], bs = "cr", fx= !pen)
# system.time(fit2 <- gamBiCopEst(dataset, family = fam, parFrhs = formula))
# 
# ## Extract the gamBiCop object and use various methods
# res <- fit$res
# res2 <- fit2$res
# show(res)
# show(res2)
# logLik(res)
# logLik(res2)
# AIC(res)
# AIC(res2)
# BIC(res)
# BIC(res2)
# fitted <- gamBiCopPred(res, newdata = Xx, type="terms")$calib
# fitted2 <- gamBiCopPred(res2, newdata = Xx, type="terms")$calib
# 
# ## Spline approximation of each true smooth function for the two basis sizes
# for(i in 1:length(basis)){
#   temp <- rep(0,3)
#   temp[i] <- 1
#   temp <- data.frame("y" = calib.fct(xx*temp[1],xx*temp[2],xx*temp[3]), "x1" = xx, "x2" = 0, "x3" = 0)  
#   form <- ~s(x1, k=basis[i], bs = "cr", fx= !pen)
#   temp <- gam(update(form,y~.), data = temp)
#   true.approx[,i] <- predict.gam(temp, type = "terms")
#   
#   temp <- rep(0,3)
#   temp[i] <- 1
#   temp <- data.frame("y" = calib.fct(xx*temp[1],xx*temp[2],xx*temp[3]), "x1" = xx, "x2" = 0, "x3" = 0)  
#   form <- ~s(x1, k=basis2[i], bs = "cr", fx= !pen)
#   temp <- gam(update(form,y~.), data = temp)
#   true.approx2[,i] <- predict.gam(temp, type = "terms")
# }
# 
# ## Compute approximation and fitted biases
# bias.approx <- true.approx - true
# bias.approx2 <- true.approx2 - true
# bias.fitted2 <- fitted2 - true
# bias.fitted <- fitted - true
# 
# ## Display results
# yy <- range(true, true.approx, fitted, true.approx2, fitted2)
# yy[1] <- yy[1]*1.5
# yy.bias <- range(bias.approx, bias.fitted, bias.approx2, bias.fitted2)
# par(mfrow = c(2,3))
# plot(xx, true[,1], type = "l", ylim = yy, xlab = "Covariate 1", ylab = "Smooth 1")
# lines(xx, fitted[,1], col = "red")
# lines(xx, true.approx[,1], col = "red", lty = 2)
# lines(xx, fitted2[,1], col = "green")
# lines(xx, true.approx2[,1], col = "green", lty = 2)
# legend("bottomleft", cex = 0.6, c("True", "Fitted", "Appox 1", "Fitted 2", "Approx 2"), col = c("black", "red", "red", "green", "green"), lty = c(1,1,2,1,2))
# 
# plot(xx, true[ ,2], type = "l", ylim = yy, xlab = "Covariate 1", ylab = "Smooth 1")
# lines(xx, fitted[ ,2], col = "red")
# lines(xx, true.approx[ ,2], col = "red", lty = 2)
# lines(xx, fitted2[ ,2], col = "green")
# lines(xx, true.approx2[ ,2], col = "green", lty = 2)
# 
# plot(xx, true[ ,3], type = "l", ylim = yy, xlab = "Covariate 1", ylab = "Smooth 1")
# lines(xx, fitted[ ,3], col = "red")
# lines(xx, true.approx[ ,3], col = "red", lty = 2)
# lines(xx, fitted2[ ,3], col = "green")
# lines(xx, true.approx2[ ,3], col = "green", lty = 2)
# 
# plot(xx, bias.fitted[,1], col = "red", type = "l", ylim = yy.bias, xlab = "Covariate 1", ylab = "Bias 1")
# lines(xx, bias.approx[,1], col = "red", lty = 2)
# lines(xx, bias.fitted2[,1], col = "green", lty = 1)
# lines(xx, bias.approx2[,1], col = "green", lty = 2)
# legend("bottomleft", cex = 0.6, c("Fitted", "Appox 1", "Fitted 2", "Approx 2"), col = c("red", "red", "green", "green"), lty = c(1,2,1,2))
# 
# plot(xx, bias.fitted[,2], col = "red", type = "l", ylim = yy.bias, xlab = "Covariate 1", ylab = "Bias 1")
# lines(xx, bias.approx[,2], col = "red", lty = 2)
# lines(xx, bias.fitted2[,2], col = "green", lty = 1)
# lines(xx, bias.approx2[,2], col = "green", lty = 2)
# 
# plot(xx, bias.fitted[,3], col = "red", type = "l", ylim = yy.bias, xlab = "Covariate 1", ylab = "Bias 1")
# lines(xx, bias.approx[,3], col = "red", lty = 2)
# lines(xx, bias.fitted2[,3], col = "green", lty = 1)
# lines(xx, bias.approx2[,3], col = "green", lty = 2)
