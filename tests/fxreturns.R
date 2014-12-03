# # see prepare_data.R
# rm(list=ls())
# library(gamVineCopula)
# library(xts)
# library(matlab)
# load("data/FXreturns.Rda")
# 
# n <- dim(ret)[1]
# ret <- ret[(n-1e4):n,]
# mu <- colMeans(ret)
# r <- sapply(ret, function(x) x - mean(x))
# 
# t <- index(ret)
# dt <- diff(as.numeric(t[1:2]))
# n <- length(t)
# m <- (24*3600)/dt
# 
# 
# #pdf(file = paste(figpath, "highfreq-data1.pdf", sep = ""), width = 21, height = 14)
# par(mfrow = c(2,2), mar = c(2,5,3,1))
# dom <- .indexyday(as.xts(t))
# xx <- c(1,which(diff(dom) < 0)+1)
# nn <- length(xx)
# for(i in 1:4){
#   plot(1:n, as.numeric(ret[,i]), type = "l", 
#        ylim = range(ret), xlab = "", ylab = "log-return", main = names(ret)[i],
#        xaxt = "n", cex.lab = 2, cex.axis = 2, cex.main = 2)
#   axis(1,cex.axis=2, at = (1:n)[xx], labels = format(t[xx], "%Y"))
# }
# #dev.off()
# 
# # Remove daily volatility component using noise robust estimator
# sigma <- apply(r,2,function(x) get.medRV(x, m))
# #z <- r[(m/2):(m/2+dim(sigma)[1]-1)]/sqrt(sigma/(1/m))
# z <- r[m:n]/sqrt(sigma/(1/m))
# 
# 
# # Remove periodic volatility component using Fourier Flexible Form (reestimated weekly)
# K <- 6
# sr <- apply(r,2,function(x) get.FFF(x,K, m))
# sz <- apply(z,2,function(x) get.FFF(x,K, m))
# rs <- r/sr
# zs <- z/sz
# 
# a <- lapply(list(r,z,rs,zs), function(y) apply(y,2, function(x) 
#     acf(abs(as.numeric(x)), m*5, plot = FALSE)$acf[2:(m*5+1)]))
# 
# plot.a <- function(a,k,l,m){
#   apply(a[[k]],2,function(x) 
#     plot(1:l/m,x, type = "l", 
#          ylim = range(x), xlab = "lags", ylab = "acf",
#          cex.lab = 2, cex.axis = 2))
# }
# par(mfrow=c(2,2))
# plot.a(a,1,m*5,m)
# 
# par(mfrow=c(2,2))
# plot.a(a,2,m*5,m)
# 
# par(mfrow=c(2,2))
# plot.a(a,3,m*5,m)
# 
# par(mfrow=c(2,2))
# plot.a(a,4,m*5,m)