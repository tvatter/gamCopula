# library(copula)
# library(VineCopula)
# library(numDeriv)
# library(microbenchmark)
# 
# FrankPar2Tau <- function(par) {
#   f <- function(x) x/(exp(x) - 1)
#   fu <- function(x) integrate(f, lower = 0, upper = x)$value
#   fl <- function(x) integrate(f, lower = x, upper = 0)$value
#   if (any(par > 0)) {
#     tau <- 1 - 4/par + 4/par^2 * sapply(par, fu)
#   }
#   else {
#     tau <- 1 - 4/par - 4/par^2 * sapply(par, fl)
#   }
#   return(tau)
# }
# 
# FrankPar2Tau.new <- function(par) {
#   tau <- 1 - 4/par + 4/par * debye1(par)
#   return(tau)
# }
# 
# FrankTau2Par <- function(tau) {
#   a <- 1
#   if (tau < 0) {
#     a <- -1
#     tau <- -tau
#   }
#   f <- function(x) {
#     x/(exp(x) - 1)
#   }
#   tauF <- function(x) 1 - 4/x + 4/x^2 * 
#     integrate(f, lower = 0 + .Machine$double.eps^0.5, upper = x)$value
#   v <- uniroot(function(x) tau - tauF(x), lower = 0, upper = 500, 
#                tol = .Machine$double.eps^0.5)$root
#   return(a * v)
# }
# 
# FrankTau2Par.new <- function(tau) {
#   a <- 1
#   if (tau < 0) {
#     a <- -1
#     tau <- -tau
#   }
#   v <- uniroot(function(x) tau - test2(x), 
#               lower = 0 + .Machine$double.eps^0.5, upper = 500, 
#               tol = .Machine$double.eps^0.5)$root
#   return(a*v)
# }
# 
# bmPar2Tau <- microbenchmark(FrankPar2Tau(0.5), FrankPar2Tau.new(0.5))
# bmTau2Par <- microbenchmark(FrankTau2Par(0.5), FrankTau2Par.new(0.5))
# 
# par(mfrow=c(1,2))
# boxplot(bmPar2Tau)
# boxplot(bmTau2Par)
# 
# FrankPar2Tau <- function(par) {
#   tau <- 1 - 4/par + 4/par * debye1(par)
#   return(tau)
# }
# 
# FrankTau2Par <- function(tau) {
#   mypar <- function(tau) {
#     a <- 1
#     if (tau < 0) {
#       a <- -1
#       tau <- -tau
#     }
#     a*safeUroot(function(x) tau - FrankPar2Tau(x), 
#               interval = c(0 + .Machine$double.eps^0.7, upper = 1e7))$root
#   }
#   return(sapply(tau,mypar))
# }
# 
# Frankdtaudpar <- function(par) {
#   tmp <- grad(FrankPar2Tau,par, method.args =
#                 list(eps=1e-4, d=0.0001, 
#                      zero.tol=sqrt(.Machine$double.eps/7e-7), 
#                      r=2, v=2, show.details=FALSE))
#   mm <- splinefun(par,tmp)
#   list(d1=tmp,d2=mm(par,deriv=1))
# }
# 
# Frankdpardtau <- function(tau) {
#   tmp <- grad(FrankTau2Par,tau, method.args =
#                 list(eps=1e-4, d=0.0001, 
#                      zero.tol=sqrt(.Machine$double.eps/7e-7), 
#                      r=2, v=2, show.details=FALSE))
#   mm <- splinefun(tau,tmp)
#   list(d1=tmp,d2=mm(tau,deriv=1))
# }
# 
# par <- seq(-4,4,1e-1)
# par <- par[-which(par == 0)]
# tau <- seq(-0.99,0.99,1e-2)
# tau <- tau[-which(tau == 0)]
# dp <- Frankdtaudpar(par)
# dt <- Frankdpardtau(tau)
# par(mfrow=c(1,3))
# plot(par,FrankPar2Tau(par), type = "l")
# lines(par,dp$d1, col = "red")
# lines(par,dp$d2, col = "green")
# plot(tau,FrankTau2Par(tau), type = "l")
# lines(tau,dt$d1, col = "red")
# lines(tau,dt$d2, col = "green")
# 
# bm <- microbenchmark(p2t=FrankPar2Tau(0.5), dtdp=Frankdtaudpar(0.5), 
#                             t2p=FrankTau2Par(0.5), dpdt=Frankdpardtau(0.5))
# boxplot(bm)