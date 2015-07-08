# library(VineCopula)
# library(copula)
# 
# param <- 3
# cops <- list(cl00 = claytonCopula(param),
#              cl09 = r90ClaytonCopula(-param),
#              cl18 = surClaytonCopula(param),
#              cl27 = r270ClaytonCopula(-param))
# 
# myplot <- function(cop, nn, rr = c(-3,3)) {
#   mv <- mvdc(cop, c("norm", "norm"),
#              list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))
#   contour(mv, dMvdc, xlim = rr, ylim = rr, main = nn)
# }
# 
# par(mfrow=c(2,2))
# sapply(1:length(cops), function(x) myplot(cops[[x]], names(cops)[x]))
# 
# 
# fams <- c(3,23,13,33)
# 
# myplot2 <- function(fam, par, nn, rr = c(-3,3), dr = 1e-1) {
#   x <- seq(rr[1],rr[2],dr)
#   xx <- outer(x,x)
#   f <- outer(x, x, function(x, y) dnorm(x)*dnorm(y)*
#                BiCopPDF(pnorm(x),pnorm(y), fam, par))
#   contour(x,x,f)
# }
# par(mfrow=c(2,2), pts = "y")
# sapply(1:length(cops), function(x) myplot2(fams[x], cops[[x]]@parameters, 
#                                           names(cops)[x]))