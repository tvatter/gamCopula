\name{gamCopula-package}
\alias{gamCopula-package}
\alias{gamCopula}
\alias{gamCopula.package}
\docType{package}
\title{
Generalized additive models for bivariate dependence structures and vine copulas
}
\description{
This package implements various inference and simulation tools to apply generalized additive models to bivariate dependence structures and vine copulas.

More references in Vatter and Chavez-Demoulin (2014) and Vatter (2014).
}
\details{
\tabular{ll}{
Package: \tab gamCopula\cr
Type: \tab Package\cr
Version: \tab 0.0-1\cr
Date: \tab 2014-08-05\cr
License: \tab GPL-3\cr
}
}
\author{
Thibault Vatter

Maintainer: Thibault Vatter <thibault.vatter@gmail.com>
}
\references{
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009)
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.

Brechmann, E. C., C. Czado, and K. Aas (2012)
Truncated regular vines in high dimensions with applications to financial data.
Canadian Journal of Statistics 40 (1), 68-85.

Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013)
Selecting and estimating regular vine copulae and application to financial returns.
Computational Statistics & Data Analysis, 59 (1), 52-69.

Joe, H. (1997)
Multivariate Models and Dependence Concepts.
London: Chapman and Hall.

Schepsmeier, U. and J. Stoeber (2012).
Derivatives and Fisher information of bivariate copulas.
\url{http://link.springer.com/article/10.1007/s00362-013-0498-x}.

Vatter, T. and V. Chavez-Demoulin (2014)
Generalized additive models for conditional dependence structures.
\url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2265534}

Vatter, T. (2014)
Generalized additive models for non-simplified pair-copula constructions.

Wood, S.N. (2004)
Stable and efficient multiple smoothing parameter estimation for generalized additive models. 
J. Amer. Statist. Ass. 99:673-686. 

Wood, S.N. (2006) 
eneralized Additive Models: an introduction with R.
CRC.

}
\keyword{ package }
\seealso{
The present package is heavly relying on the \code{\link[mgcv:mgcv-package]{mgcv}} and \code{\link{VineCopula}} packages, as it basically extends and mix both of them.
}
\examples{
set.seed(1)

## Simulation parameters (sample size, correlation between covariates,
## Clayton copula family)
n <- 2e2
rho <- 0.5
fam <- 3


## A calibration surface depending on three variables
eta0 <- 1
calib.surf <- list(
  calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
    Tm <- (Tf - Ti)/2
    a <- -(b/3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
  return(a + b * (t - Tm)^2)},
  calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1) {
    a <- b * (1 - 2 * Tf * pi/(f * Tf * pi +
                                 cos(2 * f * pi * (Tf - Ti))
                               - cos(2 * f * pi * Ti)))
    return((a + b)/2 + (b - a) * sin(2 * f * pi * (t - Ti))/2)},
  calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf/8) {
    Tm <- (Tf - Ti)/2
    a <- (b * s * sqrt(2 * pi)/Tf) * (pnorm(0, Tm, s) - pnorm(Tf, Tm, s))
    return(a + b * exp(-(t - Tm)^2/(2 * s^2)))})

## Display the calibration surface
par(mfrow = c(1, 3), pty = "s", mar = c(1, 1, 4, 1))
u <- seq(0, 1, length.out = 100)
sel <- matrix(c(1, 1, 2, 2, 3, 3), ncol = 2)
jet.colors <- colorRamp(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                          "yellow", "#FF7F00", "red", "#7F0000"))
jet <- function(x) rgb(jet.colors(exp(x/3)/(1 + exp(x/3))),
                       maxColorValue = 255)
for (k in 1:3) {
    tmp <- outer(u, u, function(x, y)
      eta0 + calib.surf[[sel[k,1]]](x) + calib.surf[[sel[k, 2]]](y))
    persp(u, u, tmp, border = NA, theta = 60, phi = 30, zlab = "",
          col = matrix(jet(tmp), nrow = 100),
          xlab = paste("X", sel[k, 1], sep = ""),
          ylab = paste("X", sel[k,2], sep = ""),
          main = paste("eta0+f", sel[k, 1],
                       "(X", sel[k, 1], ") +f",sel[k, 2],
                       "(X", sel[k, 2], ")", sep = ""))
}

## 3-dimensional matrix X of covariates
covariates.distr <- copula::mvdc(copula::normalCopula(rho, dim = 3),
                                 c("unif"), list(list(min = 0, max = 1)),
                                 marginsIdentical = TRUE)
X <- copula::rMvdc(n, covariates.distr)

## U in [0,1]x[0,1] with copula parameter depending on X
U <- CondBiCopSim(fam, function(x1,x2,x3) {eta0+sum(mapply(function(f,x)
  f(x), calib.surf, c(x1,x2,x3)))}, X[,1:3], par2 = 6, return.par = TRUE)

## Merge U and X
data <- data.frame(U$data,X)
names(data) <- c(paste("u",1:2,sep=""),paste("x",1:3,sep=""))

## Display the data
dev.off()
plot(data[, "u1"], data[, "u2"], xlab = "U1", ylab = "U2")

## Model fit with a basis size (arguably) too small
## and unpenalized cubic spines
pen <- FALSE
basis0 <- c(3, 4, 4)
formula <- ~s(x1, k = basis0[1], bs = "cr", fx = !pen) +
  s(x2, k = basis0[2], bs = "cr", fx = !pen) +
  s(x3, k = basis0[3], bs = "cr", fx = !pen)
system.time(fit0 <- gamBiCopEst(data, formula, fam, method = "FS"))

## Model fit with a better basis size and penalized cubic splines (via min GCV)
pen <- TRUE
basis1 <- c(3, 10, 10)
formula <- ~s(x1, k = basis1[1], bs = "cr", fx = !pen) +
  s(x2, k = basis1[2], bs = "cr", fx = !pen) +
  s(x3, k = basis1[3], bs = "cr", fx = !pen)
system.time(fit1 <- gamBiCopEst(data, formula, fam, method = "FS"))

## Extract the gamBiCop objects and show various methods
(res <- sapply(list(fit0,fit1), function(fit){fit$res}))
metds <- list('logLik'=logLik,'AIC'=AIC,'BIC'=BIC,'EDF'=EDF)
lapply(res, function(x) sapply(metds, function(f) f(x)))


## Comparison between fitted, true smooth and spline approximation for each
## true smooth function for the two basis sizes
fitted <- lapply(res, function(x) gamBiCopPred(x, data.frame(x1=u,x2=u,x3=u),
                                               type = "terms")$calib)
true <- vector("list", 3)
for (i in 1:3) {
    y <- eta0+calib.surf[[i]](u)
    true[[i]]$true <- y - eta0
    temp <- mgcv::gam(y ~ s(u, k = basis0[i], bs = "cr", fx = T))
    true[[i]]$approx <- mgcv::predict.gam(temp, type = "terms")
    temp <- mgcv::gam(y ~s(u, k = basis1[i], bs = "cr", fx = F))
    true[[i]]$approx2 <- mgcv::predict.gam(temp, type = "terms")
}

## Display results
par(mfrow = c(1, 3), pty = "s")
yy <- range(true, fitted)
yy[1] <- yy[1] * 1.5
for(k in 1:3){
  plot(u, true[[k]]$true, type = "l", ylim = yy,
       xlab = paste("Covariate",k), ylab = paste("Smooth",k))
  lines(u, true[[k]]$approx, col = "red", lty = 2)
  lines(u, fitted[[1]][, k], col = "red")
  lines(u, fitted[[2]][, k], col = "green")
  lines(u, true[[k]]$approx2, col = "green", lty = 2)
  legend("bottomleft", cex = 0.6, lty = c(1, 1, 2, 1, 2),
         c("True", "Fitted", "Appox 1", "Fitted 2", "Approx 2"),
         col = c("black", "red", "red", "green", "green"))
}
}