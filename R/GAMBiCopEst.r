#' Maximum penalized likelihood estimation of a Generalized Additive model 
#' for the copula parameter or Kendall's tau.
#' 
#' This function estimates the parameter(s) of a Generalized Additive model 
#' (gam) for the copula parameter or Kendall's tau.
#' It solves the maximum penalized likelihood estimation for the copula families 
#' supported in this package by reformulating each Newton-Raphson iteration as 
#' a generalized ridge regression, which is solved using 
#' the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param data A list or data frame containing the model responses, (u1,u2) in 
#' [0,1]x[0,1], and covariates required by the formula.
#' @param formula A gam formula (see \code{\link{gam}}, 
#' \code{\link{formula.gam}} and \code{\link{gam.models}} 
#' from \code{\link[mgcv:mgcv-package]{mgcv}}).
#' @param family A copula family: \code{1} Gaussian, \code{2} Student t, 
#' \code{3} Clayton, \code{4} Gumbel, \code{13} Survival Clayton, 
#' \code{14} Survival Gumbel,  \code{23} Rotated (90 degrees) Clayton, 
#' \code{24} Rotated (90 degrees) Gumbel, 
#' \code{33} Rotated (270 degrees) Clayton and 
#' \code{34} Rotated (270 degrees) Gumbel.
#' @param tau \code{FALSE} (default) for a calibration fonction specified for 
#' the Copula parameter or \code{TRUE} for a calibration function specified 
#' for Kendall's tau.
#' @param method \code{'NR'} for Newton-Raphson
#' and  \code{'FS'} for Fisher-scoring (default).
#' @param tol.rel Relative tolerance for \code{'FS'}/\code{'NR'} algorithm.
#' @param n.iters Maximal number of iterations for 
#' \code{'FS'}/\code{'NR'} algorithm.
#' @param verbose \code{TRUE} if informations should be printed during the 
#' estimation and \code{FALSE} (default) for a silent version.
#' @param ... Additional parameters to be passed to \code{\link{gam}} 
#' from \code{\link[mgcv:mgcv-package]{mgcv}}.
#' @return \code{gamBiCopEst} returns a list consisting of
#' \item{res}{S4 \code{\link{gamBiCop-class}} object.}
#' \item{method}{\code{'FS'} for Fisher-scoring and 
#' \code{'NR'} for Newton-Raphson.}
#' \item{tol.rel}{relative tolerance for \code{'FS'}/\code{'NR'} algorithm.}
#' \item{n.iters}{maximal number of iterations for 
#' \code{'FS'}/\code{'NR'} algorithm.}
#' \item{trace}{the estimation procedure's trace.}
#' \item{conv}{\code{0} if the algorithm converged and \code{1} otherwise.}
#' @seealso \code{\link{gamBiCop}} and \code{\link{gamBiCopSim}}.
#' @examples
#' set.seed(1)
#' 
#' ## Simulation parameters (sample size, correlation between covariates,
#' ## Clayton copula family)
#' n <- 2e2
#' rho <- 0.5
#' fam <- 3
#' 
#' 
#' ## A calibration surface depending on three variables
#' eta0 <- 1
#' calib.surf <- list(
#'   calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
#'     Tm <- (Tf - Ti)/2
#'     a <- -(b/3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
#'   return(a + b * (t - Tm)^2)},
#'   calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1) {
#'     a <- b * (1 - 2 * Tf * pi/(f * Tf * pi +
#'                                  cos(2 * f * pi * (Tf - Ti))
#'                                - cos(2 * f * pi * Ti)))
#'     return((a + b)/2 + (b - a) * sin(2 * f * pi * (t - Ti))/2)},
#'   calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf/8) {
#'     Tm <- (Tf - Ti)/2
#'     a <- (b * s * sqrt(2 * pi)/Tf) * (pnorm(0, Tm, s) - pnorm(Tf, Tm, s))
#'     return(a + b * exp(-(t - Tm)^2/(2 * s^2)))})
#' 
#' ## Display the calibration surface
#' par(mfrow = c(1, 3), pty = "s", mar = c(1, 1, 4, 1))
#' u <- seq(0, 1, length.out = 100)
#' sel <- matrix(c(1, 1, 2, 2, 3, 3), ncol = 2)
#' jet.colors <- colorRamp(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
#'                           "yellow", "#FF7F00", "red", "#7F0000"))
#' jet <- function(x) rgb(jet.colors(exp(x/3)/(1 + exp(x/3))), 
#'                        maxColorValue = 255)
#' for (k in 1:3) {
#'     tmp <- outer(u, u, function(x, y) 
#'       eta0 + calib.surf[[sel[k,1]]](x) + calib.surf[[sel[k, 2]]](y))
#'     persp(u, u, tmp, border = NA, theta = 60, phi = 30, zlab = "", 
#'           col = matrix(jet(tmp), nrow = 100), 
#'           xlab = paste("X", sel[k, 1], sep = ""), 
#'           ylab = paste("X", sel[k,2], sep = ""), 
#'           main = paste("eta0+f", sel[k, 1], 
#'                        "(X", sel[k, 1], ") +f",sel[k, 2], 
#'                        "(X", sel[k, 2], ")", sep = ""))
#' }
#' 
#' ## 3-dimensional matrix X of covariates
#' covariates.distr <- copula::mvdc(copula::normalCopula(rho, dim = 3),
#'                                  c("unif"), list(list(min = 0, max = 1)),
#'                                  marginsIdentical = TRUE)
#' X <- copula::rMvdc(n, covariates.distr)
#' 
#' ## U in [0,1]x[0,1] with copula parameter depending on X
#' U <- CondBiCopSim(fam, function(x1,x2,x3) {eta0+sum(mapply(function(f,x)
#'   f(x), calib.surf, c(x1,x2,x3)))}, X[,1:3], par2 = 6, return.par = TRUE)
#' 
#' ## Merge U and X
#' data <- data.frame(U$data,X)
#' names(data) <- c(paste("u",1:2,sep=""),paste("x",1:3,sep=""))
#' 
#' ## Display the data
#' dev.off()
#' plot(data[, "u1"], data[, "u2"], xlab = "U1", ylab = "U2")
#' 
#' ## Model fit with a basis size (arguably) too small 
#' ## and unpenalized cubic spines
#' pen <- FALSE
#' basis0 <- c(3, 4, 4)
#' formula <- ~s(x1, k = basis0[1], bs = "cr", fx = !pen) + 
#'   s(x2, k = basis0[2], bs = "cr", fx = !pen) + 
#'   s(x3, k = basis0[3], bs = "cr", fx = !pen)
#' system.time(fit0 <- gamBiCopEst(data, formula, fam, method = "FS"))
#' 
#' ## Model fit with a better basis size and penalized cubic splines (via min GCV)
#' pen <- TRUE
#' basis1 <- c(3, 10, 10)
#' formula <- ~s(x1, k = basis1[1], bs = "cr", fx = !pen) + 
#'   s(x2, k = basis1[2], bs = "cr", fx = !pen) + 
#'   s(x3, k = basis1[3], bs = "cr", fx = !pen)
#' system.time(fit1 <- gamBiCopEst(data, formula, fam, method = "FS"))
#' 
#' ## Extract the gamBiCop objects and show various methods
#' (res <- sapply(list(fit0,fit1), function(fit){fit$res}))
#' metds <- list('logLik'=logLik,'AIC'=AIC,'BIC'=BIC,'EDF'=EDF)
#' lapply(res, function(x) sapply(metds, function(f) f(x)))
#' 
#' 
#' ## Comparison between fitted, true smooth and spline approximation for each
#' ## true smooth function for the two basis sizes
#' fitted <- lapply(res, function(x) gamBiCopPred(x, data.frame(x1=u,x2=u,x3=u), 
#'                                                type = "terms")$calib)
#' true <- vector("list", 3)
#' for (i in 1:3) {
#'     y <- eta0+calib.surf[[i]](u)   
#'     true[[i]]$true <- y - eta0   
#'     temp <- mgcv::gam(y ~ s(u, k = basis0[i], bs = "cr", fx = T))
#'     true[[i]]$approx <- mgcv::predict.gam(temp, type = "terms")
#'     temp <- mgcv::gam(y ~s(u, k = basis1[i], bs = "cr", fx = F))
#'     true[[i]]$approx2 <- mgcv::predict.gam(temp, type = "terms")
#' }
#' 
#' ## Display results
#' par(mfrow = c(1, 3), pty = "s")
#' yy <- range(true, fitted)
#' yy[1] <- yy[1] * 1.5
#' for(k in 1:3){
#'   plot(u, true[[k]]$true, type = "l", ylim = yy, 
#'        xlab = paste("Covariate",k), ylab = paste("Smooth",k))
#'   lines(u, true[[k]]$approx, col = "red", lty = 2)
#'   lines(u, fitted[[1]][, k], col = "red")
#'   lines(u, fitted[[2]][, k], col = "green")
#'   lines(u, true[[k]]$approx2, col = "green", lty = 2)
#'   legend("bottomleft", cex = 0.6, lty = c(1, 1, 2, 1, 2),
#'          c("True", "Fitted", "Appox 1", "Fitted 2", "Approx 2"), 
#'          col = c("black", "red", "red", "green", "green"))
#' }
#' @export
gamBiCopEst <- function(data, formula = ~1, family = 1, 
                        tau = FALSE, method = "FS", 
                        tol.rel = 0.001, n.iters = 10, 
                        verbose = FALSE, ...) {

  if (!(is.list(data) || is.data.frame(data))) {
    stop("data has to be either a list or a data frame")
  } else {
    if(is.list(data)){
      if(!is.null(data$xt)){
        xt <- data$xt
        data <- as.data.frame(data[-which(names(data) == "xt")])
      }else{
        data <- as.data.frame(data)
      }     
    }
    n <- dim(data)[1]
    m <- dim(data)[2]
    u1 <- data$u1
    u2 <- data$u2
    u <- cbind(u1,u2)
  }
  
  if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
    stop("u1 and/or u2 are missing.")
  if (n < 2) 
    stop("Number of observations has to be at least 2.")
  if (any(u > 1) || any(u < 0)) 
    stop("(u1, u2) have to be in [0,1]x[0,1].")
  
  options(warn = -1)
  if (is.null(n.iters) || is.na(as.integer(n.iters)) || 
        (as.integer(n.iters) < 1) || 
        (as.integer(n.iters) !=  as.numeric(n.iters))) {
    stop("N.iters should be a positive integer.")
  } else {
    n.iters <- as.integer(n.iters)
  }
  
  if (is.null(tau) ||  is.na(tau) ||  
        !(is.logical(tau) || (tau == 0) || (tau == 1))) {
    stop("Tau should takes 0/1 or FALSE/TRUE to specify a model for 
         the copula parameter/Kendall's tau.")
  }
  
  if (is.null(tol.rel) ||  is.na(as.numeric(tol.rel)) || 
        (as.numeric(tol.rel) < 0) || (as.numeric(tol.rel) > 1)) {
    stop("Tol.rel should be a real number in [0,1].")
  } else {
    tol.rel <- as.numeric(tol.rel)
  }
  
  if (is.null(method) || !is.element(method, c("FS", "NR"))) {
    stop("Method should be a string, either NR (Newton-Raphson) 
         or FS (Fisher-scoring, faster but unstable).")
  }
  
  if (is.null(verbose) ||  is.na(verbose) || 
        !(is.logical(verbose) || (verbose == 0) || (verbose == 1))) {
    stop("Verbose should takes 0/1 or FALSE/TRUE.")
  } else {
    verbose <- as.logical(verbose)
  }
  options(warn = 0)

  temp <- VineCopula:::fasttau(u1, u2)
  rotated <- family
  
  if (!is.element(family, c(1, 2, 3, 4, 13, 14, 23, 24, 33, 34))) {
    stop("Copula family not yet implemented!")
  } else if (is.element(family, c(3, 4, 13, 14)) && (temp < 0)) {
    stop("This copula family cannot be used for negatively dependent data.")
  } else if (is.element(family, c(23, 24, 33, 34)) && (temp > 0)) {
    stop("This copula family cannot be used for positively dependent data.")
  } else if (is.element(family, c(13, 14))) {
    u <- rotate.data(u, 180)    
    if (family == 13) {
      family <- 3
    } else {
      family <- 4
    }
  } else if (is.element(family, c(23, 24))) {
    u <- rotate.data(u, 90)
    if (family == 23) {
      family <- 3
    } else {
      family <- 4
    }
  } else if (is.element(family, c(33, 34))) {
    u <- rotate.data(u, 270)
    if (family == 33) {
      family <- 3
    } else {
      family <- 4
    }
  } 

  init <- BiCopEst(u1, u2, family = family, method = "mle")
  new.pars <- list()
  new.pars$par <- rep(init$par, n)
  u <- cbind(u, new.pars$par)
  
  if (family == 2) {
    new.pars$par2 <- rep(init$par2, n)
    u <- cbind(u, new.pars$par2)
  }
  
  if (tau == FALSE) {
    new.pars$partrans <- sapply(new.pars$par, function(x) 
      BiCopPar2Eta(family, x))
  } else {
    if (family == 2) {
      new.pars$tau <- sapply(new.pars$par, function(x) 
        BiCopPar2Tau(family, x, init$par2))
    } else {
      new.pars$tau <- sapply(new.pars$par, function(x) 
        BiCopPar2Tau(family, x))
    }
    new.pars$partrans <- 2 * atanh(new.pars$tau)
  }
  
  old.pars <- new.pars
  
  if (verbose == 1) {
    t <- Sys.time()
    print(paste("gam iteration", 1))
  }

  temp <- derivatives.par(u, new.pars, family, method, tau)
  temp <- as.data.frame(wz.update(temp, new.pars, family, method, tau))
  temp <- cbind(temp, data)
  par.formula <- update(formula, z ~ .)
  
  w <- NULL
  res <- tryCatch({
    mm <- gam(par.formula, data = temp, weights = w, ...)
  }, error = function(err) {
    print(paste("A problem occured at the first iteration of the ", 
                method, "algorithm. The ERROR comming from", 
                "mgcv's gam function is:"))
    stop(err)
    print("......The results should not be trusted!")
    return(NULL)
  })
  if (verbose == 1) {
    print(Sys.time() - t)
  }
  stopifnot(!is.null(res))
  
  temp <- pars.update(mm, family, temp, tau)
  
  new.pars$par <- u[, 3] <- temp$par
  new.pars$partrans <- temp$partrans
  if (tau) {
    new.pars$tau <- temp$tau
  }
  
  if (family == 2) {
    if (verbose == 1) {
      print(paste("DF iteration", 1))
      t <- Sys.time()
    }
    
    LL <- function(nu) {
      if (2 + 1e-08 + exp(nu) == Inf) {
        nu <- log(30)
      }
      -sum(log(apply(data, 1, function(x) BiCopPDF(x[1], x[2], family = 2, 
        x[3], 2 + 1e-08 + exp(nu)))), na.rm = TRUE)
    }
    
    nu <- optimize(LL, c(log(2), log(30)))
    u[, 4] <- new.pars$par2 <- rep(2 + 1e-08 + exp(nu$minimum), n)
    
    if (verbose == 1) {
      print(u[1, 4])
      print(Sys.time() - t)
    }
  }
  
  trace <- numeric(n.iters)
  tt <- trace.update(old.pars$partrans, new.pars$partrans)
  trace[1] <- tt$trace

  eps <- tt$eps
  k <- 1
  conv <- 0
  while ((k < n.iters) & (eps > tol.rel)) {
    
    k <- k + 1
    if(k == n.iters){
      conv <- 1
    }
    old.pars <- new.pars
    
    if (verbose == 1) {
      t <- Sys.time()
      print(paste("gam iteration", k))
    }
    
    temp <- derivatives.par(u, new.pars, family, method, tau)
    temp <- as.data.frame(wz.update(temp, new.pars, family, method, tau))
    temp <- cbind(temp, data)
    
    res <- tryCatch({
      mm <- gam(par.formula, data = temp, weights = w, 
                control = gam.control(keepData = TRUE), ...)
    }, error = function(err) {
      print(paste("A problem occured at the ", k, "th iteration of the ", 
                  method, "algorithm. The ERROR comming from", 
                  "mgcv's gam function is:"))
      print(err)
      print("......The results should not be trusted!")
      return(NULL)
    })
    if (is.null(res)) {
      conv <- 1
      break
    }
    if (verbose == 1) {
      print(Sys.time() - t)
    }
    
    temp2 <- pars.update(mm, family, temp, tau)
    new.pars$par <- u[, 3] <- temp2$par
    new.pars$partrans <- temp2$partrans
    if (tau) {
      new.pars$tau <- temp2$tau
    }
    
    if ((family == 2)) {
      # && (k %% 2 == 0)){
      if (verbose == 1) {
        print(paste("DF iteration", k))
        t <- Sys.time()
      }
      
      LL <- function(nu) {
        if (2 + 1e-08 + exp(nu) == Inf) {
          nu <- log(30)
        }
        -sum(log(apply(data, 1, function(x) BiCopPDF(x[1], x[2], family = 2, 
          x[3], 2 + 1e-08 + exp(nu)))), na.rm = TRUE)
      }
      
      nu <- optimize(LL, c(log(2), log(30)))
      u[, 4] <- new.pars$par2 <- rep(2 + 1e-08 + exp(nu$minimum), n)
      
      if (verbose == 1) {
        print(u[1, 4])
        print(Sys.time() - t)
      }
    }
    tt <- trace.update(old.pars$partrans, new.pars$partrans)
    if (is.na(tt$eps)) {
      print(paste("A problem occured at the ", k, "th iteration of the ", 
                  method, "algorithm... The results should not be trusted!"))
      conv <- 1
      break
    } else {
      trace[k] <- abs(tt$trace)
      eps <- tt$eps
    }
  }
  
  if (family == 2) {
    res <- gamBiCop(rotated, mm, u[1, 4], tau)
  } else {
    res <- gamBiCop(rotated, mm, 0, tau)
  }
  out <- list(res = res, method = method, tol.rel = tol.rel, 
              n.iters = n.iters, trace = trace[1:k], conv = conv)
  return(out)
} 
