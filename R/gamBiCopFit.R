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
#' @param data A list, data frame or matrix containing the model responses,
#'  (u1,u2) in [0,1]x[0,1], and covariates required by the formula.
#' @param formula A gam formula (see \code{\link{gam}}, 
#' \code{\link{formula.gam}} and \code{\link{gam.models}} 
#' from \code{\link[mgcv:mgcv-package]{mgcv}}).
#' @param family A copula family: \code{1} Gaussian, 
#' \code{2} Student t, 
#' \code{5} Frank, 
#' \code{301} Double Clayton type I (standard and rotated 90 degrees), 
#' \code{302} Double Clayton type II (standard and rotated 270 degrees), 
#' \code{303} Double Clayton type III (survival and rotated 90 degrees), 
#' \code{304} Double Clayton type IV (survival and rotated 270 degrees), 
#' \code{401} Double Gumbel type I (standard and rotated 90 degrees), 
#' \code{402} Double Gumbel type II (standard and rotated 270 degrees), 
#' \code{403} Double Gumbel type III (survival and rotated 90 degrees), 
#' \code{404} Double Gumbel type IV (survival and rotated 270 degrees).
#' @param tau \code{FALSE} (default) for a calibration function specified for 
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
#' @return \code{gamBiCopFit} returns a list consisting of
#' \item{res}{S4 \code{\link{gamBiCop-class}} object.}
#' \item{method}{\code{'FS'} for Fisher-scoring (default) and 
#' \code{'NR'} for Newton-Raphson.}
#' \item{tol.rel}{relative tolerance for \code{'FS'}/\code{'NR'} algorithm.}
#' \item{n.iters}{maximal number of iterations for 
#' \code{'FS'}/\code{'NR'} algorithm.}
#' \item{trace}{the estimation procedure's trace.}
#' \item{conv}{\code{0} if the algorithm converged and \code{1} otherwise.}
#' @seealso \code{\link{gamBiCop}} and \code{\link{gamBiCopSimulate}}.
#' @examples
#' require(copula)
#' require(mgcv)
#' set.seed(0)
#' 
#' ## Simulation parameters (sample size, correlation between covariates,
#' ## Gaussian copula family)
#' n <- 5e2
#' rho <- 0.5
#' fam <- 1
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
#' covariates.distr <- mvdc(normalCopula(rho, dim = 3),
#'                                  c("unif"), list(list(min = 0, max = 1)),
#'                                  marginsIdentical = TRUE)
#' X <- rMvdc(n, covariates.distr)
#' 
#' ## U in [0,1]x[0,1] with copula parameter depending on X
#' U <- condBiCopSim(fam, function(x1,x2,x3) {eta0+sum(mapply(function(f,x)
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
#' system.time(fit0 <- gamBiCopFit(data, formula, fam))
#' 
#' ## Model fit with a better basis size and penalized cubic splines (via min GCV)
#' pen <- TRUE
#' basis1 <- c(3, 10, 10)
#' formula <- ~s(x1, k = basis1[1], bs = "cr", fx = !pen) + 
#'   s(x2, k = basis1[2], bs = "cr", fx = !pen) + 
#'   s(x3, k = basis1[3], bs = "cr", fx = !pen)
#' system.time(fit1 <- gamBiCopFit(data, formula, fam))
#' 
#' ## Extract the gamBiCop objects and show various methods
#' (res <- sapply(list(fit0,fit1), function(fit){fit$res}))
#' metds <- list('logLik'=logLik,'AIC'=AIC,'BIC'=BIC,'EDF'=EDF)
#' lapply(res, function(x) sapply(metds, function(f) f(x)))
#' 
#' 
#' ## Comparison between fitted, true smooth and spline approximation for each
#' ## true smooth function for the two basis sizes
#' fitted <- lapply(res, function(x) gamBiCopPredict(x, data.frame(x1=u,x2=u,x3=u), 
#'                                                type = "terms")$calib)
#' true <- vector("list", 3)
#' for (i in 1:3) {
#'     y <- eta0+calib.surf[[i]](u)   
#'     true[[i]]$true <- y - eta0   
#'     temp <- gam(y ~ s(u, k = basis0[i], bs = "cr", fx = TRUE))
#'     true[[i]]$approx <- predict.gam(temp, type = "terms")
#'     temp <- gam(y ~s(u, k = basis1[i], bs = "cr", fx = FALSE))
#'     true[[i]]$approx2 <- predict.gam(temp, type = "terms")
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
gamBiCopFit <- function(data, formula = ~1, family = 1, tau = TRUE, 
                        method = "FS", tol.rel = 0.001, n.iters = 10, 
                        verbose = FALSE, ...) {
  
  tmp <- valid.gamBiCopFit(data, n.iters, tau, tol.rel, method, verbose, family)
  if (tmp != TRUE)
    stop(tmp)
  
  if (family == 5 && method == "FS") {
    # We should put a warning here, but this is quite annoying...
    #warning("'method' switched to 'NR' for the Frank copula.")
    method <- "NR"
  }
  
  if (is.list(data)){
    if(!is.null(data$xt)){
      xt <- data$xt
      data <- data[-which(names(data) == "xt")]
    }
    data <- as.data.frame(data)
  }
  
  n <- dim(data)[1]
  m <- dim(data)[2]
  if (!all(is.element(c("u1", "u2"), names(data)))) {
    names(data)[1:2] <- c("u1", "u2")
  }

  u1 <- data$u1
  u2 <- data$u2
  u <- cbind(u1,u2)
  
  if (verbose == 1) {
    if (family %in% c(1,2)) {
      cat(paste("Initialization with the MLE\n"))
    } else {
      cat(paste("Initialization with the empirical Kendall's tau\n"))
    }
    t <- Sys.time()
  }
  
  if (family %in% c(1,2)) {
    init <- BiCopEst(u1, u2, family, method = "mle")
  } else {
    init <- list()
    init$par <- tau2par(fasttau(u1, u2),family)
  }
  
  if (verbose == 1) {
    cat(paste(Sys.time() - t,"\n"))
  }
  
  new.pars <- list()
  new.pars$par <- rep(init$par, n)
  u <- cbind(u, new.pars$par)
  
  if (family == 2) {
    new.pars$par2 <- rep(init$par2, n)
    u <- cbind(u, new.pars$par2)
  } else {
    new.pars$par2 <- rep(0, n)
  }
  
  if (tau == FALSE) {
    new.pars$partrans <- par2eta(new.pars$par,family)
  } else {
    new.pars$tau <- par2tau(new.pars$par,family)
    new.pars$partrans <- 2 * atanh(new.pars$tau)
  }
  
  old.pars <- new.pars
  
  if (verbose == 1) {
    t <- Sys.time()
    cat(paste("gam iteration", 1, "\n"))
  }
  
  par.formula <- update(formula, z ~ .)
  
  w <- NULL
  res <- tryCatch({
    tmp <- derivatives.par(u, new.pars, family, method, tau)
    tmp <- as.data.frame(wz.update(tmp, new.pars, family, method, tau))
    tmp <- cbind(tmp, data)
    mm <- gam(par.formula, data = tmp, weights = w, 
              control = gam.control(keepData = TRUE), ...)
  }, error = function(err) {
    msg <- paste("A problem occured at the first iteration of the ", 
                 method, "algorithm:\n")
    msg <- c(msg, paste(err))
#     msg <- c(msg, paste("...... switching to", 
#                         switch(method, NR = "FS", FS = "NR"), "instead!\n"))
    message(msg)
    return(NULL)
  })
  
  
#   if (is.null(res)) {
#     method <- switch(method, NR = "FS", FS = "NR")
#     n.iters <- 2*n.iters
#     res <- tryCatch({
#       tmp <- derivatives.par(u, new.pars, family, method, tau)
#       tmp <- as.data.frame(wz.update(tmp, new.pars, family, method, tau))
#       tmp <- cbind(tmp, data)
#       mm <- gam(par.formula, data = tmp, weights = w, 
#                 control = gam.control(keepData = TRUE), ...)
#     }, error = function(err) {
#       msg <- paste("A problem occured at the first iteration of the ", 
#                    method, "algorithm:\n")
#       msg <- c(msg, paste(err))
#       msg <- c(msg, paste("...... The results should not be trusted!\n"))
#       message(msg)
#       return(NULL)
#     })
#   }
  if (verbose == 1) {
    print(Sys.time() - t)
  }
  stopifnot(!is.null(res))
  
  tmp <- pars.update(mm, family, tmp, tau)
  
  new.pars$par <- u[, 3] <- tmp$par
  new.pars$partrans <- tmp$partrans
  if (tau) {
    new.pars$tau <- tmp$tau
  }
  
  if (family == 2) {
    if (verbose == 1) {
      cat(paste("DF iteration", 1, "\n"))
      t <- Sys.time()
    }
    link <- function(nu) {
      2 + 1e-08 + exp(nu)
    }
    nllvec <- function(nu, u) {
      sel <- link(nu) == Inf
      if (any(sel)) {
        nu[sel] <- log(30)
      }
      nu <- link(nu)
      nll <- function(nu) -sum(log(bicoppd1d2(cbind(u[,1:3],nu),2)))
      return(sapply(nu, nll))
    }
    nu <- optimize(nllvec, c(log(2), log(30)), u)
    u[, 4] <- new.pars$par2 <- rep(link(nu$minimum), n)
    
    if (verbose == 1) {
      #print(u[1, 4])
      print(Sys.time() - t)
    }
  }
  
  trace <- numeric(n.iters)
  tt <- trace.update(old.pars$partrans, new.pars$partrans)
  trace[1] <- tt$trace
  
  eps <- tt$eps
  k <- 1
  conv <- 0
  while ((k < n.iters) & 
         (eps > tol.rel) & 
         (k == 1 || 
          abs(diff(trace[(k-1):k])) > min(.Machine$double.eps^.75,tol.rel))) {
    
    k <- k + 1
    if(k == n.iters){
      conv <- 1
    }
    old.pars <- new.pars
    
    if (verbose == 1) {
      t <- Sys.time()
      cat(paste("gam iteration", k, "\n"))
    }
    
    res <- tryCatch({
      tmp <- derivatives.par(u, new.pars, family, method, tau)
      tmp <- as.data.frame(wz.update(tmp, new.pars, family, method, tau))
      tmp <- cbind(tmp, data)
      mm <- gam(par.formula, data = tmp, weights = w, 
                control = gam.control(keepData = TRUE), ...)
    }, error = function(err) {
      msg <- paste("A problem occured at iteration", k, "of the ", 
                   method, "algorithm:\n")
      msg <- c(msg, paste(err))
#       msg <- c(msg, paste("...... switching to", 
#                           switch(method, NR = "FS", FS = "NR"), "instead!\n"))
      message(msg)
      return(NULL)
    })

#     if (is.null(res)) {
#       method <- switch(method, NR = "FS", FS = "NR")
#       n.iters <- 2*n.iters
#       res <- tryCatch({
#         tmp <- derivatives.par(u, new.pars, family, method, tau)
#         tmp <- as.data.frame(wz.update(tmp, new.pars, family, method, tau))
#         tmp <- cbind(tmp, data)
#         print("before gam")
#         mm <- gam(par.formula, data = tmp, weights = w, 
#                   control = gam.control(keepData = TRUE, trace = TRUE), ...)
#         print("after gam")
#       }, error = function(err) {
#         msg <- paste("A problem also occured at iteration", k, "of the ", 
#                      method, "algorithm:\n")
#         msg <- c(msg, paste(err))
#         msg <- c(msg, paste("...... The results should not be trusted!\n"))
#         message(msg)
#         return(NULL)
#       })
#     }

    
    if (is.null(res)) {
      conv <- 1
      break
    }
    if (verbose == 1) {
      print(Sys.time() - t)
    }
    
    tmp2 <- pars.update(mm, family, tmp, tau)
    new.pars$par <- u[, 3] <- tmp2$par
    new.pars$partrans <- tmp2$partrans
    if (tau) {
      new.pars$tau <- tmp2$tau
    }
    
    #     if ((family == 2)) {
    #       # && (k %% 2 == 0)){
    #       if (verbose == 1) {
    #         print(paste("DF iteration", k))
    #         t <- Sys.time()
    #       }
    #       
    #       nu <- optimize(nllvec, c(log(2), log(30)))
    #       u[, 4] <- new.pars$par2 <- rep(link(nu$minimum), n)
    #       
    #       if (verbose == 1) {
    #         print(u[1, 4])
    #         print(Sys.time() - t)
    #       }
    #     }
    tt <- trace.update(old.pars$partrans, new.pars$partrans)
    if (is.na(tt$eps)) {
      message(paste("A problem occured at the ", k, "th iteration of the ", 
                    method, "algorithm... The results should not be trusted!\n"))
      conv <- 1
      break
    } else {
      trace[k] <- abs(tt$trace)
      eps <- tt$eps
    }
  }
  
  if (family == 2) {
    if (verbose == 1) {
      cat("DF final iteration\n")
      t <- Sys.time()
    }
    nu <- optimize(nllvec, c(log(2), log(30)), u)
    
    if (verbose == 1) {
      print(Sys.time() - t)
    }
    u[, 4] <- new.pars$par2 <- rep(link(nu$minimum), n)
    par2 <- u[1, 4]
  } else {
    par2 <- 0
  }
  
  out <- list(res = gamBiCop(family, mm, par2, tau), 
              method = method, tol.rel = tol.rel, 
              n.iters = n.iters, trace = trace[1:k], conv = conv)
  return(out)
} 

valid.gamBiCopFit <- function(data, n.iters, tau, tol.rel, method, verbose, 
                              family) {
  if (!(is.list(data) || is.data.frame(data) || is.matrix(data))) {
    return("data has to be either a list, a data frame or a matrix.")
  } 
  
  if (is.list(data)){
    if(!is.null(data$xt)){
      data <- data[-which(names(data) == "xt")]
    }
    data <- tryCatch(as.data.frame(data), error = function(e) e$message)
    if (is.character(data)) {
      return(data)
    }
  } else {
    data <- as.data.frame(data)
  }

  if (is.null(names(data)) || !all(is.element(c("u1", "u2"), names(data)))) {
    warnings(paste0("u1 and u2 not found in data. The first two columns are",
                    " assumed to be the response variables"))
    names(data)[1:2] <- c("u1", "u2")
  }
  n <- dim(data)[1]
  m <- dim(data)[2]
  u1 <- data$u1
  u2 <- data$u2
  u <- cbind(u1,u2)
  
  if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
    return("u1 and/or u2 are missing.")
  if (n < 2) 
    return("Number of observations has to be at least 2.")
  if (any(u > 1) || any(u < 0)) 
    return("(u1, u2) have to be in [0,1]x[0,1].")
  
  options(warn = -1)
  if (!valid.posint(n.iters)) {
    return(msg.posint(var2char(n.iters)))
  }
  
  if (!valid.logical(tau)) {
    return(msg.logical(var2char(tau)))
  }
  
  if (!valid.unif(tol.rel)) {
    return(msg.unif(var2char(tol.rel)))
  } 
  
  if (is.null(method) || length(method) != 1 || 
      !is.element(method, c("FS", "NR"))) {
    return("Method should be a string, either NR (Newton-Raphson) 
           or FS (Fisher-scoring, faster but unstable).")
  }
  
  if (!valid.logical(verbose)) {
    return(msg.logical(var2char(verbose)))
  } 
  
  #tmp <- fasttau(u1, u2)
  
  if (!valid.family(family)) {
    return(msg.family(var2char(family)))
  }  
  
  #   if (tmp > 0) {
  #     if (!valid.familypos(family, tmp)) {
  #       return(msg.familypos(var2char(family)))
  #     }
  #   } else {
  #     if (!valid.familyneg(family, tmp)) {
  #       return(msg.familyneg(var2char(family)))
  #     }
  #   } 
  
  options(warn = 0)
  
  return(TRUE)
}

## Update parameters with respect to fited gam model
"pars.update" <- function(mm, family, data, tau) {
  
  out <- list()
  fitted <- predict(mm, data)
  
  out$partrans <- as.numeric(fitted)
  
  # Define links between Kendall's tau, copula parameter and calibration 
  # function... the cst/cstinv make sure that the boundaries are never attained
  if (family %in% c(1, 2)) {
    cstpar <- csttau <- function(x) x*(1-1e-5) 
  } else {
    csttau <- function(x) x*(1-1e-5)
    cstpar <- function(x) {
      sign(x)*pmin(abs(x),200)
    }
  }
  
  par2tau.fun <- function(x) csttau(par2tau(cstpar(x),family))
  tau2par.fun <- function(x) cstpar(tau2par(csttau(x),family))
  eta2par.fun <- function(x) cstpar(eta2par(x, family))
  par2eta.fun <- function(x) par2eta(cstpar(x),family)
  eta2tau.fun <- function(x) csttau(eta2par(x, 1))
  tau2eta.fun <- function(x) par2eta(csttau(x),1)
  
  if (tau == TRUE) {
    # From transformed parameters to Kendall's tau and 
    # from Kendall's tau to copula parameter
    out$tau <- eta2tau.fun(out$partrans)
    out$partrans <- tau2eta.fun(out$tau)
    out$par <- tau2par.fun(out$tau)
  } else {
    # From transformed parameters to copula parameter
    out$par <- eta2par.fun(out$partrans)
    out$partrans <- par2eta.fun(out$par)
  }
  out$par <- pmin(out$par, 1e5)
  out$par <- pmax(out$par, -1e5)
  
  return(out)
}

## Update trace
"trace.update" <- function(old.par, new.par) {
  
  traces <- median(abs((old.par - new.par)))
  eps <- traces
  
  out <- list()
  out$traces <- traces
  out$eps <- eps
  
  return(out)
}

## Newton-Raphson/Fisher-scoring step
"wz.update" <- function(dd, new.pars, family, method, tau) {
  
  u <- dd$d1 * dd$dpar
  if (tau == TRUE) {
    u <- u * dd$dtau
  }
  if (method == "NR") {
    if (tau == TRUE) {
      w <- dd$dtau^2 * (dd$dpar^2 * (dd$d2/dd$p - dd$d1^2) + dd$dpar2 * dd$d1) + 
        dd$dtau2 * dd$d1 * dd$dpar
    } else {
      w <- dd$dpar^2 * (dd$d2/dd$p - dd$d1^2) + dd$dpar2 * dd$d1
    }
    sel <- !is.na(w) & w > 0
    if (all(!sel)) {
          w <- FisherBiCop(family, new.pars$par, new.pars$par2)
          if (tau == TRUE) {
            w <- (dd$dpar * dd$dtau)^2 * w
          } else {
            w <- dd$dpar^2 * w
          }
    } else {
      w[!sel] <- mean(w[sel])
    }
  } else {
    w <- FisherBiCop(family, new.pars$par, new.pars$par2)
    if (tau == TRUE) {
      w <- (dd$dpar * dd$dtau)^2 * w
    } else {
      w <- dd$dpar^2 * w
    }
  }
  
  qq <- quantile(u/w, c(0.025,0.975))
  z <- new.pars$partrans + pmin(pmax(u/w, qq[1]), qq[2])

  out <- list()
  out$w <- w
  out$z <- z
  
  return(out)
}

## Compute first derivatives of the copula, copula parameters transformations
## and dependence measure transformations
"derivatives.par" <- function(data, new.pars, family, method, tau) {
  
  # Derivatives of the copula with respect to its own parameter
  if (method == "NR") {
    tmp <- bicoppd1d2(data, family, d1 = TRUE, d2 = TRUE)
  } else {
    tmp <- bicoppd1d2(data, family, d1 = TRUE)
  } 
  tmp[1, which(tmp[1, ] == 0)] <- 1e-16
  
  out <- list()
  out$p <- tmp[1, ]
  out$d1 <- tmp[2, ] 

  if (method == "NR") {
    out$d2 <- tmp[3, ]
  }
  
  if (tau == TRUE) {
    # Derivatives of the copula parameter with respect to the dependence measure
    tmp <- dpardtau(new.pars$tau,family)     
    out$dtau <- tmp$d1
    if (method == "NR") {
      out$dtau2 <- tmp$d2
    }
    
    # Derivatives of the dependence measure with respect to the model parameters
    # (equivalent to the derivative of the copula parameter with respect to
    # the calibration function in the Gaussian case, as the link is the same)
    tmp <- dpardeta(new.pars$partrans, 1)
  } else {
    # Derivatives of the copula parameter with respect to the model parameters
    tmp <- dpardeta(new.pars$partrans, family)   
  }
  
  out$dpar <- tmp$d1
  if (method == "NR") {
    out$dpar2 <- tmp$d2
  }
  
  return(out)
}
