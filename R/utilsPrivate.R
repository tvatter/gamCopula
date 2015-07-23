## Update parameters with respect to fited gam model
"pars.update" <- function(mm, family, data, tau = TRUE) {
  
  out <- list()
  fitted <- predict(mm, data)
  
  out$partrans <- as.numeric(fitted)
  
  # Define links between Kendall's tau, copula parameter and calibration 
  # function... the cst/cstinv make sure that the boundaries are never attained
  if (family %in% c(3, 4)) {
    cst <- function(x) (1-1e-8)*(1e-8+x)
    cstinv <- function(x) 1e-8+x
    par2tau.fun <- function(x) cst(BiCopPar2Tau(family, cstinv(x)))
    tau2par.fun <- function(x) cstinv(BiCopTau2Par(family, cst(x)))
    eta2par.fun <- function(x) cstinv(BiCopEta2Par(family, x))
    par2eta.fun <- function(x) BiCopPar2Eta(family, x)
    eta2tau.fun <- function(x) cstinv((1+tanh(x/2))/2)
    tau2eta.fun <- function(x) 2*atanh(2*x - 1)
  } else if (family %in% c(1, 2)) {
    cst <- (1-1e-8)
    par2tau.fun <- function(x) BiCopPar2Tau(1, x*cst)*cst
    tau2par.fun <- function(x) BiCopTau2Par(1, x*cst)*cst
    eta2par.fun <- function(x) BiCopEta2Par(1, x)*cst
    par2eta.fun <- function(x) BiCopPar2Eta(1, x)
    eta2tau.fun <- function(x) cst*tanh(x/2)
    tau2eta.fun <- function(x) 2*atanh(x)
  } else {
    cst <- (1-1e-8)
    par2tau.fun <- function(x) BiCopPar2Tau(5, x*cst)*cst
    tau2par.fun <- function(x) BiCopTau2Par(5, x*cst)
    eta2par.fun <- function(x) BiCopEta2Par(5, x)
    par2eta.fun <- function(x) BiCopPar2Eta(5, x)
    eta2tau.fun <- function(x) cst*tanh(x/2)
    tau2eta.fun <- function(x) 2*atanh(x)
  }
  
  if (tau == TRUE) {
    # From transformed parameters to Kendall's tau and 
    # from Kendall's tau to copula parameter
    out$tau <- eta2tau.fun(out$partrans)
    out$partrans <- tau2eta.fun(out$tau)
    out$par <- tau2par.fun(out$tau)
  } else {
    # From transformed parameters to copula parameter
    out$par <- eta2par.fun(out$partrans)
    out$par <- pmin(out$par, 1e5)
    out$partrans <- par2eta.fun(out$par)
  }
  
  return(out)
}

## Update trace
"trace.update" <- function(old.par, new.par) {
  
  traces <- max(abs((old.par - new.par)/old.par))
  eps <- traces
  
  out <- list()
  out$traces <- traces
  out$eps <- eps
  
  return(out)
}

## Newton-Raphson/Fisher-scoring step
"wz.update" <- function(dd, new.pars, family, method, tau = TRUE) {
  
  u <- dd$d1 * dd$dpar
  if (tau == TRUE) {
    u <- u * dd$dtau
  }
  if (method == "NR") {
    if (tau == TRUE) {
      w <- dd$dtau^2 * (dd$dpar^2 * (dd$d2/dd$p - dd$d1^2) + dd$dpar2 * dd$d1) + 
        dd$dtau2 * dd$d1 * dd$dpar * dd$dtau2
    } else {
      w <- dd$dpar^2 * (dd$d2/dd$p - dd$d1^2) + dd$dpar2 * dd$d1
    }
    w <- -w
    w[w<0] <- median(w[w>0])
  } else {
    if (family == 2) {
      if (tau == TRUE) {
        w <- (dd$dpar * dd$dtau)^2 * 
          FisherBiCop(family, new.pars$par, new.pars$par2)
      } else {
        w <- dd$dpar^2 * 
          FisherBiCop(family, new.pars$par, new.pars$par2)
      }
    } else {
      if (tau == TRUE) {
        w <- (dd$dpar * dd$dtau)^2 * 
          FisherBiCop(family, new.pars$par)
      } else {
        w <- dd$dpar^2 * 
          FisherBiCop(family, new.pars$par)
      }
    }
  }
  
  z <- new.pars$partrans + u/w
  
  out <- list()
  out$w <- w
  out$z <- z
  
  return(out)
}

## Compute first derivatives of the copula, copula parameters transformations
## and dependence measure transformations
"derivatives.par" <- function(data, new.pars, family, method, tau = TRUE) {
  
  #cat(paste(paste(range(data[,3]),collapse = ""), "\n"))
  # Derivatives of the copula with respect to its own parameter
  if (family == 2) {
    dp <- apply(data, 1, function(x) c(BiCopPDF(x[1], x[2], family = 2, x[3], 
      x[4]), BiCopDeriv(x[1], x[2], family = 2, x[3], x[4], deriv = "par", 
      log = TRUE)))
    if (method == "NR") {
      temp <- apply(data, 1, function(x) BiCopDeriv2(x[1], x[2], family = 2, 
        x[3], x[4], deriv = "par"))
      dp <- rbind(dp, temp)
    }
  } else {
    dp <- apply(data, 1, function(x) 
      c(BiCopPDF(x[1], x[2], family = family, x[3]), 
        BiCopDeriv(x[1], x[2], family = family, x[3], 
                   deriv = "par", log = TRUE)))
    if (method == "NR") {
      temp <- apply(data, 1, function(x) 
        BiCopDeriv2(x[1], x[2], family = family, x[3], deriv = "par"))
      dp <- rbind(dp, temp)
    }
  }
  
  dp[1, which(dp[1, ] == 0)] <- 1e-16
  
  if (family == 4) {
    sel <- is.nan(dp[2, ])
    if (sum(sel) > 1) {
      dp[2, sel] <- apply(data[sel, ], 1, function(x) derivGumbel(x[3], x[1], 
        x[2]))/dp[1, sel]
    } else if (sum(sel) == 1) {
      dp[2, sel] <- derivGumbel(data[sel, 3], data[sel, 1], 
                                data[sel, 2])/dp[1,sel]
    }
  }
  
  out <- list()
  out$p <- dp[1, ]
  out$d1 <- dp[2, ]
  
  if (method == "NR") {
    out$d2 <- dp[3, ]
  }
  
  if (tau == TRUE) {
    # Derivatives of the dependence measure with respect to the model parameters
    dtau <- 1/(1 + cosh(new.pars$partrans))
    if (method == "NR") {
      dtau2 <- -4 * sinh(new.pars$partrans/2)^4 * (1/sinh(new.pars$partrans))^3
    }
    
    # Derivatives of the copula parameter with respect to the dependence measure
    if (is.element(family, c(1, 2))) {
      dpar <- cos(new.pars$tau * pi/2) * pi/2
      if (method == "NR") {
        dpar2 <- -sin(new.pars$tau * pi/2) * (pi/2)^2
      }
    } else if (family == 3) {
      dpar <- 2/(new.pars$tau - 1)^2
      if (method == "NR") {
        dpar2 <- -4/(new.pars$tau - 1)^3
      }
    } else if (family == 4) {
      dpar <- 1/(new.pars$tau - 1)^2
      if (method == "NR") {
        dpar2 <- -2/(new.pars$tau - 1)^3
      }
    } else if (family == 5) {
      tmp <- Frankdpardtau(new.pars$tau)
      dpar <- tmp$d1
      dpar2 <- tmp$d2
    }
    
    out$dpar <- dpar
    out$dtau <- dtau
    if (method == "NR") {
      out$dtau2 <- dtau2
      out$dpar2 <- dpar2
    }
  } else {
    # Derivatives of the copula parameter with respect to the model parameters
    if (is.element(family, c(1, 2))) {
      dpar <- 1/(1 + cosh(new.pars$partrans))
      if (method == "NR") {
        dpar2 <- -4 * sinh(new.pars$partrans/2)^4 * 
          (1/sinh(new.pars$partrans))^3
      }
    } else if (is.element(family, c(3, 4))) {
      dpar <- exp(new.pars$partrans)
      if (method == "NR") {
        dpar2 <- exp(new.pars$partrans)
      }
    } else {
      dpar <- rep(1, length(new.pars$partrans))
      if (method == "NR") {
        dpar2 <- rep(0, length(new.pars$partrans))
      }
    }
    
    out$dpar <- dpar
    if (method == "NR") {
      out$dpar2 <- dpar2
    }
  }
  
  return(out)
}


fasttau <- function(x, y, weights = NA) {
  m <- length(x)
  n <- length(y)
  if (m == 0 || n == 0) 
    stop("both 'x' and 'y' must be non-empty")
  if (m != n) 
    stop("'x' and 'y' must have the same length")
  ktau <- TauMatrix(matrix(c(x, y), length(x), 2), weights)[2, 1]
  return(ktau)
}

get.modelCount <- function(d) {
  model.count <- rep(0, d^2)
  temp <- 1:(d * (d - 1)/2)
  t1 <- 1
  sel <- seq(d, d^2 - d, by = d)
  for (i in 1:(d - 1)) {
    t2 <- t1 + d - i - 1
    model.count[sel[1:(d - i)] - i + 1] <- temp[t1:t2]
    t1 <- t2 + 1
  }
  model.count <- matrix(model.count, d, d)
}