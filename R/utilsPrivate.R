
## Fisher information with respect to the Copula parameter for a
## bivariate copula
"FisherBiCop" <- function(family, par, par2 = NULL) {
  
  if (family == 1) {
    out <- FisherGaussian(par)
  } else if (family == 2) {
    out <- FisherStudentRho(par, par2)
  } else if (family == 3) {
    out <- FisherClayton(par)
  } else if (family == 4) {
    out <- FisherGumbel(par)
  }
  
  return(out)
}


## Fisher information with respect to Kendall's tau for a bivariate copula
"FisherBiCop2" <- function(family, tau, par2 = NULL) {
  
  if (family == 1) {
    par <- sin(pi/2 * tau)
    out <- FisherGaussian(par) * (pi * cos(pi * tau/2)/2)^2
  } else if (family == 2) {
    par <- sin(pi/2 * tau)
    out <- FisherStudentRho(par, par2) * (pi * cos(pi * tau/2)/2)^2
  } else if (family == 3) {
    par <- 2 * tau/(1 - tau)
    out <- FisherClayton(par) * (2/(tau - 1)^2)^2
  } else if (family == 4) {
    par <- 1/(1 - tau)
    out <- FisherGumbel(par) * (1/(tau - 1)^2)^2
  }
  
  return(out)
}

## Fisher information for the Gumbel copula See Schepsmeier & Stober (2012)
## (contains several typos)
"FisherGumbel" <- function(par) {
  k0 <- 5/6 - pi^2/18
  E1 <- sapply(par - 1, gsl::expint_E1)
  
  out <- (1/par^4) * (par^2 * (-2/3 + pi^2/9) - par + 2 * k0/par + 
                        (par^3 + par^2 + (k0 - 1) * par - 2 * k0 + k0/par) *
                        E1 * exp(par - 1))
  
  return(out)
}
## Fisher information for the Clayton copula See Oakes (1982) or Schepsmeier &
## Stober (2012) (the latter contains several typos)
"FisherClayton" <- function(par) {
  # browser()
  par <- par + 1
  v1 <- trigamma(1/(2 * (par - 1)))
  v2 <- trigamma(par/(2 * (par - 1)))
  v3 <- trigamma((2 * par - 1)/(2 * (par - 1)))
  pp <- 1/((3 * par - 2) * (2 * par - 1)) * 
    (1 + par/(2 * (par - 1)) * (v1 - v2) + 1/(2 * (par - 1)) * (v2 - v3))
  out <- 1/par^2 + 2/(par * (par - 1) * (2 * par - 1)) + 4 * par/(3 * par - 2) - 
    2 * (2 * par - 1) * pp/(par - 1)
  
  return(out)
}

## Fisher information for the Gaussian copula
"FisherGaussian" <- function(par) {
  
  out <- (1 + par^2)/(1 - par^2)^2
  return(out)
}

## Fisher information for the Student copula
"FisherStudentRho" <- function(par, par2) {
  
  out <- (2 + par2 + par2 * par^2)/((4 + par2) * (1 - par^2)^2)
  return(out)
}

## Update parameters with respect to fited gam model
"pars.update" <- function(mm, family, data, tau = TRUE) {
  
  out <- list()
  fitted <- predict(mm, data)
  
  out$partrans <- as.numeric(fitted)
  
  if (tau == TRUE) {
    # From transformed parameters to Kendall's tau
    if (is.element(family, c(1, 2, 5))) {
      out$tau <- tanh(out$partrans/2)
      out$tau <- pmax(-1 + 1e-08, out$tau)
      out$tau <- pmin(1 - 1e-08, out$tau)
      out$partrans <- 2 * atanh(out$tau)
    } else {
      out$tau <- (1 + tanh(out$partrans/2))/2
      out$tau <- pmax(1e-08, out$tau)
      out$tau <- pmin(1 - 1e-08, out$tau)
      out$partrans <- 2 * atanh(2 * out$tau - 1)
    }
    
    # From Kendall's tau to copula parameter
    if (is.element(family, c(1, 2, 5))) {
      if (is.element(family, c(1, 2))) {
        out$par <- sapply(out$tau, function(x) BiCopTau2Par(1, x))
      } else {
        out$par <- sapply(out$tau, function(x) BiCopTau2Par(5, x))
      }
      sel1 <- out$par == 1
      sel2 <- out$par == -1
      out$par[sel1] <- 1 - 1e-08
      out$par[sel2] <- -1 + 1e-08
    } else {
      out$par <- sapply(out$tau, function(x) BiCopTau2Par(family, x))
      if (family == 3) {
        out$par <- pmax(1e-08, out$par)
      }
      if (family == 4) {
        out$par <- pmax(1 + 1e-08, out$par)
      }
    }
  } else {
    # From transformed parameters to copula parameter
    if (is.element(family, c(1, 2))) {
      out$par <- tanh(out$partrans/2)
      sel1 <- out$par == 1
      sel2 <- out$par == -1
      out$par[sel1] <- 1 - 1e-08
      out$par[sel2] <- -1 + 1e-08
      out$partrans <- 2 * atanh(out$par)
    } else if (family == 3) {
      out$par <- exp(out$partrans)
      sel1 <- out$par == Inf
      sel2 <- out$par == 0
      out$par[sel1] <- 1e+05
      out$par[sel2] <- 1e-08
      out$partrans <- log(out$par)
    } else if (family == 4) {
      out$par <- 1 + exp(out$partrans)
      sel1 <- out$par == Inf
      sel2 <- out$par == 1
      out$par[sel1] <- 1e+05
      out$par[sel2] <- 1 + 1e-08
      out$partrans <- log(out$par - 1)
    } else {
      out$par <- out$partrans
      sel1 <- out$par == Inf
      sel2 <- out$par == 0
      sel3 <- out$par == -Inf
      out$par[sel1] <- 1e+05
      out$par[sel2] <- 1e-08
      out$par[sel3] <- -1e+05
      out$partrans <- out$par
    }
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
        w <- (dd$dpar * dd$dtau)^2 * FisherBiCop(family, new.pars$par, new.pars$par2)
      } else {
        w <- dd$dpar^2 * FisherBiCop(family, new.pars$par, new.pars$par2)
      }
    } else {
      if (tau == TRUE) {
        w <- (dd$dpar * dd$dtau)^2 * FisherBiCop(family, new.pars$par)
      } else {
        w <- dd$dpar^2 * FisherBiCop(family, new.pars$par)
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

## Derivative of the Gumbel copula density with respect to the parameter (Higher
## precision than in the VineCopula package)
derivGumbel <- function(s, u, v) {
  a <- -log(u)
  b <- -log(v)
  A <- a^s + b^s
  x <- -1 + s
  y <- (-3 + 1/s)
  
  temp <- 1/(u * v * s^2) * exp(-A^(1/s))
  temp <- (a^(x/y) * A * b^(x/y))^y * temp
  
  temp2 <- s * (-a^s * ((-1 + s)^2 + (-3 + 2 * s) * A^(1/s) + A^(2/s)) + 
                  s * (-1 +  s + A^(1/s)) * b^s) * log(a)
  temp2 <- c(temp2, (1 - s + (-3 + s) * A^(1/s) + A^(2/s)) * A * log(A))
  temp2 <- c(temp2, s * (s * a^s * (1 + (-1 + s + A^(1/s)) * log(b)) + 
                           b^s * (s -  ((-1 + s)^2 + (-3 + 2 * s) * A^(1/s) +
                                          A^(2/s)) * log(b))))
  
  return(temp * sum(temp2))
}

# Check parameter(s) for a given copula family
family.check <- function(family, par, par2 = 0) {
  if (!(family %in% c(0, 1, 2, 3, 4, 13, 14, 23, 24, 33, 34))) {
    return("Copula family not yet implemented.")
  } else if ((family == 1 || family == 2) && abs(par) >= 1) {
    return(paste("The parameter of the Gaussian and",
                 "t-copula has to be in the interval (-1,1)."))
  } else if (family == 2 && par2 <= 2) {
    return(paste("The degrees of freedom parameter of the", 
                 "t-copula has to be larger than 2."))
  } else if ((family == 3 || family == 13) && par <= 0) {
    return("The parameter of the Clayton copula has to be positive.")
  } else if ((family == 4 || family == 14) && par < 1) {
    return(paste("The parameter of the Gumbel copula",
                 "has to be in the interval [1,oo)."))
  } else if ((family == 23 || family == 33) && par >= 0) {
    return("The parameter of the rotated Clayton copula has to be negative.")
  } else if ((family == 24 || family == 34) && par > -1) {
    return(paste("The parameter of the rotated Gumbel",
                 "copula has to be in the interval (-oo,-1]."))
  } else {
    return(TRUE)
  }
}

getRotations <- function (i) {
  out <- i
  if (i %in% c(3, 13, 23, 33)) 
    out <- c(3, 13, 23, 33)
  if (i %in% c(4, 14, 24, 34)) 
    out <- c(4, 14, 24, 34)
  if (i %in% c(6, 16, 26, 36)) 
    out <- c(6, 16, 26, 36)
  if (i %in% c(7, 17, 27, 37)) 
    out <- c(7, 17, 27, 37)
  if (i %in% c(8, 18, 28, 38)) 
    out <- c(8, 18, 28, 38)
  if (i %in% c(9, 19, 29, 39)) 
    out <- c(9, 19, 29, 39)
  if (i %in% c(10, 20, 30, 40)) 
    out <- c(10, 20, 30, 40)
  if (i %in% c(104, 114, 124, 134)) 
    out <- c(104, 114, 124, 134)
  if (i %in% c(204, 214, 224, 234)) 
    out <- c(204, 214, 224, 234)
  out
}



## Rotate data
rotate.data <- function(data, angle = 0) {
  
  valid <- c(0, 90, 180, 270)
  if (!is.element(angle, valid)) {
    print("Angle must be 0, 90, 180 or 270.")
    return(data)
  }
  
  n <- dim(data)
  if (n[2] != 2) {
    print("Data must be bivariate.")
    return(data)
  }
  
  if (angle == 90) {
    return(cbind(1 - data[, 1], data[, 2]))
  } else if (angle == 180) {
    return(cbind(1 - data[, 1], 1 - data[, 2]))
  } else if (angle == 270) {
    return(cbind(data[, 1], 1 - data[, 2]))
  } else {
    return(data)
  }
}

## Rotate copula parameter
rotate.param <- function(param, angle = 0) {
  
  valid <- c(0, 90, 180, 270)
  if (!is.element(angle, valid)) {
    print("Angle must be 0, 90, 180 or 270.")
    return(param)
  }
  
  if (angle == 90) {
    return(-param)
  } else if (angle == 180) {
    return(param)
  } else if (angle == 270) {
    return(-param)
  } else {
    return(param)
  }
}


withRotations <- function(nums) {
  unique(unlist(lapply(nums, getRotations)))
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


## Link and inverse link functions for all VineCopula package's copula families
links <- function(family, inv = FALSE) {
  if (!inv) {
    f1 <- function(x) tanh(x/2)
    f2 <- function(x) exp(x)
    f3 <- function(x) 1 + f2(x)
    f4 <- function(x) (tanh(x/2)+1)/2
    f5 <- function(x) 2 + f2(x)    
    f2r <- function(x) -f2(x)
    f3r <- function(x) -f3(x)
    f4r <- function(x) -f4(x)
    id <- function(x) x
  } else {
    f1 <- function(x) 2*atanh(x)
    f2 <- function(x) log(x)
    f3 <- function(x) f2(x-1)
    f4 <- function(x) 2*atanh(2*x-1)
    f5 <- function(x) f2(x-2)
    f2r <- function(x) -f2(x)
    f3r <- function(x) -f3(x)
    f4r <- function(x) -f4(x)
    id <- function(x) x
  }
  
  switch(BiCopName(family),
         N = list(par = f1, par2 = NULL),
         t = list(par = f1, par2 = f5),
         C = list(par = f2, par2 = NULL),
         G = list(par = f3, par2 = NULL),
         F = list(par = id, par2 = NULL),
         J = list(par = f3, par2 = NULL),
         BB1 = list(par = f2, par2 = f3),
         BB6 = list(par = f3, par2 = f3),
         BB7 = list(par = f3, par2 = f2),
         BB8 = list(par = f3, par2 = f4),
         SC = list(par = f2, par2 = NULL),
         SG = list(par = f3, par2 = NULL),
         SJ = list(par = f3, par2 = NULL),
         SBB1 = list(par = f2, par2 = f3),
         SBB6 = list(par = f3, par2 = f3),
         SBB7 = list(par = f3, par2 = f2),
         SBB8 = list(par = f3, par2 = f4),
         C90 = list(par = f2r, par2 = NULL),
         G90 = list(par = f3r, par2 = NULL),
         F90 = list(par = f3r, par2 = NULL),
         J90 = list(par = f3r, par2 = NULL),
         BB1_90 = list(par = f2r, par2 = f3r),
         BB6_90 = list(par = f3r, par2 = f3r),
         BB7_90 = list(par = f3r, par2 = f2r),
         BB8_90 = list(par = f3r, par2 = f4r),
         C270 = list(par = f2r, par2 = NULL),
         G270 = list(par = f3r, par2 = NULL),
         F270 = list(par = f3r, par2 = NULL),
         J270 = list(par = f3r, par2 = NULL),
         BB1_270 = list(par = f2r, par2 = f3r),
         BB6_270 = list(par = f3r, par2 = f3r),
         BB7_270 = list(par = f3r, par2 = f2r),
         BB8_270 = list(par = f3r, par2 = f4r),
         Tawn = list(par = f3, par2 = f4),
         Tawn180 = list(par = f3, par2 = f4),
         Tawn90 = list(par = f3r, par2 = f4),
         Tawn270 = list(par = f3r, par2 = f4),
         Tawn2 = list(par = f3, par2 = f4),
         Tawn2_180 = list(par = f3, par2 = f4),
         Tawn2_90 = list(par = f3r, par2 = f4),
         Tawn2_270 = list(par = f3r, par2 = f4))
}