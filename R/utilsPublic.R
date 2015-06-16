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

## From calibration function to copula parameter
BiCopEta2Par <- function(family, eta) {
  
  if (!is.element(family, c(1, 2, 3, 4, 5, 13, 14, 23, 24, 33, 34))) 
    stop("Copula family not yet implemented.")
  
  if (!is.double(eta)) 
    stop("Eta should be real.")
  
  if (is.element(family, c(1, 2))) {
    return(tanh(eta/2))
  } else if (is.element(family, c(3, 13))) {
    return(exp(eta))
  } else if (is.element(family, c(4, 14))) {
    return(1 + exp(eta))
  } else if (is.element(family, c(23, 33))) {
    return(-exp(eta))
  } else if (is.element(family, c(24, 34))) {
    return(-1 - exp(eta))
  } else {
    return(eta)
  }
}

## From copula parameter to calibration function
BiCopPar2Eta <- function(family, par) {
  
  if (!is.element(family, c(1, 2, 3, 4, 5, 13, 14, 23, 24, 33, 34))) 
    stop("Copula family not yet implemented.")
  
  if (!is.double(par)) 
    stop("Par should be real.")
  
  if (is.element(family, c(1, 2))) {
    if (any(abs(par) >= 1)) {
      stop("Par should be in (-1,1) for the Gaussian and t copulas.")
    } else {
      return(2 * atanh(par))
    }
  } else if (is.element(family, c(3, 13))) {
    if (any(par <= 0)) {
      stop("Par should be greater than zero for the Clayton and survival Clayton copulas.")
    } else {
      return(log(par))
    }
  } else if (is.element(family, c(4, 14))) {
    if (any(par <= 1)) {
      stop("Par should be greater than one for the Gumbel and survival Gumbel copulas.")
    } else {
      return(log(par - 1))
    }
  } else if (is.element(family, c(23, 33))) {
    if (any(par >= 0)) {
      stop("Par should be SMALLER than zero for the 90 and 270 rotated Clayton copulas.")
    } else {
      return(log(-par))
    }
  } else if (is.element(family, c(24, 34))) {
    if (any(par >= -1)) {
      stop("Par should be smaller than minus one for the 90 and 270 Gumbel copulas.")
    } else {
      return(log(-1 - par))
    }
  } else {
    return(par)
  }
}


## Fisher information with respect to the Copula parameter for a bivariate copula
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

## Fisher information for the Gumbel copula See Schepsmeier & Stöber (2012)
## (contains several typos)
"FisherGumbel" <- function(par) {
  k0 <- 5/6 - pi^2/18
  E1 <- sapply(par - 1, gsl::expint_E1)
  
  out <- (1/par^4) * (par^2 * (-2/3 + pi^2/9) - par + 2 * k0/par + (par^3 + par^2 + 
    (k0 - 1) * par - 2 * k0 + k0/par) * E1 * exp(par - 1))
  
  return(out)
}


## Fisher information for the Clayton copula See Oakes (1982) or Schepsmeier &
## Stöber (2012) (the latter contains several typos)
"FisherClayton" <- function(par) {
  # browser()
  par <- par + 1
  v1 <- trigamma(1/(2 * (par - 1)))
  v2 <- trigamma(par/(2 * (par - 1)))
  v3 <- trigamma((2 * par - 1)/(2 * (par - 1)))
  pp <- 1/((3 * par - 2) * (2 * par - 1)) * (1 + par/(2 * (par - 1)) * (v1 - v2) + 
    1/(2 * (par - 1)) * (v2 - v3))
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


## Simulation from a conditional copula given a calibration function of covariates
## (modified from Elif Acar's code)
CondBiCopSim <- function(family, calib.fnc, X, 
                         par2 = 0, return.par = TRUE, tau = FALSE) {
  
  if (!(family %in% c(1, 2, 3, 4, 13, 14, 23, 24, 33, 34))) 
    stop("Copula family not yet implemented.")
  
  ## univariate covariate:
  if (is.vector(X)) {
    dim.x <- 1
    dim.arg <- length(formals(calib.fnc))
    
    if (dim.x != dim.arg) 
      stop(paste("Calibration function requires", dim.arg, "covariates."))
    
    eta <- sapply(X, function(x) calib.fnc(x))
  }
  
  ## multivariate covariate:
  if (is.matrix(X) || is.data.frame(X)) {
    dim.x <- ncol(X)
    dim.arg <- length(formals(calib.fnc))
    if (dim.x != dim.arg) 
      stop(paste("Calibration function requires", dim.arg, "covariates."))
    
    fnc.expr <- paste("calib.fnc(", paste(paste("x[", 1:dim.arg, "]", sep = ""), 
      collapse = ","), ")")
    fnc.expr <- parse(text = fnc.expr)
    eta <- apply(X, 1, function(x) eval(fnc.expr))
  }
  
  if (tau != TRUE) {
    par <- BiCopEta2Par(family, eta)
  } else {
    if (is.element(family, c(1, 2))) {
      par <- sapply(tanh(eta/2), function(x) BiCopTau2Par(1, x))
    } else {
      par <- sapply((1 + tanh(eta/2))/2, function(x) BiCopTau2Par(family, x))
    }
  }

  data <- t(mapply(BiCopSim, 1, family, par, par2))
  if (return.par == TRUE) {
    return(list(data = data, par = par, eta = eta))
  } else {
    return(data)
  }
} 

## Link and inverse link functions for all VineCopula package's copula families
links <- function(family, inv = FALSE) {
  if (!inv) {
    f1 <- function(x) tanh(x/2)
    f2 <- function(x) exp(x)
    f3 <- function(x) 1 + f2(x)
    f4 <- function(x) (tanh(x/2)+1)/2
    f5 <- function(x) 2 + f2(x)    
    f2r <- function(x) -f2(-x)
    f3r <- function(x) -f3(-x)
    f4r <- function(x) -f4(-x)
    id <- function(x) x
  } else {
    f1 <- function(x) 2*atanh(x)
    f2 <- function(x) log(x)
    f3 <- function(x) f2(x-1)
    f4 <- function(x) 2*atanh(2*x-1)
    f5 <- function(x) f2(x-2)
    f2r <- function(x) -f2(-x)
    f3r <- function(x) -f3(-x)
    f4r <- function(x) -f4(-x)
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