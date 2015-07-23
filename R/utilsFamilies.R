## Faster than the implementation from the VineCopula package
FrankPar2Tau <- function(par) {
  tau <- 1 - 4/par + 4/par * debye1(par)
  return(tau)
}

## Faster than the implementation from the VineCopula package
FrankTau2Par <- function(tau) {
  mypar <- function(tau) {
    a <- 1
    if (tau < 0) {
      a <- -1
      tau <- -tau
    }
    a*safeUroot(function(x) tau - FrankPar2Tau(x), 
                interval = c(0 + .Machine$double.eps^0.7, upper = 1e7))$root
  }
  return(sapply(tau,mypar))
}

## Derivative and hessian of the link between Kendall's tau and the Frank
## copula parameter
Frankdtaudpar <- function(par) {
  tmp <- grad(FrankPar2Tau,par, method.args =
                list(eps=1e-4, d=0.0001, 
                     zero.tol=sqrt(.Machine$double.eps/7e-7), 
                     r=2, v=2, show.details=FALSE))
  mm <- splinefun(par,tmp)
  list(d1=tmp,d2=mm(par,deriv=1))
}

## Derivative and hessian of the inverse link between Kendall's tau and the 
## Frank copula parameter
Frankdpardtau <- function(tau) {
  tmp <- grad(FrankTau2Par,tau, method.args =
                list(eps=1e-4, d=0.0001, 
                     zero.tol=sqrt(.Machine$double.eps/7e-7), 
                     r=2, v=2, show.details=FALSE))
  mm <- splinefun(tau,tmp)
  list(d1=tmp,d2=mm(tau,deriv=1))
}

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
  if (!(family %in% c(0,get.familyset()))) {
    return("Copula family not yet implemented.")
  } else if ((family == 1 || family == 2) && abs(par) >= 1) {
    return(paste("The parameter of the Gaussian and",
                 "t-copula has to be in the interval (-1,1)."))
  } else if (family == 5 && par == 0) {
    return("The parameter of the Frank copula has to be different from 0.")
  } else if (family == 2 && par2 <= 2) {
    return(paste("The degrees of freedom parameter of the", 
                 "t-copula has to be larger than 2."))
  } else if ((family == 3 || family == 13) && par <= 0) {
    return("The parameter of the Clayton copula has to be positive.")
  } else if ((family == 4 || family == 14) && par <= 1) {
    return(paste("The parameter of the Gumbel copula",
                 "has to be in the interval (1,oo)."))
  } else if ((family == 23 || family == 33) && par >= 0) {
    return("The parameter of the rotated Clayton copula has to be negative.")
  } else if ((family == 24 || family == 34) && par >= -1) {
    return(paste("The parameter of the rotated Gumbel",
                 "copula has to be in the interval (-oo,-1)."))
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