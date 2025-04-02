## The current familyset of the gamCopula package
get.familyset <- function() {
  # the Frank is available, but it works poorly.... still don't know why
  c(1, 2, 5, 301:304, 401:404)
}

## Fisher information with respect to the Copula parameter for a
## bivariate copula
"FisherBiCop" <- function(family, par, par2 = NULL) {
  if (family == 1) {
    out <- FisherGaussian(par)
  } else if (family == 2) {
    out <- FisherStudentRho(par, par2)
  } else if (family %in% c(301:304)) {
    out <- FisherClayton(abs(par))
  } else if (family %in% c(401:404)) {
    out <- FisherGumbel(abs(par))
  }
  return(out)
}

## Fisher information for the Gumbel copula See Schepsmeier & Stober (2012)
## (contains several typos)
"FisherGumbel" <- function(par) {
  k0 <- 5 / 6 - pi^2 / 18
  E1 <- sapply(par - 1, gsl::expint_E1)

  out <- (1 / par^4) * (par^2 * (-2 / 3 + pi^2 / 9) - par + 2 * k0 / par +
    (par^3 + par^2 + (k0 - 1) * par - 2 * k0 + k0 / par) *
      E1 * exp(par - 1))

  return(out)
}
## Fisher information for the Clayton copula See Oakes (1982) or Schepsmeier &
## Stober (2012) (the latter contains several typos)
"FisherClayton" <- function(par) {
  # browser()
  par <- par + 1
  v1 <- trigamma(1 / (2 * (par - 1)))
  v2 <- trigamma(par / (2 * (par - 1)))
  v3 <- trigamma((2 * par - 1) / (2 * (par - 1)))
  pp <- 1 / ((3 * par - 2) * (2 * par - 1)) *
    (1 + par / (2 * (par - 1)) * (v1 - v2) + 1 / (2 * (par - 1)) * (v2 - v3))
  out <- 1 / par^2 + 2 / (par * (par - 1) * (2 * par - 1)) + 4 * par / (3 * par - 2) -
    2 * (2 * par - 1) * pp / (par - 1)

  return(out)
}

## Fisher information for the Gaussian copula
"FisherGaussian" <- function(par) {
  out <- (1 + par^2) / (1 - par^2)^2
  return(out)
}

## Fisher information for the Student copula
"FisherStudentRho" <- function(par, par2) {
  out <- (2 + par2 + par2 * par^2) / ((4 + par2) * (1 - par^2)^2)
  return(out)
}

# Check parameter(s) for a given copula family
#' @noRd
family.check <- function(family, par, par2 = 0) {
  if (!(family %in% c(0, get.familyset()))) {
    return("Copula family not yet implemented.")
  } else if ((family == 1 || family == 2) && abs(par) >= 1) {
    return(paste(
      "The parameter of the Gaussian and",
      "t-copula has to be in the interval (-1,1)."
    ))
  } else if (family == 5 && par == 0) {
    return("The parameter of the Frank copula has to be different from 0.")
  } else if (family == 2 && par2 <= 2) {
    return(paste(
      "The degrees of freedom parameter of the",
      "t-copula has to be larger than 2."
    ))
  } else if ((family == 3 || family == 13) && par <= 0) {
    return("The parameter of the Clayton copula has to be positive.")
  } else if ((family == 4 || family == 14) && par <= 1) {
    return(paste(
      "The parameter of the Gumbel copula",
      "has to be in the interval (1,oo)."
    ))
  } else if ((family == 23 || family == 33) && par >= 0) {
    return("The parameter of the rotated Clayton copula has to be negative.")
  } else if ((family == 24 || family == 34) && par >= -1) {
    return(paste(
      "The parameter of the rotated Gumbel",
      "copula has to be in the interval (-oo,-1)."
    ))
  } else {
    return(TRUE)
  }
}

## Faster than the implementation from the VineCopula package
FrankPar2Tau <- function(par) {
  tau <- 1 - 4 / par + 4 / par * debye1(par)
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
    a * safeUroot(function(x) tau - FrankPar2Tau(x),
      interval = c(0 + .Machine$double.eps^0.7, upper = 1e7)
    )$root
  }
  return(sapply(tau, mypar))
}

## Derivative and hessian of the link between Kendall's tau and the Frank
## copula parameter
Frankdtaudpar <- function(par) {
  tmp <- grad(FrankPar2Tau, par,
    method.args =
      list(
        eps = 1e-4, d = 0.0001,
        zero.tol = sqrt(.Machine$double.eps / 7e-7),
        r = 2, v = 2, show.details = FALSE
      )
  )
  mm <- splinefun(par, tmp)
  list(d1 = tmp, d2 = mm(par, deriv = 1))
}

## Derivative and hessian of the inverse link between Kendall's tau and the
## Frank copula parameter
Frankdpardtau <- function(tau) {
  tmp <- grad(FrankTau2Par, tau,
    method.args =
      list(
        eps = 1e-4, d = 0.0001,
        zero.tol = sqrt(.Machine$double.eps / 7e-7),
        r = 2, v = 2, show.details = FALSE
      )
  )
  mm <- splinefun(tau, tmp)
  list(d1 = tmp, d2 = mm(tau, deriv = 1))
}

## Link and inverse link functions for all VineCopula package's copula families
links <- function(family, inv = FALSE) {
  if (!inv) {
    f1 <- function(x) tanh(x / 2)
    f2 <- function(x) exp(x)
    f3 <- function(x) 1 + f2(x)
    f4 <- function(x) (tanh(x / 2) + 1) / 2
    f5 <- function(x) 2 + f2(x)
    f2r <- function(x) -f2(x)
    f3r <- function(x) -f3(x)
    f4r <- function(x) -f4(x)
    id <- function(x) x
  } else {
    f1 <- function(x) 2 * atanh(x)
    f2 <- function(x) log(x)
    f3 <- function(x) f2(x - 1)
    f4 <- function(x) 2 * atanh(2 * x - 1)
    f5 <- function(x) f2(x - 2)
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
    Tawn2_270 = list(par = f3r, par2 = f4)
  )
}

## Kendall's tau as a function of the copula parameter
par2tau <- function(x, family) {
  if (family %in% c(1, 2)) {
    return(2 / pi * asin(x))
  }
  if (family %in% c(301:304)) {
    return(x / (2 + abs(x)))
  }
  if (family %in% c(401:404)) {
    return(x / abs(x) - 1 / x)
  }
  if (family == 5) {
    if (is.matrix(x)) {
      return(apply(x, 2, FrankPar2Tau))
    } else {
      return(FrankPar2Tau(x))
    }
  }
}

## Copula parameter as a function of the Kendall's tau
tau2par <- function(x, family) {
  if (family %in% c(1, 2)) {
    return(sin(x * pi / 2))
  }
  if (family %in% c(301:304)) {
    return(2 * x / (1 - abs(x)))
  }
  if (family %in% c(401:404)) {
    y <- abs(x)
    return(x / ((1 - y) * y))
  }
  if (family == 5) {
    if (is.matrix(x)) {
      return(apply(x, 2, FrankTau2Par))
    } else {
      return(FrankTau2Par(x))
    }
  }
}

# First and second derivative of the copula parameter
# with respect to Kendall's tau
dpardtau <- function(x, family) {
  if (family %in% c(1, 2)) {
    d1 <- cos(x * pi / 2) * pi / 2
    d2 <- -sin(x * pi / 2) * (pi / 2)^2
  }
  if (family %in% c(301:304)) {
    y <- abs(x)
    d1 <- 2 / (1 - y)^2
    d2 <- -4 * sign(x) / (y - 1)^3
  }
  if (family %in% c(401:404)) {
    y <- abs(x)
    d1 <- 1 / (1 - y)^2
    d2 <- -2 * sign(x) / (y - 1)^3
  }
  if (family == 5) {
    tmp <- Frankdpardtau(x)
    d1 <- tmp$d1
    d2 <- tmp$d2
  }
  return(list(d1 = d1, d2 = d2))
}

## TODO: change to vectorized versions
# Copula density and its first/second derivatives with respect to the parameter
bicoppd1d2 <- function(x, family, p = TRUE, d1 = FALSE,
                       d2 = FALSE, h = FALSE, hinv = FALSE, cdf = FALSE) {
  if (dim(x)[2] != 4) {
    x <- cbind(x, 0)
  }
  x <- as.matrix(x)

  out <- matrix(0, 0, dim(x)[1])
  myfuns <- list()

  if (p == TRUE) {
    out <- rbind(out, 0)
    myfuns$p <- function(x, family) BiCopPDF(
        u1 = x[, 1], u2 = x[, 2],
        par = x[, 3], par2 = x[, 4],
        family = family, check.pars = FALSE
      )
  }

  if (d1 == TRUE) {
    out <- rbind(out, 0)
    myfuns$d1 <- function(x, family) {
      tmp <- BiCopDeriv(
        u1 = x[, 1], u2 = x[, 2], par = x[, 3], par2 = x[, 4],
        family = family, log = TRUE, check.pars = FALSE
      )
      sel <- is.na(tmp)
      if (sum(sel) != 0) {
        myLogPDF <- function(par)
          log(BiCopPDF(
            u1 = x[sel, 1], u2 = x[sel, 2],
            par = par, par2 = x[sel, 4],
            family = family, check.pars = FALSE
          ))
        tmp[sel] <- grad(myLogPDF, x[sel, 3])
      }
      return(tmp)
    }
  }

  if (d2 == TRUE) {
    out <- rbind(out, 0)
    myfuns$d2 <- function(x, family) {
      tmp <- BiCopDeriv2(
        u1 = x[, 1], u2 = x[, 2], par = x[, 3], par2 = x[, 4],
        family = family, check.pars = FALSE
      )
      sel <- is.na(tmp)
      if (sum(sel) != 0) {
        myPDF <- function(par) BiCopPDF(
            u1 = x[sel, 1], u2 = x[sel, 2],
            par = par, par2 = x[sel, 4],
            family = family, check.pars = FALSE
          )
        tmp[sel] <- grad(function(x) grad(myPDF, x), x[sel, 3])
      }
      return(tmp)
    }
  }

  if (h == TRUE) {
    out <- rbind(out, 0, 0)
    myfuns$h <- function(x, family) do.call(
        rbind,
        BiCopHfunc(
          u1 = x[, 1], u2 = x[, 2],
          par = x[, 3], par2 = x[, 4],
          family = family,
          check.pars = FALSE
        )
      )
  }

  if (hinv == TRUE) {
    out <- rbind(out, 0, 0)
    myfuns$h <- function(x, family) do.call(
        rbind,
        BiCopHinv(
          u1 = x[, 1], u2 = x[, 2],
          par = x[, 3], par2 = x[, 4],
          family = family,
          check.pars = FALSE
        )
      )
  }

  if (cdf == TRUE) {
    out <- rbind(out, 0)
    myfuns$c <- function(x, family) BiCopCDF(
        u1 = x[, 1], u2 = x[, 2],
        par = x[, 3], par2 = x[, 4],
        family = family
      )
  }

  if (family %in% c(0, 1, 2, 5)) {
    out[, 1:dim(out)[2]] <- t(sapply(myfuns, function(f) f(x, family)))
  } else {
    fam <- getFams(family)
    sel <- x[, 3] > 0
    if (sum(sel) > 0) {
      out[, sel] <- t(sapply(myfuns, function(f)
        f(matrix(x[sel, ], nrow = sum(sel), ncol = 4), fam[1])))
    }
    if (sum(!sel) > 0) {
      out[, !sel] <- t(sapply(myfuns, function(f)
        f(matrix(x[!sel, ], nrow = sum(!sel), ncol = 4), fam[2])))
    }
  }

  return(out)
}

## Copula parameter as a function of the calibration function
eta2par <- function(x, family) {
  if (family %in% c(1, 2)) {
    return(tanh(x / 2))
  }
  if (family %in% c(301:304, 5)) {
    return(x)
  }
  if (family %in% c(401:404)) {
    y <- abs(x)
    return(x * (1 + y) / y)
  }
}

## Calibration function as a function of the copula parameter
par2eta <- function(x, family) {
  if (family %in% c(1, 2)) {
    return(2 * atanh(x))
  }
  if (family %in% c(301:304, 5)) {
    return(x)
  }
  if (family %in% c(401:404)) {
    return(x * (1 - 1 / abs(x)))
  }
}

# First and second derivative of the copula parameter
# with respect to the calibration function
dpardeta <- function(x, family) {
  if (family %in% c(1, 2)) {
    d1 <- 1 / (1 + cosh(x))
    d2 <- -4 * sinh(x / 2)^4 * (1 / sinh(x))^3
  } else {
    d1 <- rep(1, length(x))
    d2 <- rep(0, length(x))
  }
  return(list(d1 = d1, d2 = d2))
}

## Return the two component of the double Clayton/Gumbel considered
getFams <- function(family) {
  if (family %in% c(301:304)) {
    fam <- as.numeric(rev(expand.grid(c(23, 33), c(3, 13)))[family - 300, ])
  } else {
    fam <- as.numeric(rev(expand.grid(c(24, 34), c(4, 14)))[family - 400, ])
  }
  return(fam)
}

## Transformation between the Clayton/Gumbel codes from VineCopula and
## the codes from gamCopula
famTrans <- function(fam, inv = TRUE, par = 0, set = FALSE, familyset = NULL) {
  if (inv) {
    if (is.element(fam, c(0, get.familyset()))) {
      return(fam)
    } else {
      if (is.null(familyset)) {
        familyset <- c(301:304, 401:404)
      }
      sel <- (familyset %in% c(301:304, 401:404))
      familyset <- familyset[sel]
      fams <- sapply(familyset, getFams)
      sel <- which(apply(fams, 2, function(x) is.element(fam, x)))
      return(familyset[sel[1]])

      #       if (fam == 3) {
      #         return(301)
      #       } else if (fam == 13) {
      #         return(304)
      #       } else if (fam == 23) {
      #         return(303)
      #       } else if (fam == 33) {
      #         return(302)
      #       } else if (fam == 4) {
      #         return(401)
      #       } else if (fam == 14) {
      #         return(404)
      #       } else if (fam == 24) {
      #         return(403)
      #       } else if (fam == 34) {
      #         return(402)
      #       }
    }
  } else {
    if (set == FALSE) {
      if (par > 0 && fam %in% c(301, 302)) {
        return(3)
      } else if (par > 0 && fam %in% c(303, 304)) {
        return(13)
      } else if (par < 0 && fam %in% c(301, 303)) {
        return(23)
      } else if (par < 0 && fam %in% c(302, 304)) {
        return(33)
      } else if (par > 0 && fam %in% c(401, 402)) {
        return(4)
      } else if (par > 0 && fam %in% c(403, 404)) {
        return(14)
      } else if (par < 0 && fam %in% c(401, 403)) {
        return(24)
      } else if (par < 0 && fam %in% c(402, 404)) {
        return(34)
      } else {
        return(fam)
      }
    } else {
      sel <- fam %in% c(301:304, 401:404)
      if (sum(sel) != 0) {
        return(sort(c(fam[!sel], unique(as.numeric(sapply(fam[sel], getFams))))))
      } else {
        return(fam)
      }
    }
  }
}

## Find rotations for the double Clayton/Gumbel
getRotations <- function(i) {
  out <- i
  if (i %in% c(301, 302, 303, 304)) {
    out <- c(301, 302, 303, 304)
  }
  if (i %in% c(401, 402, 403, 404)) {
    out <- c(401, 402, 403, 404)
  }
  out
}

## Helper function to obtain all unique rotations for the Clayton/Gumbel
withRotations <- function(nums) {
  unique(unlist(lapply(nums, getRotations)))
}

## Alternative to BiCopName from the VineCopula package for the double
## Clayton/Gumbel
bicopname <- function(family) {
  tmp <- c("standard", "survival", "90 degrees rotated", "270 degrees rotated")
  if (family %in% c(301:304)) {
    nn <- paste("Clayton type", family - 300, collapse = "")
    tmp <- tmp[as.numeric(rev(expand.grid(c(3, 4), c(1, 2)))[family - 300, ])]
    nn <- paste0(nn, " (", paste(tmp, collapse = " and "), ")",
      collapse = ""
    )
  } else if (family %in% c(401:404)) {
    nn <- paste("Gumbel type", family - 400, collapse = "")
    tmp <- tmp[as.numeric(rev(expand.grid(c(3, 4), c(1, 2)))[family - 400, ])]
    nn <- paste0(nn, " (", paste(tmp, collapse = " and "), ")",
      collapse = ""
    )
  } else {
    nn <- BiCopName(family, short = FALSE)
  }
  return(nn)
}
