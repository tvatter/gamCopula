#' Copula Parameter of a Bivariate Copula for a Given Value of the Calibration
#' Function
#'
#' Computes the (first) copula parameter of a bivariate copula for a given
#' value of the calibration function (eta).
#'
#' @param family A copula family:
#' \code{1} Gaussian,
#' \code{2} Student t,
#' \code{301} Double Clayton type I (standard and rotated 90 degrees),
#' \code{302} Double Clayton type II (standard and rotated 270 degrees),
#' \code{303} Double Clayton type III (survival and rotated 90 degrees),
#' \code{304} Double Clayton type IV (survival and rotated 270 degrees),
#' \code{401} Double Gumbel type I (standard and rotated 90 degrees),
#' \code{402} Double Gumbel type II (standard and rotated 270 degrees),
#' \code{403} Double Gumbel type III (survival and rotated 90 degrees),
#' \code{404} Double Gumbel type IV (survival and rotated 270 degrees).
#' @param eta The calibration function.
#' @return The value of the first copula parameter, depending on the copula
#' parameter and family as:
#' \itemize{
#' \item \code{1} Gaussian, \code{f(x) = tanh(x/2)}
#' \item \code{2} Student t, \code{f(x) = tanh(x/2)}
#' \item \code{301} Double Clayton type I (standard and rotated 90 degrees),
#' \code{f(x) = x}
#' \item \code{302} Double Clayton type II (standard and rotated 270 degrees),
#' \code{f(x) = x}
#' \item \code{303} Double Clayton type III (survival and rotated 90 degrees),
#' \code{f(x) = x}
#' \item \code{304} Double Clayton type IV (survival and rotated 270 degrees),
#' \code{f(x) = x}
#' \item \code{401} Double Gumbel type I (standard and rotated 90 degrees),
#' \code{f(x) = x*(1+abs(x))/abs(x)}
#' \item \code{402} Double Gumbel type II (standard and rotated 270 degrees),
#' \code{f(x) = x*(1+abs(x))/abs(x)}
#' \item \code{403} Double Gumbel type III (survival and rotated 90 degrees),
#' \code{f(x) = x*(1+abs(x))/abs(x)}
#' \item \code{404} Double Gumbel type IV (survival and rotated 270 degrees)
#' \code{f(x) = x*(1+abs(x))/abs(x)}.
#' }
#' @seealso \code{\link{BiCopPar2Eta}}, \code{\link[VineCopula]{BiCopPar2Tau}},
#' and \code{\link[VineCopula]{BiCopTau2Par}} from the \pkg{VineCopula} package.
#' @name BiCopEta2Par
#' @rdname BiCopEta2Par
#' @export
BiCopEta2Par <- function(family, eta) {
  if (length(family) != 1) {
    stop("Family has to be a scalar/integer.")
  }

  if (!is.element(family, c(1, 2, 3, 4, 5, 13, 14, 23, 24, 33, 34))) {
    stop("Copula family not yet implemented.")
  }

  if (missing(eta)) {
    stop("Eta not set.")
  }

  if (!is.double(eta)) {
    stop("Eta should be real.")
  }

  ll <- links(family)
  return(ll$par(eta))
}

#' Calibration Function of a Bivariate Copula for a Given Parameter's Value
#'
#' Computes the calibration function (eta) of a bivariate copula for a given
#' value of the (first) copula parameter.
#'
#' @param family A copula family:
#' \code{1} Gaussian,
#' \code{2} Student t,
#' \code{301} Double Clayton type I (standard and rotated 90 degrees),
#' \code{302} Double Clayton type II (standard and rotated 270 degrees),
#' \code{303} Double Clayton type III (survival and rotated 90 degrees),
#' \code{304} Double Clayton type IV (survival and rotated 270 degrees),
#' \code{401} Double Gumbel type I (standard and rotated 90 degrees),
#' \code{402} Double Gumbel type II (standard and rotated 270 degrees),
#' \code{403} Double Gumbel type III (survival and rotated 90 degrees),
#' \code{404} Double Gumbel type IV (survival and rotated 270 degrees).
#' @param par The (first) copula parameter
#' @return The value of the calibration function, depending on the copula
#' parameter and family as:
#' \itemize{
#' \item \code{1} Gaussian, \code{f(x) = 2*atanh(x)}
#' \item \code{2} Student t, \code{f(x) = 2*atanh(x)}
#' \item \code{301} Double Clayton type I (standard and rotated 90 degrees),
#' \code{f(x) = x}
#' \item \code{302} Double Clayton type II (standard and rotated 270 degrees),
#' \code{f(x) = x}
#' \item \code{303} Double Clayton type III (survival and rotated 90 degrees),
#' \code{f(x) = x}
#' \item \code{304} Double Clayton type IV (survival and rotated 270 degrees),
#' \code{f(x) = x}
#' \item \code{401} Double Gumbel type I (standard and rotated 90 degrees),
#' \code{f(x) = x*(1-1/abs(x))}
#' \item \code{402} Double Gumbel type II (standard and rotated 270 degrees),
#' \code{f(x) = x*(1-1/abs(x))}
#' \item \code{403} Double Gumbel type III (survival and rotated 90 degrees),
#' \code{f(x) = x*(1-1/abs(x))}
#' \item \code{404} Double Gumbel type IV (survival and rotated 270 degrees)
#' \code{f(x) = x*(1-1/abs(x))}.
#' }
#' @seealso  \code{\link{BiCopEta2Par}},
#' \code{\link[VineCopula]{BiCopPar2Tau}},
#' and \code{\link[VineCopula]{BiCopTau2Par}} from the \pkg{VineCopula} package.
#' @name BiCopPar2Eta
#' @rdname BiCopPar2Eta
#' @export
BiCopPar2Eta <- function(family, par) {
  if (length(family) != 1) {
    stop("Input for family has to be a scalar/integer.")
  }

  if (!is.element(family, c(1, 2, 3, 4, 5, 13, 14, 23, 24, 33, 34))) {
    stop("Copula family not yet implemented.")
  }

  if (missing(par)) {
    stop("Par not set.")
  }

  if (!is.double(par)) {
    stop("Par should be real.")
  }

  if (is.element(family, c(1, 2)) && any(abs(par) >= 1)) {
    stop("Par should be in (-1,1) for the Gaussian and t copulas.")
  }

  if (is.element(family, c(3, 13)) && any(par <= 0)) {
    stop(paste(
      "Par should be greater than zero for the Clayton",
      "and survival Clayton copulas."
    ))
  }

  if (is.element(family, c(4, 14)) && any(par <= 1)) {
    stop(paste(
      "Par should be greater than one for the Gumbel",
      "and survival Gumbel copulas."
    ))
  }

  if (is.element(family, c(23, 33)) && any(par >= 0)) {
    stop(paste(
      "Par should be smaller than zero for",
      "the 90 and 270 rotated Clayton copulas."
    ))
  }

  if (is.element(family, c(24, 34)) && any(par >= -1)) {
    stop(paste(
      "Par should be smaller than minus one",
      "for the 90 and 270 Gumbel copulas."
    ))
  }

  ll <- links(family, TRUE)
  return(ll$par(par))
}


#' Simulation from a Conditional Bivariate Copula
#'
#' Simulates from a conditional bivariate copula, where each copula parameter
#' takes a different value, depending on the calibration
#' function and covariates.
#'
#' @param family  family A copula family:
#' \code{1} Gaussian,
#' \code{2} Student t,
#' \code{3} Clayton,
#' \code{4} Gumbel,
#' \code{5} Frank,
#' \code{13} Survival Clayton,
#' \code{14} Survival Gumbel,
#' \code{23} Rotated (90 degrees) Clayton,
#' \code{24} Rotated (90 degrees) Gumbel,
#' \code{33} Rotated (270 degrees) Clayton and
#'  \code{34} Rotated (270 degrees) Gumbel.
#' @param calib.fnc A calibration function.
#' @param X A vector (if \code{calib.fnc} takes a single argument) or matrix
#' (if \code{calib.fnc} takes multiple arguments) of covariates values.
#' @param par2 The second copula parameter (for the Student t), default
#' \code{par2 = 0}.
#' @param return.par Should the parameter (and calibration function) be returned
#' as well (default \code{return.par = TRUE})?
#' @param tau Should the calibration function (and the model) be specified for
#' the copula parameter or Kendall's tau (default \code{tau = TRUE})?
#' @return If \code{return.par = TRUE}, then the function returns a list with:
#' \itemize{
#' \item \code{data}, a matrix with two columns containing the simulated data,
#' \item \code{par}, a vector containing the values of the copula parameter,
#' \item and \code{eta}, a vector containing the values of the
#' calibration function.
#' }
#' If \code{return.par = FALSE}, then the function simply returns \code{data},
#' a matrix with two columns containing the simulated data.
#' @seealso \code{\link{gamBiCopFit}} and \code{\link{gamBiCopSimulate}}.
#' @examples
#' require(copula)
#' set.seed(0)
#'
#' ## Simulation parameters (sample size, correlation between covariates,
#' ## Gaussian copula family)
#' n <- 2e2
#' rho <- 0.5
#' fam <- 1
#'
#'
#' ## A calibration surface depending on three variables
#' eta0 <- 1
#' calib.surf <- list(
#'   calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
#'     Tm <- (Tf - Ti) / 2
#'     a <- -(b / 3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
#'     return(a + b * (t - Tm)^2)
#'   },
#'   calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1) {
#'     a <- b * (1 - 2 * Tf * pi / (f * Tf * pi +
#'       cos(2 * f * pi * (Tf - Ti))
#'       - cos(2 * f * pi * Ti)))
#'     return((a + b) / 2 + (b - a) * sin(2 * f * pi * (t - Ti)) / 2)
#'   },
#'   calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf / 8) {
#'     Tm <- (Tf - Ti) / 2
#'     a <- (b * s * sqrt(2 * pi) / Tf) * (pnorm(0, Tm, s) - pnorm(Tf, Tm, s))
#'     return(a + b * exp(-(t - Tm)^2 / (2 * s^2)))
#'   }
#' )
#'
#' ## Display the calibration surface
#' par(mfrow = c(1, 3), pty = "s", mar = c(1, 1, 4, 1))
#' u <- seq(0, 1, length.out = 100)
#' sel <- matrix(c(1, 1, 2, 2, 3, 3), ncol = 2)
#' jet.colors <- colorRamp(c(
#'   "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
#'   "yellow", "#FF7F00", "red", "#7F0000"
#' ))
#' jet <- function(x) rgb(jet.colors(exp(x / 3) / (1 + exp(x / 3))),
#'     maxColorValue = 255
#'   )
#' for (k in 1:3) {
#'   tmp <- outer(u, u, function(x, y)
#'     eta0 + calib.surf[[sel[k, 1]]](x) + calib.surf[[sel[k, 2]]](y))
#'   persp(u, u, tmp,
#'     border = NA, theta = 60, phi = 30, zlab = "",
#'     col = matrix(jet(tmp), nrow = 100),
#'     xlab = paste("X", sel[k, 1], sep = ""),
#'     ylab = paste("X", sel[k, 2], sep = ""),
#'     main = paste("eta0+f", sel[k, 1],
#'       "(X", sel[k, 1], ") +f", sel[k, 2],
#'       "(X", sel[k, 2], ")",
#'       sep = ""
#'     )
#'   )
#' }
#'
#' ## 3-dimensional matrix X of covariates
#' covariates.distr <- mvdc(normalCopula(rho, dim = 3),
#'   c("unif"), list(list(min = 0, max = 1)),
#'   marginsIdentical = TRUE
#' )
#' X <- rMvdc(n, covariates.distr)
#'
#' ## U in [0,1]x[0,1] with copula parameter depending on X
#' U <- condBiCopSim(fam, function(x1, x2, x3) {
#'   eta0 + sum(mapply(function(f, x)
#'     f(x), calib.surf, c(x1, x2, x3)))
#' }, X[, 1:3], par2 = 6, return.par = TRUE)
#'
#' ## Merge U and X
#' data <- data.frame(U$data, X)
#' names(data) <- c(paste("u", 1:2, sep = ""), paste("x", 1:3, sep = ""))
#'
#' ## Display the data
#' dev.off()
#' plot(data[, "u1"], data[, "u2"], xlab = "U1", ylab = "U2")
#' @name condBiCopSim
#' @rdname condBiCopSim
#' @export
condBiCopSim <- function(family, calib.fnc, X,
                         par2 = 0, return.par = TRUE, tau = TRUE) {
  if (!is.numeric(family) || length(family) != 1) {
    stop("Input for family has to be a scalar/integer.")
  }

  if (!(family %in% get.familyset())) {
    stop("Copula family not yet implemented.")
  }

  if (!is.function(calib.fnc)) {
    stop("The calibration function should be a function.")
  }

  if (!is.numeric(X)) {
    stop("The covariates should be reals.")
  }

  if (!is.numeric(par2) && length(par2) != 1) {
    stop("par2 should be a real of length 1.")
  }

  if (!is.logical(return.par) || is.na(return.par)) {
    stop("return.par should be TRUE or FALSE.")
  }

  if (!is.logical(tau) || is.na(tau)) {
    stop("tau should be TRUE or FALSE.")
  }

  ## univariate covariate:
  if (is.vector(X)) {
    dim.x <- 1
    dim.arg <- length(formals(calib.fnc))

    if (dim.x != dim.arg) {
      stop(paste("Calibration function requires", dim.arg, "covariates."))
    }

    eta <- sapply(X, function(x) calib.fnc(x))
  }

  ## multivariate covariate:
  if (is.matrix(X) || is.data.frame(X)) {
    dim.x <- ncol(X)
    dim.arg <- length(formals(calib.fnc))
    if (dim.x != dim.arg) {
      stop(paste("Calibration function requires", dim.arg, "covariates."))
    }

    fnc.expr <- paste("calib.fnc(", paste(paste("x[", 1:dim.arg, "]", sep = ""),
      collapse = ","
    ), ")")
    fnc.expr <- parse(text = fnc.expr)
    eta <- apply(X, 1, function(x) eval(fnc.expr))
  }

  if (tau != TRUE) {
    par <- eta2par(eta, family)
  } else {
    par <- tau2par(eta2par(eta, 1), family)
  }

  if (family %in% c(1, 2, 5)) {
    family <- rep(family, length(par))
  } else {
    if (family %in% c(301:304)) {
      fam <- rev(expand.grid(c(23, 33), c(3, 13)))[family - 300, ]
    } else {
      fam <- rev(expand.grid(c(24, 34), c(4, 14)))[family - 400, ]
    }
    sel <- par > 0
    family <- rep(0, length(par))
    family[sel] <- fam[1]
    family[!sel] <- fam[2]
  }

  data <- t(mapply(BiCopSim, 1, family, par, par2))
  if (return.par == TRUE) {
    return(list(data = data, par = par, eta = eta))
  } else {
    return(data)
  }
}
