#' Sequential pair-copula selection and maximum penalized likelihood estimation
#' of a GAM-Vine model.
#'
#' This function select the copula family and estimates the parameter(s) of a
#' Generalized Additive model
#' (GAM) Vine model, where GAMs for individual edges are specified either for
#' the copula parameter or Kendall's tau.
#' It solves the maximum penalized likelihood estimation for the copula families
#' supported in this package by reformulating each Newton-Raphson iteration as
#' a generalized ridge regression, which is solved using
#' the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param data A matrix or data frame containing the data in [0,1]^d.
#' @param Matrix Lower triangular \code{d x d} matrix that defines the R-vine
#' tree structure.
#' @param lin.covs A matrix or data frame containing the parametric (i.e.,
#' linear) covariates (default: \code{lin.covs = NULL}).
#' @param smooth.covs A matrix or data frame containing the non-parametric
#' (i.e., smooth) covariates (default: \code{smooth.covs = NULL}).
#' @param simplified If \code{TRUE}, then a simplified vine is fitted (which is
#' possible only if there are exogenous covariates). If \code{FALSE} (default),
#' then a non-simplified vine is fitted.
#' @param familyset An integer vector of pair-copula families to select from
#' (the independence copula MUST NOT be specified in this vector unless one
#' wants to fit an independence vine!). The vector has to include at least one
#' pair-copula family that allows for positive and one that allows for negative
#' dependence. Not listed copula families might be included to better handle
#' limit cases. If \code{familyset = NA} (default), selection among all
#' possible families is performed. Coding of pair-copula families:
#' \code{1} Gaussian, \code{2} Student t,
#' \code{3} Clayton, \code{4} Gumbel, \code{13} Survival Clayton,
#' \code{14} Survival Gumbel,  \code{23} Rotated (90 degrees) Clayton,
#' \code{24} Rotated (90 degrees) Gumbel,
#' \code{33} Rotated (270 degrees) Clayton and
#' \code{34} Rotated (270 degrees) Gumbel.
#' @param rotations If \code{TRUE}, all rotations of the families in familyset
#' are included.
#' @param familycrit Character indicating the criterion for bivariate copula
#' selection. Possible choices: \code{familycrit = 'AIC'} (default) or
#' \code{'BIC'}, as in \code{\link{BiCopSelect}} from the
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package.
#' @param level Numerical; Passed to \code{\link{gamBiCopSelect}}, it is the
#' significance level of the test for removing individual
#' predictors (default: \code{level = 0.05}) for each conditional pair-copula.
#' @param trunclevel Integer; level of truncation.
#' @param tau \code{TRUE} (default) for a calibration function specified for
#' Kendall's tau or \code{FALSE} for a calibration function specified
#' for the Copula parameter.
#' @param method \code{'NR'} for Newton-Raphson
#' and  \code{'FS'} for Fisher-scoring (default).
#' @param tol.rel Relative tolerance for \code{'FS'}/\code{'NR'} algorithm.
#' @param n.iters Maximal number of iterations for
#' \code{'FS'}/\code{'NR'} algorithm.
#' @param parallel \code{TRUE} (default) for parallel selection of copula
#' family at each edge or \code{FALSE} for the sequential version.
#' for the Copula parameter.
#' @param verbose \code{TRUE} if informations should be printed during the
#' estimation and \code{FALSE} (default) for a silent version.
#' from \code{\link[mgcv:mgcv-package]{mgcv}}.
#' @param select.once if \code{TRUE} the GAM structure is only selected once,
#'   for the family that appears first in \code{familyset}.
#' @param ... Additional parameters to be passed to \code{\link{gam}}
#' from \code{\link[mgcv:mgcv-package]{mgcv}}.
#' @return \code{gamVineCopSelect} returns a \code{\link{gamVine-class}} object.
#' @examples
#' require(mgcv)
#' set.seed(0)
#'
#' ##  Simulation parameters
#' # Sample size
#' n <- 1e3
#' # Copula families
#' familyset <- c(1:2, 301:304, 401:404)
#' # Define a 4-dimensional R-vine tree structure matrix
#' d <- 4
#' Matrix <- c(2, 3, 4, 1, 0, 3, 4, 1, 0, 0, 4, 1, 0, 0, 0, 1)
#' Matrix <- matrix(Matrix, d, d)
#' nnames <- paste("X", 1:d, sep = "")
#'
#' ## A function factory
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
#' ##  Create the model
#' # Define gam-vine model list
#' count <- 1
#' model <- vector(mode = "list", length = d * (d - 1) / 2)
#' sel <- seq(d, d^2 - d, by = d)
#'
#' # First tree
#' for (i in 1:(d - 1)) {
#'   # Select a copula family
#'   family <- sample(familyset, 1)
#'   model[[count]]$family <- family
#'
#'   # Use the canonical link and a randomly generated parameter
#'   if (is.element(family, c(1, 2))) {
#'     model[[count]]$par <- tanh(rnorm(1) / 2)
#'     if (family == 2) {
#'       model[[count]]$par2 <- 2 + exp(rnorm(1))
#'     }
#'   } else {
#'     if (is.element(family, c(401:404))) {
#'       rr <- rnorm(1)
#'       model[[count]]$par <- sign(rr) * (1 + abs(rr))
#'     } else {
#'       model[[count]]$par <- rnorm(1)
#'     }
#'     model[[count]]$par2 <- 0
#'   }
#'   count <- count + 1
#' }
#'
#' # A dummy dataset
#' data <- data.frame(u1 = runif(1e2), u2 = runif(1e2), matrix(runif(1e2 * d), 1e2, d))
#'
#' # Trees 2 to (d-1)
#' for (j in 2:(d - 1)) {
#'   for (i in 1:(d - j)) {
#'     # Select a copula family
#'     family <- sample(familyset, 1)
#'
#'     # Select the conditiong set and create a model formula
#'     cond <- nnames[sort(Matrix[(d - j + 2):d, i])]
#'     tmpform <- paste("~", paste(paste("s(", cond, ", k=10, bs='cr')",
#'       sep = ""
#'     ), collapse = " + "))
#'     l <- length(cond)
#'     temp <- sample(3, l, replace = TRUE)
#'
#'     # Spline approximation of the true function
#'     m <- 1e2
#'     x <- matrix(seq(0, 1, length.out = m), nrow = m, ncol = 1)
#'     if (l != 1) {
#'       tmp.fct <- paste("function(x){eta0+",
#'         paste(sapply(1:l, function(x)
#'           paste("calib.surf[[", temp[x], "]](x[", x, "])",
#'             sep = ""
#'           )), collapse = "+"), "}",
#'         sep = ""
#'       )
#'       tmp.fct <- eval(parse(text = tmp.fct))
#'       x <- eval(parse(text = paste0("expand.grid(",
#'         paste0(rep("x", l), collapse = ","), ")",
#'         collapse = ""
#'       )))
#'       y <- apply(x, 1, tmp.fct)
#'     } else {
#'       tmp.fct <- function(x) eta0 + calib.surf[[temp]](x)
#'       colnames(x) <- cond
#'       y <- tmp.fct(x)
#'     }
#'
#'     # Estimate the gam model
#'     form <- as.formula(paste0("y", tmpform))
#'     dd <- data.frame(y, x)
#'     names(dd) <- c("y", cond)
#'     b <- gam(form, data = dd)
#'     # plot(x[,1],(y-fitted(b))/y)
#'
#'     # Create a dummy gamBiCop object
#'     tmp <- gamBiCopFit(data = data, formula = form, family = 1, n.iters = 1)$res
#'
#'     # Update the copula family and the model coefficients
#'     attr(tmp, "model")$coefficients <- coefficients(b)
#'     attr(tmp, "model")$smooth <- b$smooth
#'     attr(tmp, "family") <- family
#'     if (family == 2) {
#'       attr(tmp, "par2") <- 2 + exp(rnorm(1))
#'     }
#'     model[[count]] <- tmp
#'     count <- count + 1
#'   }
#' }
#'
#' # Create the gamVineCopula object
#' GVC <- gamVine(Matrix = Matrix, model = model, names = nnames)
#' print(GVC)
#' \dontrun{
#' ## Simulate and fit the model
#' sim <- gamVineSimulate(n, GVC)
#' fitGVC <- gamVineSeqFit(sim, GVC, verbose = TRUE)
#' fitGVC2 <- gamVineCopSelect(sim, Matrix, verbose = TRUE)
#'
#' ## Plot the results
#' par(mfrow = c(3, 4))
#' plot(GVC, ylim = c(-2.5, 2.5))
#'
#' plot(fitGVC, ylim = c(-2.5, 2.5))
#'
#' plot(fitGVC2, ylim = c(-2.5, 2.5))
#' }
#'
#' @seealso  \code{\link{gamVineSeqFit}},\code{\link{gamVineStructureSelect}},
#'  \code{\link{gamVine-class}}, \code{\link{gamVineSimulate}} and
#'  \code{\link{gamBiCopFit}}.
gamVineCopSelect <- function(data, Matrix,
                             lin.covs = NULL, smooth.covs = NULL,
                             simplified = FALSE,
                             familyset = NA, rotations = TRUE,
                             familycrit = "AIC", level = 0.05,
                             trunclevel = NA, tau = TRUE, method = "FS",
                             tol.rel = 0.001, n.iters = 10,
                             parallel = FALSE, verbose = FALSE,
                             select.once = TRUE, ...) {
  tmp <- valid.gamVineCopSelect(
    data, Matrix, lin.covs, smooth.covs, simplified,
    familyset, rotations, familycrit, level,
    trunclevel, tau, method, tol.rel, n.iters,
    parallel, verbose, select.once
  )
  if (tmp != TRUE) {
    stop(tmp)
  }

  ## Transform to dataframe, get dimensions, etc (see in utilsPrivate)
  tmp <- prepare.data2(
    data, lin.covs, smooth.covs,
    trunclevel, familyset, rotations
  )
  n <- tmp$n
  d <- tmp$d
  l <- tmp$l
  nn <- tmp$nn
  data <- tmp$data
  covariates <- tmp$covariates
  trunclevel <- tmp$trunclevel
  familyset <- tmp$familyset

  oldMat <- Matrix
  o <- diag(oldMat)
  oo <- o[length(o):1]
  Mat <- reorderRVineMatrix(Matrix)
  data[, 1:d] <- data[, oo]

  MaxMat <- createMaxMat(Mat)
  CondDistr <- neededCondDistr(Mat)

  V <- list()
  V$direct <- array(NA, dim = c(d, d, n))
  V$indirect <- array(NA, dim = c(d, d, n))
  V$direct[d, , ] <- t(data[, d:1])

  model.count <- get.modelCount(d)
  model <- vector("list", d * (d - 1) / 2)
  for (i in (d - 1):1) {
    for (k in d:(i + 1)) {
      mki <- model.count[k, i]
      m <- MaxMat[k, i]
      zr1 <- V$direct[k, i, ]

      if (m == Mat[k, i]) {
        zr2 <- V$direct[k, (d - m + 1), ]
      } else {
        zr2 <- V$indirect[k, (d - m + 1), ]
      }

      if (verbose == TRUE) {
        if (k == d) {
          message(oldMat[i, i], ",", oldMat[k, i])
        } else {
          message(
            oldMat[i, i], ",", oldMat[k, i], "|",
            paste(oldMat[(k + 1):d, i], collapse = ",")
          )
        }
      }

      if (d + 1 - k > trunclevel) {
        tmp <- fitACopula(zr2, zr1, 0, familycrit, level)
      } else {
        if (k == d && l == 0) {
          tmp <- fitACopula(zr2, zr1, familyset, familycrit, level, FALSE)
        } else {
          tmp <- list(
            cbind(
              as.numeric(zr2),
              as.numeric(zr1)
            ),
            lin.covs, smooth.covs
          )
          if (k != d && simplified == FALSE) {
            cond <- Mat[(k + 1):d, i]
            if (is.null(smooth.covs)) {
              tmp[[3]] <- as.data.frame(data[, cond])
              names(tmp[[3]]) <- nn[oo[cond]]
            } else {
              tmp[[3]] <- as.data.frame(do.call(
                cbind,
                list(data[, cond], smooth.covs)
              ))
              names(tmp[[3]]) <- c(nn[oo[cond]], colnames(smooth.covs))
            }
          }

          tmp <- fitAGAMCopula(
            tmp, familyset, familycrit, level,
            tau, method, tol.rel, n.iters, FALSE,
            rotations, select.once, ...
          )
        }
      }
      model[[mki]] <- tmp$model
      if (CondDistr$direct[k - 1, i]) {
        V$direct[k - 1, i, ] <- as.numeric(tmp$CondOn2)
      }

      if (CondDistr$indirect[k - 1, i]) {
        V$indirect[k - 1, i, ] <- as.numeric(tmp$CondOn1)
      }
    }
  }
  return(gamVine(Matrix, model, nn[1:d], covariates))
}

valid.gamVineCopSelect <- function(data, Matrix, lin.covs, smooth.covs,
                                   simplified, familyset, rotations, familycrit,
                                   level, trunclevel, tau, method, tol.rel,
                                   n.iters, parallel, verbose, select.once) {
  tmp <- valid.gamVineStructureSelect(
    data, lin.covs, smooth.covs, simplified,
    0, familyset, rotations, familycrit,
    "tau", level, trunclevel,
    tau, method, tol.rel, n.iters,
    parallel, verbose, select.once
  )
  if (tmp != TRUE) {
    return(tmp)
  }

  d <- dim(data)[2]
  if (!is.matrix(Matrix) || dim(Matrix)[1] != dim(Matrix)[2] ||
    max(Matrix) > d || RVineMatrixCheck(Matrix) != 1) {
    return("Matrix is not a valid R-vine matrix or its dimension is incorrect.")
  }

  return(TRUE)
}
