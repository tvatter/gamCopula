#' Selection and Maximum penalized likelihood estimation of a Generalized
#' Additive model (gam) for the copula parameter or Kendall's tau.
#'
#' This function selects an appropriate bivariate copula family for given
#' bivariate copula data using one of a range of methods. The corresponding
#' parameter estimates are obtained by maximum penalized  likelihood estimation,
#' where each Newton-Raphson iteration is reformulated as a generalized ridge
#' regression solved using the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param udata A matrix or data frame containing the model responses, (u1,u2) in
#' [0,1]x[0,1]
#' @param lin.covs A matrix or data frame containing the parametric (i.e.,
#' linear) covariates.
#' @param smooth.covs A matrix or data frame containing the non-parametric
#' (i.e., smooth) covariates.
#' @param familyset (Similar to \code{\link{BiCopSelect}} from the
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package)
#' Vector of bivariate copula families to select from.
#' If \code{familyset = NA} (default), selection
#' among all possible families is performed.   Coding of bivariate copula
#' families:
#' \code{1} Gaussian,
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
#' @param rotations If \code{TRUE}, all rotations of the families in familyset
#' are included.
#' @param familycrit Character indicating the criterion for bivariate copula
#' selection. Possible choices: \code{familycrit = 'AIC'} (default) or
#' \code{'BIC'}, as in \code{\link{BiCopSelect}} from the
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package.
#' @param level Numerical; significance level of the test for removing individual
#' predictors (default: \code{level = 0.05}).
#' @param edf Numerical; if the estimated EDF for individual predictors is
#' smaller than \code{edf} but the predictor is still significant, then
#' it is set as linear (default: \code{edf = 1.5}).
#' @param tau \code{FALSE} for a calibration function specified for
#' the Copula parameter or \code{TRUE} (default) for a calibration function
#' specified for Kendall's tau.
#' @param method \code{'FS'} for Fisher-scoring (default) and
#' \code{'NR'} for Newton-Raphson.
#' @param tol.rel Relative tolerance for \code{'FS'}/\code{'NR'} algorithm.
#' @param n.iters Maximal number of iterations for
#' \code{'FS'}/\code{'NR'} algorithm.
#' @param parallel \code{TRUE} for a parallel estimation across copula families.
#' @param verbose \code{TRUE} prints informations during the estimation.
#' @param select.once if \code{TRUE} the GAM structure is only selected once,
#'   for the family that appears first in \code{familyset}.
#' @param ... Additional parameters to be passed to \code{\link{gam}}
#' @return \code{gamBiCopFit} returns a list consisting of
#' \item{res}{S4 \code{\link{gamBiCop-class}} object.}
#' \item{method}{\code{'FS'} for Fisher-scoring and
#' \code{'NR'} for Newton-Raphson.}
#' \item{tol.rel}{relative tolerance for \code{'FS'}/\code{'NR'} algorithm.}
#' \item{n.iters}{maximal number of iterations for
#' \code{'FS'}/\code{'NR'} algorithm.}
#' \item{trace}{the estimation procedure's trace.}
#' \item{conv}{\code{0} if the algorithm converged and \code{1} otherwise.}
#' @seealso \code{\link{gamBiCop}} and \code{\link{gamBiCopFit}}.
#' @examples
#' require(copula)
#' set.seed(0)
#'
#' ## Simulation parameters (sample size, correlation between covariates,
#' ## Student copula with 4 degrees of freedom)
#' n <- 5e2
#' rho <- 0.9
#' fam <- 2
#' par2 <- 4
#'
#' ## A calibration surface depending on four variables
#' eta0 <- 1
#' calib.surf <- list(
#'   calib.lin <- function(t, Ti = 0, Tf = 1, b = 2) {
#'     return(-2 + 4 * t)
#'   },
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
#' ## 6-dimensional matrix X of covariates
#' covariates.distr <- mvdc(normalCopula(rho, dim = 6),
#'   c("unif"), list(list(min = 0, max = 1)),
#'   marginsIdentical = TRUE
#' )
#' X <- rMvdc(n, covariates.distr)
#' colnames(X) <- paste("x", 1:6, sep = "")
#'
#' ## U in [0,1]x[0,1] depending on the four first columns of X
#' U <- condBiCopSim(fam, function(x1, x2, x3, x4) {
#'   eta0 + sum(mapply(function(f, x)
#'     f(x), calib.surf, c(x1, x2, x3, x4)))
#' }, X[, 1:4], par2 = 4, return.par = TRUE)
#' \dontrun{
#' ## Selection using AIC (about 30sec on single core)
#' ## Use parallel = TRUE to speed-up....
#' system.time(best <- gamBiCopSelect(U$data, smooth.covs = X))
#' print(best$res)
#' EDF(best$res) ## The first function is linear
#' ## Plot only the smooth component
#' par(mfrow = c(2, 2))
#' plot(best$res)
#' }
#' @export
gamBiCopSelect <- function(udata, lin.covs = NULL, smooth.covs = NULL,
                           familyset = NA, rotations = TRUE,
                           familycrit = "AIC", level = 5e-2, edf = 1.5, tau = TRUE,
                           method = "FS", tol.rel = 1e-3, n.iters = 10,
                           parallel = FALSE, verbose = FALSE,
                           select.once = TRUE, ...) {
  tmp <- valid.gamBiCopSelect(
    udata, lin.covs, smooth.covs, rotations, familyset,
    familycrit, level, edf, tau, method, tol.rel,
    n.iters, parallel, verbose, select.once
  )
  if (tmp != TRUE) {
    stop(tmp)
  }

  if (length(familyset) == 1 && is.na(familyset)) {
    familyset <- get.familyset()
  }
  if (rotations) {
    familyset <- withRotations(familyset)
  }
  familyset <- unique(familyset)

  parallel <- as.logical(parallel)
  if (parallel) {
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl, cores = detectCores() - 1)
    on.exit(stopCluster(cl))
  }

  x <- NULL
  if (select.once) {
    res <- vector("list", length(familyset))
    # select GAM structure for first family in familyset
    res[[1]] <- tryCatch(gamBiCopVarSel(
      udata, lin.covs, smooth.covs,
      familyset[1], tau, method, tol.rel,
      n.iters, level, edf, verbose, ...
    ),
    error = function(e) e
    )
    # fit with same structure for all other families
    if (length(familyset) > 1) {
      # combine data and covariates into one matrix
      tmpdata <- udata
      if (!is.null(lin.covs)) {
        tmpdata <- cbind(tmpdata, lin.covs)
      }
      if (!is.null(smooth.covs)) {
        tmpdata <- cbind(tmpdata, smooth.covs)
      }

      if (!parallel) {
        res[-1] <- foreach(x = familyset[-1]) %do%
          tryCatch(gamBiCopFit(
            tmpdata,
            res[[1]]$res@model$formula,
            x, tau, method, tol.rel, n.iters,
            verbose, ...
          ),
          error = function(e) e
          )
      } else {
        res[-1] <- foreach(x = familyset[-1]) %dopar%
          tryCatch(gamBiCopFit(
            tmpdata,
            res[[1]]$res@model$formula,
            x, tau, method, tol.rel, n.iters,
            verbose, ...
          ),
          error = function(e) e
          )
      }
    }
  } else {
    # select GAM structure independently for all families
    if (!parallel) {
      res <- foreach(x = familyset) %do%
        tryCatch(gamBiCopVarSel(
          udata, lin.covs, smooth.covs,
          x, tau, method, tol.rel, n.iters, # max.knots,
          level, edf, verbose, ...
        ),
        error = function(e) e
        )
    } else {
      res <- foreach(x = familyset) %dopar%
        tryCatch(gamBiCopVarSel(
          udata, lin.covs, smooth.covs,
          x, tau, method, tol.rel, n.iters,
          level, edf, verbose, ...
        ),
        error = function(e) e
        )
    }
  }

  res <- res[sapply(res, length) == 6]
  if (length(res) == 0) { #|| sum(sapply(res, function(x) x$conv) == 0) == 0) {
    return(paste(
      "No convergence of the estimation for any copula family.",
      "Try modifying the parameters (FS/NR, tol.rel, n.iters)..."
    ))
  }

  if (familycrit == "AIC") {
    return(res[[which.min(sapply(res, function(x) AIC(x$res)))]])
  } else {
    return(res[[which.min(sapply(res, function(x) BIC(x$res)))]])
  }
}

gamBiCopVarSel <- function(udata, lin.covs, smooth.covs,
                           family, tau = TRUE, method = "FS",
                           tol.rel = 1e-3, n.iters = 10, # max.knots,
                           level = 5e-2, edf = 1.5,
                           verbose = FALSE, ...) {
  udata <- as.data.frame(udata)
  colnames(udata) <- c("u1", "u2")
  data <- udata

  if (!is.null(lin.covs)) {
    if (is.null(colnames(lin.covs))) {
      lin.covs <- as.data.frame(lin.covs)
      colnames(lin.covs) <- paste0("l", 1:ncol(lin.covs))
    }
    data <- cbind(data, lin.covs)
  }
  if (!is.null(smooth.covs)) {
    if (is.null(colnames(smooth.covs))) {
      smooth.covs <- as.data.frame(smooth.covs)
      colnames(smooth.covs) <- paste0("s", 1:ncol(smooth.covs))
    }
    data <- cbind(data, smooth.covs)
  }

  if (verbose == TRUE) {
    cat(paste("Model selection for family", family, "\n"))
  }

  myGamBiCopEst <- function(formula) gamBiCopFit(
      data, formula, family, tau,
      method, tol.rel, n.iters
    )

  n <- nrow(data)
  d <- ncol(data) - 2
  ## Create a list with formulas of smooth terms corresponding to the covariates
  lin.nms <- colnames(lin.covs)

  nn <- smooth.nms <- colnames(smooth.covs)

  if (!is.null(ncol(smooth.covs))) {
    # if (length(max.knots) == 1) {
    #   basis <- rep(max.knots, ncol(smooth.covs))
    # } else {
    #   stopifnot(length(max.knots) == ncol(smooth.covs))
    #   basis <- max.knots
    # }
    basis <- rep(10, ncol(smooth.covs))
    formula.expr <- mapply(get.smooth, smooth.nms, basis)
  } else {
    basis <- NULL
    formula.expr <- NULL
  }
  if (!is.null(ncol(lin.covs))) {
    formula.lin <- lin.nms # get.linear(lin.nms, lin.nms)
  } else {
    formula.lin <- NULL
  }

  ## Update the list by removing unsignificant predictors
  if (verbose == TRUE) {
    cat("Remove unsignificant covariates.......\n")
    t1 <- Sys.time()
  }

  sel.smooth <- smooth2lin <- FALSE
  sel.lin <- NULL

  while ((!all(c(sel.lin, sel.smooth)) && length(basis) > 0) || any(smooth2lin)) {
    smooth2lin <- FALSE
    sel.lin <- NULL
    formula.tmp <- get.formula(c(formula.lin, formula.expr))
    # if (verbose == TRUE) {
    #   cat("Model formula:\n")
    #   print(formula.tmp)
    # }

    tmp <- myGamBiCopEst(formula.tmp)

    ## Remove unsignificant parametric terms
    if (length(summary(tmp$res@model)$p.coeff) > 1) {
      sel.lin <- summary(tmp$res@model)$p.pv[-1] < level
      formula.lin <- formula.lin[sel.lin]
    }

    if (summary(tmp$res@model)$m > 0) {
      ## Remove unsignificant smooth terms
      sel.smooth <- summary(tmp$res@model)$s.pv < level
      nn <- nn[sel.smooth]
      basis <- rep(10, length(nn))
      formula.expr <- mapply(get.smooth, nn, basis)

      ## Check whether some smooth need to be set to linear
      ## (and adjust the model if required)
      smooth2lin <- sel.smooth & summary(tmp$res@model)$edf < edf
      if (any(smooth2lin)) {
        sel.lin <- which(sel.smooth)
        sel.lin <- summary(tmp$res@model)$edf[sel.lin] < edf
        sel.smooth <- !sel.lin
        formula.lin <- c(formula.lin, nn[sel.lin])
        formula.expr <- formula.expr[sel.smooth]
        basis <- basis[sel.smooth]
        nn <- nn[sel.smooth]
      }
    }
  }

  ## Check if we can output a constant or linear model directly
  ## (or re-estimate the model if not)
  if (length(formula.expr) == 0) {
    if ((is.null(formula.lin) ||
      length(formula.lin) == 0)) {
      tmp <- myGamBiCopEst(~1)
    } else {
      formula.tmp <- get.formula(formula.lin)
      tmp <- myGamBiCopEst(formula.tmp)
    }
    return(tmp)
  } else {
    formula.tmp <- get.formula(c(formula.lin, formula.expr))
    tmp <- myGamBiCopEst(formula.tmp)
  }

  ## Increasing the basis size appropriately
  sel <- summary(tmp$res@model)$edf > (basis - 1) * 0.8
  if (verbose == TRUE && any(sel)) {
    t2 <- Sys.time()
    cat(paste0("Time elapsed: ", round(as.numeric(t2 - t1), 2), " seconds.\n"))
    t1 <- Sys.time()
    cat(paste("Select the basis sizes .......\n"))
  }

  while (any(sel) && all(basis < n / 30)) {

    ## Increase basis sizes
    basis[sel] <- 2 * basis[sel]

    ## Update the smooth terms, formula and fit the new model
    if (any(sel)) {
      formula.expr[sel] <- mapply(get.smooth, nn[sel], basis[sel])
      formula.tmp <- get.formula(c(formula.lin, formula.expr))
      # if (verbose == TRUE) {
      #   cat("Updated model formula:\n")
      #   print(formula.tmp)
      # }
      tmp <- myGamBiCopEst(formula.tmp)
      sel[sel] <- summary(tmp$res@model)$edf[sel] > (basis[sel] - 1) * 0.8
    }

    ## Check that the basis size is smaller than 1/2 the size of
    ## the corresponding covariate's support
    if (any(sel)) {
      if (sum(sel) == 1) {
        sel[sel] <- 2 * basis[sel] < length(unique(data[, nn[sel]]))
      } else {
        sel[sel] <- 2 * basis[sel] < apply(data[, nn[sel]], 2, function(x)
          length(unique(x)))
      }
    }
  }

  if (verbose == TRUE) {
    t2 <- Sys.time()
    cat(paste0("Time elapsed: ", round(as.numeric(t2 - t1), 2), " seconds.\n"))
  }
  return(tmp)
}

get.smooth <- function(x, k, bs = "cr") {
  paste("s(", x, ", k=", k, ", bs='", bs, "')", sep = "")
}
get.formula <- function(expr, y = FALSE) {
  if (!y) {
    as.formula(paste("~", paste(expr, collapse = " + ")))
  } else {
    as.formula(paste("y ~", paste(expr, collapse = " + ")))
  }
}

get.linear <- function(x, nn) {
  names(unlist(sapply(nn, function(z) grep(z, x))))
}

valid.gamBiCopSelect <- function(udata, lin.covs, smooth.covs, rotations,
                                 familyset, familycrit, level, edf,
                                 tau, method, tol.rel, n.iters, parallel,
                                 verbose, select.once) {
  data <- tryCatch(as.data.frame(udata), error = function(e) e$message)
  if (is.character(udata)) {
    return(udata)
  }

  tmp <- valid.gamBiCopFit(udata, n.iters, tau, tol.rel, method, verbose, 1)
  if (tmp != TRUE) {
    return(tmp)
  }
  m <- dim(udata)[2]
  n <- dim(udata)[1]
  if (m != 2) {
    return("udata should have only two columns.")
  }

  if (!is.null(lin.covs) && !(is.data.frame(lin.covs) || is.matrix(lin.covs))) {
    return("lin.covs has to be either a data frame or a matrix.")
  } else {
    if (!is.null(lin.covs) && dim(lin.covs)[1] != n) {
      return("lin.covs must have the same number of rows as udata.")
    }
  }
  if (!is.null(smooth.covs) && !(is.data.frame(smooth.covs) ||
    is.matrix(smooth.covs))) {
    return("smooth.covs has to be either a data frame or a matrix.")
  } else {
    if (!is.null(smooth.covs) && dim(smooth.covs)[1] != n) {
      return("smooth.covs must have the same number of rows as udata.")
    }
  }

  if (!valid.familyset(familyset)) {
    return(return(msg.familyset(var2char(familyset))))
  }

  if (!valid.logical(rotations)) {
    return(msg.logical(var2char(rotations)))
  }

  if (is.null(familycrit) || length(familycrit) != 1 ||
    (familycrit != "AIC" && familycrit != "BIC")) {
    return("Selection criterion not implemented.")
  }

  if (!valid.logical(parallel)) {
    return(msg.logical(var2char(parallel)))
  }

  if (!valid.unif(level)) {
    return(msg.unif(var2char(level)))
  }

  if (!valid.real(edf)) {
    return(msg.real(var2char(edf)))
  }

  if (!is.logical(select.once)) {
    return("'select.once' must be logical")
  }

  return(TRUE)
}
