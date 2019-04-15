#' Conditional density function of a gamVine
#'
#' This function returns the density of a conditional pair-copula constructions,
#' where either the copula parameters or the Kendall's taus are modeled as a function
#' of the covariates.
#'
#' @param object \code{\link{gamVine-class}} object.
#' @param data (Same as in \code{\link{predict.gam}} from the
#' \code{\link[mgcv:mgcv-package]{mgcv}} package) A matrix or data frame
#' containing the values of the model covariates at which predictions are
#' required, along with a number of additional columns corresponding to the
#' variables in the pair copula decomposition.
#' @return The conditional density.
#' @seealso \code{\link{gamVine}}
gamVinePDF <- function(object, data) {
  tmp <- valid.gamVine(object)
  if (tmp != "TRUE") {
    return(tmp)
  }
  covariates <- object@covariates

  ## Transform to dataframe, get dimensions, etc (see in utilsPrivate)
  tmp <- prepare.data(data, covariates)
  n <- tmp$n
  d <- tmp$d
  l <- tmp$l
  nn <- tmp$nn
  data <- tmp$data
  covariates <- tmp$covariates

  oldobject <- object
  oldMat <- object@Matrix
  o <- diag(oldMat)
  oo <- o[length(o):1]
  if (any(o != length(o):1)) {
    object <- gamVineNormalize(object)
    data[, 1:d] <- data[, oo]
  }

  Mat <- object@Matrix
  fam <- gamVineFamily(object)
  MaxMat <- createMaxMat(Mat)
  CondDistr <- neededCondDistr(Mat)

  V <- list()
  V$dens <- array(1, dim = c(d, d, n))
  V$direct <- array(NA, dim = c(d, d, n))
  V$indirect <- array(NA, dim = c(d, d, n))
  V$direct[d, , ] <- t(data[, d:1])

  model.count <- get.modelCount(d)
  for (i in (d - 1):1) {
    for (k in d:(i + 1)) {
      # print(model.count[k, i])
      m <- MaxMat[k, i]
      zr1 <- V$direct[k, i, ]

      if (m == Mat[k, i]) {
        zr2 <- V$direct[k, (d - m + 1), ]
      } else {
        zr2 <- V$indirect[k, (d - m + 1), ]
      }

      mki <- model.count[k, i]
      mm <- object@model[[mki]]
      if (valid.gamBiCop(mm) != TRUE) {
        par <- rep(mm$par, n)
        par2 <- mm$par2
      } else {
        par <- gamBiCopPredict(mm, target = "par")$par
        par2 <- mm@par2
      }
      fams <- vapply(
        1:length(par),
        function(j) famTrans(fam[k, i], inv = FALSE, par = par[j]),
        numeric(1)
      )
      V$dens[k, i, ] <- BiCopPDF(zr2, zr1, fams, par, par2, check.pars = FALSE)

      if (CondDistr$direct[k - 1, i]) {
        V$direct[k - 1, i, ] <- BiCopHfunc(zr2, zr1,
          fams, par, par2,
          check.pars = FALSE
        )$hfunc1
      }

      if (CondDistr$indirect[k - 1, i]) {
        V$indirect[k - 1, i, ] <- BiCopHfunc(zr2, zr1,
          fams, par, par2,
          check.pars = FALSE
        )$hfunc2
      }
    }
  }

  c(apply(V$dens, 3, prod))
}

npars.gamVine <- function(object, ...) {
  sum(sapply(object@model, pair_npar))
}

pair_npar <- function(x) {
  if (gamCopula:::valid.gamBiCop(x) != TRUE) {
    return(VineCopula::BiCop(x$family, x$par, x$par2)$npars)
  }
  l <- gamCopula:::logLik.gamBiCop(x)
  attributes(l)$df
}

AIC.gamVine <- function(object, data) {
  -2 * sum(log(gamVinePDF(data, object))) + 2 * npars.gamVine(object)
}
