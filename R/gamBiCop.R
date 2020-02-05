#' Construction of a gamBiCop Class Object
#'
#' Constructs an object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#'
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
#' @param model A \code{\link{gamObject}} as return by the
#' \code{\link{gam}} function
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#' @param par2 Second parameter for the Student t-copula.
#' @param tau \code{FALSE} for a calibration function specified
#' for the Copula parameter or \code{TRUE} (default) for a calibration
#' function specified for Kendall's tau.
#' @return An object of the class
#'  \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @seealso \code{\link[gamCopula:gamBiCop-class]{gamBiCop}},
#' \code{\link{gamBiCopFit}}, \code{\link{gamBiCopPredict}} and
#' \code{\link{gamBiCopSimulate}}.
#' @name gamBiCop
#' @rdname gamBiCop
#' @export
gamBiCop <- function(family, model, par2 = 0, tau = TRUE) {
  if (family != 2) {
    par2 <- 0
  }
  tmp <- tryCatch(as.integer(family), error = function(e) e)
  msg <- "should be or be coercisable to an integer."
  if (!is(tmp, "integer") ||
    length(family) != 1 || as.integer(family) != family) {
    stop(paste("family", msg))
  }
  new("gamBiCop",
    family = as.integer(family), model = model,
    par2 = par2, tau = tau
  )
}

valid.gamBiCop <- function(object) {
  d <- length(attributes(object))
  if (d < 2) {
    return("A gamBiCop contains at least a copula family and a mgcv model.")
  } else if (d >= 2) {
    if (!is(object@model, "gam")) {
      return("Invalid mgcv model.")
    } else if (!(object@family %in% get.familyset())) {
      return("Copula family not yet implemented.")
    }
  }
  if (!(is.logical(object@tau) || (object@tau == 0) || (object@tau == 1))) {
    return(paste(
      "Tau should takes 0/1 or FALSE/TRUE to model",
      "the copula parameter/tau's tau."
    ))
  }
  if ((object@family == 2) &&
    (is.null(object@par2) || is.na(as.numeric(object@par2)) ||
      as.numeric(object@par2) <= 2)) {
    return("Par2 greater than 2 is needed for the t-copula.")
  }
  return(TRUE)
}

show.gamBiCop <- function(object) {
  cat(bicopname(object@family), "copula with ")
  if (object@tau == TRUE) {
    cat("tau(z) = (exp(z)-1)/(exp(z)+1) where \n")
  } else {
    if (object@family %in% c(1, 2)) {
      cat("par(z) = (exp(z)-1)/(exp(z)+1) where \n")
    } else if (object@family %in% c(3, 13)) {
      cat("par(z) = exp(z) where \n")
    } else if (object@family %in% c(4, 14)) {
      cat("par(z) = 1+exp(z) where \n")
    } else if (object@family %in% c(23, 33)) {
      cat("par(z) = -exp(z) where \n")
    } else if (object@family %in% c(24, 34)) {
      cat("par(z) = -1-exp(z) where \n")
    }
  }
  show(object@model$formula)
}

summary.gamBiCop <- function(object, ...) {
  cat(bicopname(object@family), "copula with ")
  if (object@tau == TRUE) {
    cat("tau(z) = (exp(z)-1)/(exp(z)+1) where \n")
  } else {
    if (object@family %in% c(1, 2)) {
      cat("par(z) = (exp(z)-1)/(exp(z)+1) where \n")
    } else if (object@family %in% c(3, 13)) {
      cat("par(z) = exp(z) where \n")
    } else if (object@family %in% c(4, 14)) {
      cat("par(z) = 1+exp(z) where \n")
    } else if (object@family %in% c(23, 33)) {
      cat("par(z) = -exp(z) where \n")
    } else if (object@family %in% c(24, 34)) {
      cat("par(z) = -1-exp(z) where \n")
    }
  }
  tmp <- capture.output(summary(object@model))
  cat(paste(tmp[-c(1:4)], "\n"))
}

nobs.gamBiCop <- function(object, ...) {
  n <- dim(object@model$data)[1]
  return(n)
}

logLik.gamBiCop <- function(object, ...) {
  family <- object@family
  par <- gamBiCopPredict(object, target = "par")$par
  data <- cbind(na.omit(object@model$data)[, c(3, 4)], par)

  if (family == 2) {
    par2 <- rep(object@par2, length(par))
    data <- cbind(data, par2)
  }

  Li <- as.numeric(bicoppd1d2(data, family))
  val <- sum(log(Li))
  df <- sum(object@model$edf)

  if (family == 2) {
    df <- df + 1
  }

  attr(val, "df") <- df
  attr(val, "nobs") <- dim(data)[1]
  class(val) <- "logLik"
  return(val)
}

AIC.gamBiCop <- function(object, ..., k = 2) {
  l <- logLik.gamBiCop(object)
  d <- attributes(l)$df
  return(k * d - 2 * l[1])
}

BIC.gamBiCop <- function(object, ...) {
  return(AIC.gamBiCop(object, ..., k = log(nobs(object))))
}


plot.gamBiCop <- function(x, ...) {
  plot.gam(x@model, ...)
}

formula.gamBiCop <- function(x, ...) {
  return(x@model$formula)
}

setValidity("gamBiCop", valid.gamBiCop)
setMethod("show", signature("gamBiCop"), show.gamBiCop)

#' Summary for a gamBiCop Object
#'
#' Takes a \code{\link[gamCopula:gamBiCop-class]{gamBiCop}} object and produces
#' various useful summaries from it.
#'
#' @param object An object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @param ... unused in this class
#' @return A useful summary (see \code{\link{summary.gam}}
#' from \code{\link[mgcv:mgcv-package]{mgcv}} for more details).
#' @seealso \code{\link{summary.gam}}
#' from \code{\link[mgcv:mgcv-package]{mgcv}}
#' @docType methods
#' @name summary.gamBiCop
#' @rdname summary.gamBiCop-methods
#' @aliases summary,gamBiCop-method
#' @export
setMethod("summary", signature("gamBiCop"), summary.gamBiCop)

#' Model Formula of the gamBiCop Object
#'
#' Extracts the \code{\link{gam}} formula from an object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' This function is a wrapper to \code{\link{formula.gam}}
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param x An object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @param ... un-used in this class
#' @seealso \code{\link{formula.gam}} function
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#' @docType methods
#' @name formula.gamBiCop
#' @rdname formula.gamBiCop-methods
#' @aliases formula,gamBiCop-method
#' @export
setMethod("formula", signature("gamBiCop"), formula.gamBiCop)

#' Extract the Number of Observations from gamBiCop Object
#'
#' Extract the number of 'observations' from a model fit.
#' This is principally intended to be used in computing the BIC
#' (see \code{\link{AIC}}).
#'
#' @param object An object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @param ... un-used in this class
#' @return A single number, normally an integer.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @name nobs.gamBiCop
#' @rdname nobs.gamBiCop-methods
#' @aliases nobs,gamBiCop-method
#' @export
setMethod("nobs", signature("gamBiCop"), nobs.gamBiCop)

#' Extract the Log-likelihood from a gamBiCop Object
#'
#' Function to extract the log-likelihood from an object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}} (note that the models are
#' usually fitted by penalized likelihood maximization).
#' This function is used by \code{\link{AIC}} and \code{\link{BIC}}.
#'
#' @param object An object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @param ... un-used in this class
#' @return Standard \code{logLik} object: see \code{\link{logLik}}.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @name logLik.gamBiCop
#' @rdname logLik.gamBiCop-methods
#' @aliases logLik,gamBiCop-method
#' @export
setMethod("logLik", signature("gamBiCop"), logLik.gamBiCop)

#' Schwarz's Bayesian Information Criterion for a gamBiCop Object
#'
#' Function calculating the Schwarz's Bayesian Information Criterion (BIC)
#' for an object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}} (note that the models are
#' usually fitted by penalized likelihood maximization).
#'
#' @param object An object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @param ... un-used in this class
#' @return A numeric value with the corresponding BIC.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @name BIC.gamBiCop
#' @rdname BIC.gamBiCop-methods
#' @aliases BIC,gamBiCop-method
#' @export
setMethod("BIC", signature("gamBiCop"), BIC.gamBiCop)

#' Akaike's An Information Criterion for a gamBiCop Object
#'
#' Function calculating Akaike's 'An Information Criterion' (AIC) for an object
#' of the class \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}
#' (note that the models are usually fitted by penalized likelihood
#' maximization).
#'
#' @param object An object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @param ... un-used in this class
#' @param k numeric, the penalty per parameter to be used; the default
#' \code{k = 2} is the classical AIC.
#' @return A numeric value with the corresponding AIC.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @name AIC.gamBiCop
#' @rdname AIC.gamBiCop-methods
#' @aliases AIC,gamBiCop-method
#' @export
setMethod("AIC", signature("gamBiCop"), AIC.gamBiCop)

#' Plot a gamBiCop Object
#'
#' Plot from an object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' The function is based on (see \code{\link{plot.gam}}
#' from \code{\link[mgcv:mgcv-package]{mgcv}}).
#'
#' @param x An object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @param y Not used with this class.
#' @param ... additional arguments to be passed to \code{\link{plot.gam}}.
#' @return This function simply generates plots.
#' @seealso \code{\link{plot.gam}} from \code{\link[mgcv:mgcv-package]{mgcv}}).
#' @docType methods
#' @name plot.gamBiCop
#' @rdname plot.gamBiCop-methods
#' @aliases plot,gamBiCop,ANY-method
#' @export
setMethod("plot", signature(x = "gamBiCop"), plot.gamBiCop)


#' Equivalent Degrees of Freedom for an Object of the Class gamBiCop
#'
#' Function calculating the Equivalent Degrees of Freedom (EDF)
#' for a \code{\link[gamCopula:gamBiCop-class]{gamBiCop}} object.
#' It basically sums the edf of the \code{\link{gamObject}}
#' for each smooth component.
#'
#' @param object An object of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @return Estimated degrees of freedom for each smooth component.
#' @docType methods
#' @rdname EDF.gamBiCop-methods
#' @export
EDF <- function(object) {
  edf <- object@model$edf[-1]

  param.terms <- object@model$pterms
  ll.param <- dim(attributes(param.terms)$factors)[2]
  if (is.null(ll.param)) {
    ll.param <- 0
  }

  smooth.terms <- object@model$smooth
  ll.smooth <- length(smooth.terms)
  bs.dim <- unlist(lapply(smooth.terms, function(x) x$bs.dim)) - 1

  out <- rep(NA, ll.param + ll.smooth + 1)
  out[1:(ll.param + 1)] <- 1
  sel <- c(ll.param, ll.param + cumsum(bs.dim))
  if (length(sel) > 1) {
    for (i in 1:(length(sel) - 1)) {
      out[ll.param + 1 + i] <- sum(edf[(sel[i] + 1):sel[i + 1]])
    }
  }

  return(out)
}
