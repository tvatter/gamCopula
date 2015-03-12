
show.gamBiCop <- function(object) {
  cat("Family: ", object@family, "\n")
  if (object@tau == TRUE) {
    cat("Model: ")
    cat("tau(z) = (exp(z)-1)/(exp(z)+1) where \n")
  } else {
    cat("Model for the Copula parameter:\n")
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
setMethod("show", signature("gamBiCop"), show.gamBiCop)

#' Extract the number of obserations from a fitted \code{\link{gamBiCop-class}} object
#' 
#' Extract the number of 'observations' from a model fit. 
#' This is principally intended to be used in computing BIC (see \code{\link{AIC}}).
#'
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @return A single number, normally an integer.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @name nobs.gamBiCop
#' @rdname nobs-methods
#' @export
nobs.gamBiCop <- function(object, ...) {
  n <- dim(object@model$data)[1]
  return(n)
}
#' @docType methods
#' @rdname nobs-methods
setMethod("nobs", signature("gamBiCop"), nobs.gamBiCop)


#' Log-likelihood for a \code{\link{gamBiCop-class}} object
#' 
#' Function to extract the log-likelihood for a fitted \code{\link{gamBiCop-class}}
#' object (note that the models are usually fitted by penalized likelihood maximization). 
#' This function is used by \code{\link{AIC}} and \code{\link{BIC}}.
#'
#' @param object \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @return Standard \code{logLik} object: see \code{\link{logLik}}.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @rdname logLik-methods
#' @export
logLik.gamBiCop <- function(object, ...) {
  family <- object@family
  par <- gamBiCopPred(object, target = "par")$par
  data <- cbind(na.omit(object@model$data)[, c(3, 4)], par)
  
  if (family == 2) {
    par2 <- rep(object@par2, length(par))
    data <- cbind(data, par2)
    Li <- apply(data, 1, function(x) BiCopPDF(x[1], x[2], family = 2, x[3], x[4]))
  } else {
    Li <- apply(data, 1, function(x) BiCopPDF(x[1], x[2], family = family, x[3]))
  }
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
#' @docType methods
#' @rdname logLik-methods
#' @export
setMethod("logLik", signature("gamBiCop"), logLik.gamBiCop)

#' Akaike's 'An Information Criterion' for a fitted \code{\link{gamBiCop-class}}
#' 
#' Function calculating Akaike's 'An Information Criterion' (AIC) for a fitted \code{\link{gamBiCop-class}}
#' object (note that the models are usually fitted by penalized likelihood maximization). 
#'
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @param k numeric, the penalty per parameter to be used; the default \code{k = 2} is the classical AIC.
#' @return A numeric value with the corresponding AIC.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @rdname AIC-methods
AIC.gamBiCop <- function(object, ..., k = 2) {
  l <- logLik.gamBiCop(object)
  d <- attributes(l)$df
  return(k * d - 2 * l[1])
}
#' @docType methods
#' @rdname AIC-methods
#' @export
setMethod("AIC", signature("gamBiCop"), AIC.gamBiCop)

#' Schwarz's Bayesian Information Criterion for a fitted \code{\link{gamBiCop-class}}
#' 
#' Function calculating the Schwarz's Bayesian Information Criterion (BIC) 
#' for a fitted \code{\link{gamBiCop-class}} object 
#' (note that the models are usually fitted by penalized likelihood maximization). 
#'
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @return A numeric value with the corresponding BIC.
#' @seealso \code{\link{AIC}} and \code{\link{BIC}}.
#' @docType methods
#' @rdname BIC-methods
BIC.gamBiCop <- function(object, ...) {
  return(AIC.gamBiCop(object, ..., k = log(nobs(object))))
}
#' @docType methods
#' @rdname BIC-methods
#' @export
setMethod("BIC", signature("gamBiCop"), BIC.gamBiCop)


#' \code{\link{gamBiCop-class}} formula
#' 
#' Description of the \code{\link{gam}} formula for  fitted \code{\link{gamBiCop-class}} object.
#' This function is a wrapper to \code{\link{formula.gam}}
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param x fitted \code{\link{gamBiCop-class}} object.
#' @param ... un-used in this class
#' @seealso \code{\link{formula.gam}} function 
#' from the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#' @docType methods
#' @rdname formula-methods
formula.gamBiCop <- function(x, ...) {
  return(x@model$formula)
}
#' @docType methods
#' @rdname formula-methods
#' @export
setMethod("formula", signature("gamBiCop"), formula.gamBiCop)

#' Plot a fitted \code{\link{gamBiCop-class}} object
#' 
#' Plot from a model fit. 
#' The function is based on (see \code{\link{plot.gam}}
#' from \code{\link[mgcv:mgcv-package]{mgcv}}).
#' 
#' @param x fitted \code{\link{gamBiCop-class}} object.
#' @param ... additional arguments to be passed to \code{\link{plot.gam}}.
#' @return This function simply generates plots.
#' @seealso \code{\link{plot.gam}} from \code{\link[mgcv:mgcv-package]{mgcv}}).
#' @docType methods
#' @rdname plot-methods
#' @export
plot.gamBiCop <- function(x, ...) {
  plot.gam(x@model, ...)
}
setMethod("plot", signature(x="gamBiCop"), plot.gamBiCop)

#' Equivalent Degrees of Freedom for a fitted \code{\link{gamBiCop-class}}
#' 
#' Function calculating the Equivalent Degrees of Freedom (EDF) 
#' for a fitted \code{\link{gamBiCop-class}} object. 
#' It basically sums the edf of the \code{\link{gamObject}} 
#' for each smooth component.
#'
#' @param object fitted \code{\link{gamBiCop-class}} object.
#' @return Estimated degrees of freedom for each smooth component.
#' @docType methods
#' @rdname EDF-methods
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
