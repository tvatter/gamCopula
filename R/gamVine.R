#' Construction of a gamVine Class Object
#'
#' Constructs an object of the class
#' \code{\link[gamCopula:gamVine-class]{gamVine}}.
#'
#' @param Matrix lower triangular d x d matrix that defines the tree structure.
#' @param model list containing d x (d-1)/2 lists with three numeric items
#' (family, par and par2) and/or objects of the class
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}.
#' @param names vector of d names.
#' @param covariates vector of names for the covariates.
#' @return An object of the class
#'  \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @seealso \code{\link[gamCopula:gamVine-class]{gamVine}},
#' \code{\link{RVineMatrix}}, \code{\link[gamCopula:gamBiCop-class]{gamBiCop}}
#' \code{\link{gamVineSeqFit}}, \code{\link{gamVineCopSelect}},
#' \code{\link{gamVineStructureSelect}} and \code{\link{gamVineSimulate}}.
#' @name gamVine
#' @rdname gamVine
#' @export
gamVine <- function(Matrix, model, names = NA, covariates = NA) {
  tmp <- tryCatch(as.character(names), error = function(e) e)
  msg <- "should be or be coercisable to a character vector."
  if (!is(tmp, "character")) {
    stop(paste("names", msg))
  }
  tmp <- tryCatch(as.character(covariates), error = function(e) e)
  if (!is(tmp, "character")) {
    stop(paste("covariates", msg))
  }
  new("gamVine",
    Matrix = Matrix, model = model,
    names = as.character(names), covariates = as.character(covariates)
  )
}

valid.gamVine <- function(object) {
  d <- length(attributes(object))
  if ((d < 3) || any(names(attributes(object))[1:3] !=
    c("Matrix", "model", "names"))) {
    msg <- paste("A gamVine contains at least 1) a R-Vine matrix, 2) a list of",
      "lists with three numeric items (family, par and par2) and/or",
      " gamBiCop objects and 3) a vector of names.",
      sep = ""
    )
    return(msg)
  }
  Matrix <- object@Matrix
  names <- object@names
  Matrix[is.na(Matrix)] <- 0
  d <- dim(Matrix)[1]
  if (dim(Matrix)[2] != d) {
    return("Structure matrix has to be quadratic.")
  }
  if (max(Matrix) > d) {
    return("Error in the structure matrix.")
  }
  if (RVineMatrixCheck(Matrix) != 1) {
    return("'Matrix' is not a valid R-vine matrix")
  }
  if (length(names) > 1 & length(names) != d) {
    return("Length of the vector 'names' is not correct.")
  } else if (length(names) == 0) {
    names <- paste("x", 1:d, sep = "")
  }

  model <- object@model
  covariates <- object@covariates
  if (length(model) != d * (d - 1) / 2) {
    return("Length of the list 'model' is not correct.")
  }
  count <- 1

  # Trees 1 to (d-1)
  for (j in 1:(d - 1)) {
    for (i in 1:(d - j)) {
      # if (i == 2) {
      #  browser()
      # }
      mm <- model[[count]]
      if (valid.gamBiCop(mm) == TRUE) {
        cond <- sort(unique(all.vars(mm@model$pred.formula)))
        cond2 <- covariates
        if (j != 1) {
          cond2 <- c(cond2, names[sort(Matrix[(d - j + 2):d, i])])
        }

        if (!all(is.element(cond, cond2))) {
          msg <- paste(
            "The formula of element", count, "of model (tree", j, ")",
            " does not not contain the correct conditioning",
            " variables."
          )
          return(msg)
        }
      } else if (is.list(mm)) {
        if (any(!is.element(names(mm), c("family", "par", "par2"))) ||
          !is.numeric(unlist(mm))) {
          msg <- paste(
            "Element", count, "of model, (tree", j, ") should be a",
            " valid gamBiCop object or a list containing three",
            " items (family, par, par2)."
          )
          return(msg)
        }
        chk <- family.check(mm$family, mm$par, mm$par2)
        if (chk != TRUE) {
          return(paste("In element", count, "of model (tree", j, "):", chk))
        }
      } else {
        msg <- paste(
          "Element", count, "of model (tree", j, ") should be a",
          " valid gamBiCop object or a list containing three",
          " items (family, par, par2)."
        )
        return(msg)
      }
      count <- count + 1
    }
  }
  return(TRUE)
}

dim.gamVine <- function(x) {
  GVC <- x
  return(dim(GVC@Matrix)[1])
}

show.gamVine <- function(object) {
  detail <- FALSE
  GVC <- object
  cat("GAM-Vine matrix:", "\n")
  show(GVC@Matrix)
  cat("\n", "Where", "\n")
  for (i in 1:length(GVC@names)) {
    cat(i, " <-> ", GVC@names[[i]], "\n")
  }

  d <- dim(GVC)
  count <- 1
  cat("\n", "Tree 1:", "\n")
  for (i in 1:(d - 1)) {
    a <- paste(GVC@names[[GVC@Matrix[i, i]]], ",",
      GVC@names[[GVC@Matrix[d, i]]],
      sep = ""
    )
    mm <- GVC@model[[count]]
    if (valid.gamBiCop(mm) != TRUE) {
      a <- paste(a, ": ", bicopname(mm$family), sep = "")
      cat(a, "\n")
    } else {
      cat(a, ": ")
      show(mm)
    }
    # cat(a, "\n")
    # if(GVC@model[d,i]!=0) { show(GVC@model[d,i]) }
    count <- count + 1
  }
  for (j in 2:(d - 1)) {
    a <- paste("Tree ", j, ":", sep = "")
    cat("\n", a, "\n")
    for (i in 1:(d - j)) {
      a <- paste(GVC@names[[GVC@Matrix[i, i]]], ",",
        GVC@names[[GVC@Matrix[d - j + 1, i]]],
        sep = ""
      )
      a <- paste(a, "|", sep = "")
      conditioningSet <- (d - j + 2):d
      for (k in 1:length(conditioningSet)) {
        if (k > 1) {
          a <- paste(a, ",", sep = "")
        }
        a <- paste(a, GVC@names[[GVC@Matrix[conditioningSet[k], i]]], sep = "")
      }
      mm <- GVC@model[[count]]
      if (valid.gamBiCop(mm) != TRUE) {
        a <- paste(a, ": ", bicopname(mm$family), sep = "")
        cat(a, "\n")
      } else {
        cat(a, ": ")
        show(mm)
      }
      count <- count + 1
    }
  }
}

summary.gamVine <- function(object, ...) {
  detail <- FALSE
  GVC <- object
  cat("GAM-Vine matrix:", "\n")
  show(GVC@Matrix)
  cat("\n", "Where", "\n")
  for (i in 1:length(GVC@names)) {
    cat(i, " <-> ", GVC@names[[i]], "\n")
  }

  d <- dim(GVC)
  count <- 1
  cat("\n", "Tree 1:", "\n")
  for (i in 1:(d - 1)) {
    mm <- GVC@model[[count]]
    a <- paste(GVC@names[[GVC@Matrix[i, i]]], ",",
      GVC@names[[GVC@Matrix[d, i]]],
      sep = ""
    )
    if (valid.gamBiCop(mm) != TRUE) {
      a <- paste(a, ": ", bicopname(mm$family), sep = "")
      if (mm$family != 0) {
        a <- paste(a, " with par=", round(mm$par, 2), sep = "")
        if (mm$family %in% c(
          2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40, 104,
          114, 124, 134, 204, 214, 224, 234
        )) {
          a <- paste(a, " and par2=", round(mm$par2, 2), sep = "")
        }
        a <- paste(a, " (tau=", round(par2tau(mm$par, mm$family), 2), ")", sep = "")
      }
      cat(a, "\n")
    } else {
      cat(a, ": ")
      summary(mm)
    }
    count <- count + 1
  }
  for (j in 2:(d - 1)) {
    a <- paste("Tree ", j, ":", sep = "")
    cat("\n", a, "\n")
    for (i in 1:(d - j)) {
      a <- paste(GVC@names[[GVC@Matrix[i, i]]], ",",
        GVC@names[[GVC@Matrix[d - j + 1, i]]],
        sep = ""
      )
      a <- paste(a, "|", sep = "")
      conditioningSet <- (d - j + 2):d
      for (k in 1:length(conditioningSet)) {
        if (k > 1) {
          a <- paste(a, ",", sep = "")
        }
        a <- paste(a, GVC@names[[GVC@Matrix[conditioningSet[k], i]]], sep = "")
      }
      mm <- GVC@model[[count]]
      if (valid.gamBiCop(mm) != TRUE) {
        a <- paste(a, ": ", bicopname(mm$family), sep = "")
        if (mm$family != 0) {
          a <- paste(a, " with par=", round(mm$par, 2), sep = "")
          if (mm$family %in% c(
            2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40,
            104, 114, 124, 134, 204, 214, 224, 234
          )) {
            a <- paste(a, " and par2=", round(mm$par2, 2), sep = "")
          }
          a <- paste(a, " (tau=", round(par2tau(mm$par, mm$family), 2), ")", sep = "")
        }
        cat(a, "\n")
      } else {
        cat(a, ": ")
        summary(mm)
      }
      count <- count + 1
    }
  }
}

plot.gamVine <- function(x, ...) {
  sel <- sapply(x@model, function(y) valid.gamBiCop(y) == "TRUE" &&
      length(summary(y@model)$s.pv) > 0)
  covariates <- x@covariates
  if (!any(is.na(covariates))) {
    l <- length(covariates)
  } else {
    l <- 0
  }
  if (any(sel)) {
    d <- (1 + sqrt(1 + 8 * length(x@model))) / 2
    model.count <- get.modelCount(d)
    M <- x@Matrix
    nn <- x@names
    par(ask = T)
    sel <- which(sel)
    for (j in sel) {
      mm <- x@model[[j]]
      k <- arrayInd(which(model.count == j), dim(M))
      con1 <- paste0(sort(c(
        nn[M[k[1], k[2]]],
        nn[M[which(M[, k[2]] != 0)[1], k[2]]]
      )),
      collapse = ","
      )
      if (k[1] != d) {
        con1 <- paste0(c(con1, paste0(sort(nn[M[(k[1] + 1):d, k[2]]]),
          collapse = ","
        )), collapse = ";")
      }
      con2 <- paste0(sort(all.vars(mm@model$pred.formula)), collapse = ",")

      args <- list(
        x = x@model[[j]],
        main = paste0(c(con1, con2), collapse = "|"),
        se = F
      )
      args <- modifyList(args, list(...))
      do.call(plot, args)
    }
  }
}

setValidity("gamVine", valid.gamVine)
setMethod("show", signature("gamVine"), show.gamVine)

#' Summary for an Object of the Class gamVine
#'
#' Takes an object of the class
#' \code{\link[gamCopula:gamVine-class]{gamVine}} and produces various
#' useful summaries from it.
#'
#' @param object An object of the class
#' \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @param ... unused in this class
#' @return A useful summary (see \code{\link{summary.gam}}
#' from \code{\link[mgcv:mgcv-package]{mgcv}} for more details).
#' @seealso \code{\link{summary.gam}}
#' from \code{\link[mgcv:mgcv-package]{mgcv}}
#' @docType methods
#' @name summary.gamVine
#' @rdname summary.gamVine-methods
#' @aliases summary,gamVine-method
#' @export
setMethod("summary", signature("gamVine"), summary.gamVine)

#' Dimension of an Object of the Class gamVine
#'
#' Retrieve the dimension of an object of the class
#' \code{\link[gamCopula:gamVine-class]{gamVine}}.
#'
#' @param x An object of the class
#' \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @return Dimension of the
#' \code{\link[gamCopula:gamVine-class]{gamVine}} object.
#' @seealso \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @docType methods
#' @name dim.gamVine
#' @rdname dim.gamVine-methods
#' @aliases dim,gamVine-method
#' @export
setMethod("dim", signature("gamVine"), dim.gamVine)

#' Plot an Object of the Class gamVine
#'
#' Plot an object of the class
#' \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' The function is based on (see \code{\link{plot.gam}}
#' from \code{\link[mgcv:mgcv-package]{mgcv}}).
#'
#' @param x An object of the class
#' \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @param y Not used with this class.
#' @param ... additional arguments to be passed to \code{\link{plot.gam}}.
#' @return This function simply generates plots.
#' @seealso \code{\link{plot.gam}} from \code{\link[mgcv:mgcv-package]{mgcv}}).
#' @docType methods
#' @name plot.gamVine
#' @rdname plot.gamVine-methods
#' @aliases plot,gamVine,ANY-method
#' @export
setMethod("plot", signature(x = "gamVine"), plot.gamVine)

#' Normalize an Object of the Class gamVine
#'
#' Change the R-vine matrix in the natural order,
#' i.e. with d:1 on the diagonal
#'
#' @param GVC An object of the class
#' \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @return The normalized \code{\link[gamCopula:gamVine-class]{gamVine}} object.
#' @seealso \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @name gamVineNormalize
#' @rdname gamVineNormalize
#' @export
gamVineNormalize <- function(GVC) {
  if (any(!isS4(GVC), !is(GVC, "gamVine"))) {
    stop("'GVC' has to be an gamVine object.")
  }

  oldOrder <- diag(GVC@Matrix)
  Matrix <- reorderRVineMatrix(GVC@Matrix)

  return(gamVine(Matrix, GVC@model,
    names = rev(GVC@names[oldOrder]), covariates = GVC@covariates
  ))
}


#' Family Matrix of an Object of the Class gamVine
#'
#' Return the matrix of copula family (see
#' \code{\link[gamCopula:gamBiCop-class]{gamBiCop}})
#' corresponding to the model list in the
#' \code{\link[gamCopula:gamVine-class]{gamVine}} object.
#'
#' @param GVC An object of the class
#'  \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @return Matrix of copula families corresponding to the model list in the
#' \code{\link[gamCopula:gamVine-class]{gamVine}} object.
#' @seealso \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @name gamVineFamily
#' @rdname gamVineFamily
#' @export
gamVineFamily <- function(GVC) {
  if (any(!isS4(GVC), !is(GVC, "gamVine"))) {
    stop("'GVC' has to be an gamVine object.")
  }

  d <- dim(GVC)
  fam <- rep(0, d^2)
  temp <- sapply(GVC@model, function(x) if (isS4(x) && is(x, "gamBiCop")) {
      x@family
    } else {
      x$family
    })
  sel <- seq(d, d^2 - d, by = d)
  t1 <- 1
  for (i in 1:(d - 1)) {
    t2 <- t1 + d - i - 1
    fam[sel[1:(d - i)] - i + 1] <- temp[t1:t2]
    t1 <- t2 + 1
  }

  return(matrix(fam, d, d))
}

#' Transform an Object of the Class R-Vine into an Object of the Class gamVine
#'
#' Transform an object of the class \code{\link{RVineMatrix}}
#' into an object of the class \code{\link[gamCopula:gamVine-class]{gamVine}}.
#'
#' @param RVM An object of the class \code{\link{RVineMatrix}}.
#' @return An object of the class
#' \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @seealso \code{\link{RVineMatrix}} and
#' \code{\link[gamCopula:gamVine-class]{gamVine}}.
#' @name RVM2GVC
#' @rdname RVM2GVC
#' @export
RVM2GVC <- function(RVM) {
  if (!is(RVM, "RVineMatrix")) {
    stop("RVM has to be an object of the class RVineMatrix.")
  }

  if (!all(RVM$family[lower.tri(RVM$family)] %in% c(0:5, 13, 14, 23, 24, 33, 34))) {
    msg <- paste(
      "Only the Gaussian, t, Clayton, Gumbel and Frank copulas and",
      "rotations are allowed in gamVine objects."
    )
    stop(msg)
  }

  d <- dim(RVM$Matrix)[1]

  sel <- seq(d, d^2 - d, by = d)
  family <- RVM$family[sel]
  par <- RVM$par[sel]
  par2 <- RVM$par2[sel]
  for (j in 2:(d - 1)) {
    family <- c(family, RVM$family[sel[1:(d - j)] - j + 1])
    par <- c(par, RVM$par[sel[1:(d - j)] - j + 1])
    par2 <- c(par2, RVM$par2[sel[1:(d - j)] - j + 1])
  }

  family <- sapply(family, famTrans)
  out <- cbind(family, par, par2)
  colnames(out) <- c()
  out <- apply(out, 1, function(x) list(family = x[1], par = x[2], par2 = x[3]))
  if (!is.null(RVM$names)) {
    out <- gamVine(RVM$Matrix, out, RVM$names)
  } else {
    nnames <- paste("X", 1:d, sep = "")
    out <- gamVine(RVM$Matrix, out, nnames)
  }
  return(out)
}
