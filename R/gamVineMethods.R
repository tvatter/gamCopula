
#' Normalize a \code{\link{gamVine-class}} object
#' 
#' Change the R-vine matrix in the natural order, 
#' i.e. with d:1 on the diagonal
#'
#' @param GVC fitted \code{\link{gamVine-class}} object.
#' @return Normalized \code{\link{gamVine-class}} object.
#' @seealso \code{\link{gamVine-class}} and \code{\link{gamVine}}.
#' @export
gamVineNormalize <- function(GVC) {
  if (any(!isS4(GVC), !is(GVC, "gamVine"))) 
    stop("'GVC' has to be an gamVine object.")
  
  oldOrder <- diag(GVC@Matrix)
  Matrix <- reorderRVineMatrix(GVC@Matrix)
  
  return(gamVine(Matrix, GVC@model, names = rev(GVC@names[oldOrder])))
}

#' Dimension of a \code{\link{gamVine-class}} object
#' 
#' Retrieve the dimension of a \code{\link{gamVine-class}} object.
#'
#' @param x fitted \code{\link{gamVine-class}} object.
#' @return Dimension of the \code{\link{gamVine-class}} object.
#' @seealso \code{\link{gamVine-class}} and \code{\link{gamVine}}.
#' @docType methods
#' @rdname dim-methods
#' @export
dim.gamVine <- function(x) {
  GVC <- x
  return(dim(GVC@Matrix)[1])
}
#' @docType methods
#' @rdname dim-methods
setMethod("dim", signature("gamVine"), dim.gamVine)

show.gamVine <- function(object) {
  detail <- FALSE
  GVC <- object
  cat("GAM-Vine matrix:","\n")
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
               GVC@names[[GVC@Matrix[d,i]]], sep = "")
    a <- paste(a, ": ", BiCopName(GVC@model[[count]]$family,
                                  short = FALSE), sep = "")
    cat(a, "\n")
    # if(GVC@model[d,i]!=0) { show(GVC@model[d,i]) }
    count <- count + 1
  }
  for (j in 2:(d - 1)) {
    a <- paste("Tree ", j, ":", sep = "")
    cat("\n", a, "\n")
    for (i in 1:(d - j)) {
      a <- paste(GVC@names[[GVC@Matrix[i, i]]], ",", 
                 GVC@names[[GVC@Matrix[d - j + 1, i]]], sep = "")
      a <- paste(a, "|", sep = "")
      conditioningSet <- (d - j + 2):d
      for (k in 1:length(conditioningSet)) {
        if (k > 1) {
          a <- paste(a, ",", sep = "")
        }
        a <- paste(a, GVC@names[[GVC@Matrix[conditioningSet[k], i]]],sep = "")
      }
      mm <- GVC@model[[count]]
      if ( valid.gamBiCop(mm) != TRUE) {
        a <- paste(a, ": ", BiCopName(mm$family, short = FALSE), sep = "")
        cat(a, "\n")
      } else {
        cat(a, ": ")
        show(mm)
      }
      count <- count + 1
    }
  }
}
setMethod("show", signature("gamVine"), show.gamVine)

summary.gamVine <- function(object) {
  detail <- FALSE
  GVC <- object
  cat("GAM-Vine matrix:","\n")
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
               GVC@names[[GVC@Matrix[d,i]]], sep = "")
    a <- paste(a, ": ", BiCopName(mm$family, short = FALSE), sep = "")
    if (mm$family!=0) {
      a <- paste(a, " with par=", round(mm$par,2), sep="")
      if (mm$family %in% c(2,7,8,9,10,17,18,19,20,27,28,29,30,37,38,39,40,104,
                          114,124,134,204,214,224,234)) {
        a <- paste(a," and par2=", round(mm$par2,2), sep="")
      }
      a <- paste(a, " (tau=", 
                 round(BiCopPar2Tau(mm$family,mm$par,mm$par2),2), ")", sep="")
    }
    cat(a, "\n")
    count <- count + 1
  }
  for (j in 2:(d - 1)) {
    a <- paste("Tree ", j, ":", sep = "")
    cat("\n", a, "\n")
    for (i in 1:(d - j)) {
      a <- paste(GVC@names[[GVC@Matrix[i, i]]], ",", 
                 GVC@names[[GVC@Matrix[d - j + 1, i]]], sep = "")
      a <- paste(a, "|", sep = "")
      conditioningSet <- (d - j + 2):d
      for (k in 1:length(conditioningSet)) {
        if (k > 1) {
          a <- paste(a, ",", sep = "")
        }
        a <- paste(a, GVC@names[[GVC@Matrix[conditioningSet[k], i]]],sep = "")
      }
      mm <- GVC@model[[count]]
      if ( valid.gamBiCop(mm) != TRUE) {
        a <- paste(a, ": ", BiCopName(mm$family, short = FALSE), sep = "")
        if (mm$family!=0) {
          a <- paste(a, " with par=", round(mm$par,2), sep="")
          if (mm$family %in% c(2,7,8,9,10,17,18,19,20,27,28,29,30,37,38,39,40,
                               104,114,124,134,204,214,224,234)) {
            a <- paste(a," and par2=", round(mm$par2,2), sep="")
          }
          a <- paste(a, " (tau=", 
                     round(BiCopPar2Tau(mm$family,mm$par,mm$par2),2), 
                     ")", sep="")
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
setMethod("summary", signature("gamVine"), summary.gamVine)

#' Family matrix of \code{\link{gamVine-class}} object
#' 
#' Return the matrix of copula family (see \code{\link{gamBiCop-class}}) 
#' corresponding to the model list in the \code{\link{gamVine-class}} object.
#'
#' @param GVC \code{\link{gamVine-class}} object.
#' @return Matrix of copula family (see \code{\link{gamBiCop-class}}) 
#' corresponding to the model list in the \code{\link{gamVine-class}} object.
#' @seealso \code{\link{gamVine-class}} and \code{\link{gamVine}}.
#' @export
gamVineFamily <- function(GVC) {
  
  if (any(!isS4(GVC), !is(GVC, "gamVine"))) 
    stop("'GVC' has to be an gamVine object.")
  
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

#' \code{\link{RVineMatrix}} to \code{\link{gamVine-class}}
#' 
#' Transform a \code{\link{RVineMatrix}} object 
#' into a \code{\link{gamVine-class}} object.
#'
#' @param RVM \code{\link{RVineMatrix}} object.
#' @return A \code{\link{gamVine-class}} object.
#' @seealso \code{\link{RVineMatrix}}, \code{\link{gamVine-class}} and 
#' \code{\link{gamVine}}.
#' @export
RVM2GVC <- function(RVM) {
  
  if (class(RVM) != "RVineMatrix") {
    stop("RVM has to be a RVineMatrix class object.")
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

#' Plot a fitted \code{\link{gamVine-class}} object
#' 
#' Plot from a model fit. 
#' The function is based on (see \code{\link{plot.gam}}
#' from \code{\link[mgcv:mgcv-package]{mgcv}}).
#' 
#' @param x fitted \code{\link{gamVine-class}} object.
#' @param ... additional arguments to be passed to \code{\link{plot.gam}}.
#' @return This function simply generates plots.
#' @seealso \code{\link{plot.gam}} from \code{\link[mgcv:mgcv-package]{mgcv}}).
#' @docType methods
#' @rdname plot-methods
#' @export
plot.gamVine <- function(x, ...) {
  sel <- sapply(x@model,valid.gamBiCop) == "TRUE"
  if (any(sel)) {
    d <- (1+sqrt(1+8*length(x@model)))/2
    model.count <- get.modelCount(d)
    M <- x@Matrix
    nn <- x@names
    par(ask = T)
    sel <- which(sel)
    for(j in sel) {
      k <- which(model.count == j)
      cc <- floor(k/d)+1
      rr <- k %% d
      if (rr == 0) {
        rr <- d
      }
      con1 <- paste(nn[M[rr,cc]],nn[M[which(M[,cc] != 0)[1],cc]], 
                    sep = ",", collapse = "")
      con2 <- paste(nn[M[(rr+1):d,cc]], sep = ",", collapse = "")
      plot(x@model[[j]], main = paste(con1,con2,sep = "|", collapse = ""), ...)
    }
  }
}
setMethod("plot", signature(x="gamVine"), plot.gamVine)