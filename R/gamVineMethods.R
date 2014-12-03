
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
  Matrix <- VineCopula:::reorderRVineMatrix(GVC@Matrix)
  
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

#' Print a \code{\link{gamVine-class}} object
#'
#' @param x \code{\link{gamVine-class}} object.
#' @param ... un-used for this class.
## @param detail should additional details be printed (\code{detail=TRUE}) or
## not?
#' @seealso \code{\link{gamVine-class}} and \code{\link{gamVine}}.
#' @docType methods
#' @rdname print-methods
#' @export
print.gamVine <- function(x, ...) {
  ## TODO , detail=TRUE
  detail <- FALSE
  GVC <- x
  message("gam-vine matrix:")
  print(GVC@Matrix, ...)
  message("")
  message("Where")
  for (i in 1:length(GVC@names)) {
    message(i, " <-> ", GVC@names[[i]])
  }
  
  d <- dim(GVC)
  count <- 1
  if (detail == TRUE || detail == T) {
    message("")
    message("Tree 1:")
    for (i in 1:(d - 1)) {
      a <- paste(GVC@names[[GVC@Matrix[i, i]]], ",", GVC@names[[GVC@Matrix[d, 
        i]]], sep = "")
      if (is.numeric(GVC@model[[count]])) {
        a <- paste(a, ": ", BiCopName(GVC@model[[count]], short = FALSE), 
          sep = "")
      } else {
        a <- paste(a, ": ", BiCopName(GVC@model[[count]]@family, short = FALSE), 
          sep = "")
      }
      message(a)
      # if(GVC@model[d,i]!=0) { print(GVC@model[d,i]) }
      count <- count + 1
    }
    for (j in 2:(d - 1)) {
      message("")
      a <- paste("Tree ", j, ":", sep = "")
      message(a)
      for (i in 1:(d - j)) {
        a <- paste(GVC@names[[GVC@Matrix[i, i]]], ",", GVC@names[[GVC@Matrix[d - 
          j + 1, i]]], sep = "")
        a <- paste(a, "|", sep = "")
        conditioningSet <- (d - j + 2):d
        for (k in 1:length(conditioningSet)) {
          if (k > 1) {
          a <- paste(a, ",", sep = "")
          }
          a <- paste(a, GVC@names[[GVC@Matrix[conditioningSet[k], i]]], 
                     sep = "")
        }
        if (is.numeric(GVC@model[[count]])) {
          a <- paste(a, ": ", BiCopName(GVC@model[[count]], 
                                        short = FALSE), sep = "")
        } else {
          EDF <- EDF(GVC@model[[count]])
          a <- paste(a, ": ", BiCopName(GVC@model[[count]]@family, 
                                        short = FALSE), sep = "")
          a <- paste(a, paste("EDF:", paste(round(EDF[-1], 3), 
                                            collapse = ", ")),sep = "\n")
        }
        message(a)
        # if(GVC@model[d-j+1,i]!=0) { if(GVC@model[d-j+1,i]!=0) {
        # print(GVC@model[d-j+1,i]) } }
        count <- count + 1
      }
    }
    
  }
}
#' @docType methods
#' @rdname print-methods
setMethod("print", signature("gamVine"), print.gamVine)

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