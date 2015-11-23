#' Selection and Maximum penalized likelihood estimation of a Generalized 
#' Additive model (gam) for the copula parameter or Kendall's tau.
#'
#' This function selects an appropriate bivariate copula family for given 
#' bivariate copula data using one of a range of methods. The corresponding 
#' parameter estimates are obtained by maximum penalized  likelihood estimation,
#' where each Newton-Raphson iteration is reformulated as a generalized ridge 
#' regression solved using the \code{\link[mgcv:mgcv-package]{mgcv}} package. 
#' 
#' @param data A list or data frame containing the model responses, (u1,u2) in 
#' [0,1]x[0,1], and covariates required by the formula.
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
#' @param selcrit Character indicating the criterion for bivariate copula 
#' selection. Possible choices: \code{selcrit = 'AIC'} (default) or 
#' \code{'BIC'}, as in \code{\link{BiCopSelect}} from the 
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package. 
#' @param level Numerical; significance level of the test for removing individual
#' predictors (default: \code{level = 0.05}).  
#' @param edf Numerical; if the estimated EDF for individual predictors is 
#' smaller than \code{edf} but the predictor is still significant, then
#' it is set as linear (default: \code{edf = 1.5}). 
#' @param tau \code{FALSE} for a calibration fonction specified for 
#' the Copula parameter or \code{TRUE} (default) for a calibration function 
#' specified for Kendall's tau.  
#' @param method \code{'FS'} for Fisher-scoring (default) and 
#' \code{'NR'} for Newton-Raphson.  
#' @param tol.rel Relative tolerance for \code{'FS'}/\code{'NR'} algorithm.  
#' @param n.iters Maximal number of iterations for 
#' \code{'FS'}/\code{'NR'} algorithm.
#' @param parallel \code{TRUE} for a parallel estimation across copula families.
#' As the code is based on mclapply, this parameter has no effect on windows.
#' @param verbose \code{TRUE} prints informations during the estimation.
#' @param ... Additional parameters to be passed to \code{\link{gam}} 
#' @return \code{gamBiCopEst} returns a list consisting of 
#' \item{res}{S4 \code{\link{gamBiCop-class}} object.} 
#' \item{method}{\code{'FS'} for Fisher-scoring and 
#' \code{'NR'} for Newton-Raphson.} 
#' \item{tol.rel}{relative tolerance for \code{'FS'}/\code{'NR'} algorithm.} 
#' \item{n.iters}{maximal number of iterations for 
#' \code{'FS'}/\code{'NR'} algorithm.} 
#' \item{trace}{the estimation procedure's trace.} 
#' \item{conv}{\code{0} if the algorithm converged and \code{1} otherwise.}
#' @seealso \code{\link{gamBiCop}} and \code{\link{gamBiCopEst}}. 
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
#' calib.surf <- list(calib.lin <- function(t, Ti = 0, Tf = 1, b = 2) {
#'     return(-2+4*t)},
#'   calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
#'     Tm <- (Tf - Ti)/2
#'     a <- -(b/3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
#'   return(a + b * (t - Tm)^2)},
#'   calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1) {
#'     a <- b * (1 - 2 * Tf * pi/(f * Tf * pi +
#'                                  cos(2 * f * pi * (Tf - Ti))
#'                                - cos(2 * f * pi * Ti)))
#'     return((a + b)/2 + (b - a) * sin(2 * f * pi * (t - Ti))/2)},
#'   calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf/8) {
#'     Tm <- (Tf - Ti)/2
#'     a <- (b * s * sqrt(2 * pi)/Tf) * (pnorm(0, Tm, s) - pnorm(Tf, Tm, s))
#'     return(a + b * exp(-(t - Tm)^2/(2 * s^2)))})
#' 
#' ## 6-dimensional matrix X of covariates
#' covariates.distr <- mvdc(normalCopula(rho, dim = 6),
#'                                  c("unif"), list(list(min = 0, max = 1)),
#'                                  marginsIdentical = TRUE)
#' X <- rMvdc(n, covariates.distr)
#' 
#' ## U in [0,1]x[0,1] depending on the four first columns of X
#' U <- condBiCopSim(fam, function(x1,x2,x3,x4) {eta0+sum(mapply(function(f,x)
#'   f(x), calib.surf, c(x1,x2,x3,x4)))}, X[,1:4], par2 = 4, return.par = TRUE)
#' 
#' ## Merge U and X
#' data <- data.frame(U$data,X)
#' names(data) <- c(paste("u",1:2,sep=""),paste("x",1:6,sep=""))
#' 
#' ## Selection using AIC (take about 5mn on single core) 
#' ## Use parallel = TRUE to speed-up....
#' system.time(best <- gamBiCopSel(data))
#' print(best$res)
#' EDF(best$res)
#' plot(best$res)
#' @export
gamBiCopSel <- function(data, familyset = NA, rotations = TRUE, 
                        selcrit = "AIC", level = 5e-2, edf = 1.5, tau = TRUE, 
                        method = "FS", tol.rel = 1e-3, n.iters = 10,
                        parallel = FALSE, verbose = FALSE, ...) {
  
  tmp <- valid.gamBiCopSel(data, rotations, familyset, selcrit, level, edf, tau, 
                           method, tol.rel, n.iters, parallel, verbose)
  if (tmp != TRUE)
    stop(tmp)
  
  if (is.list(data)){
    if(!is.null(data$xt)){
      xt <- data$xt
      data <- data[-which(names(data) == "xt")]
    }
    data <- as.data.frame(data)
  }
  
  n <- dim(data)[1]
  m <- dim(data)[2]
  u1 <- data$u1
  u2 <- data$u2
  u <- cbind(u1,u2)
  
  if (length(familyset) == 1 && is.na(familyset)) {
    familyset <- get.familyset()
  }
  if (rotations) {
    familyset <- withRotations(familyset)
  }
  
  parallel <- as.logical(parallel)
  x <- NULL
  if (parallel == FALSE) {
    res <- foreach(x=familyset) %do% 
      tryCatch(gamBiCopVarSel(data, x, tau, method, tol.rel, n.iters, 
                              level, edf, verbose,...),
               error = function(e) e)
  } else {
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl, cores = detectCores() - 1)
    res <- foreach(x=familyset) %dopar% 
      tryCatch(gamBiCopVarSel(data, x, tau, method, tol.rel, n.iters, 
                              level, edf, verbose,...),
               error = function(e) e)
    stopCluster(cl)
  }
  res <- res[sapply(res,length) == 6] 
  if (length(res) == 0 ) { #|| sum(sapply(res, function(x) x$conv) == 0) == 0) {
    return(paste("No convergence of the estimation for any copula family.",
                 "Try modifying the parameters (FS/NR, tol.rel, n.iters)..."))
  }
  #   } else {
  #     res <- res[sapply(res, function(x) x$conv) == 0]
  #   }
  
  if (selcrit == "AIC") {
    return(res[[which.min(sapply(res, function(x) AIC(x$res)))]])
  } else {
    return(res[[which.min(sapply(res, function(x) BIC(x$res)))]])
  }
} 

gamBiCopVarSel <- function(data, family,
                           tau = TRUE, method = "FS",
                           tol.rel = 1e-3, n.iters = 10, 
                           level = 5e-2, edf = 1.5, 
                           verbose = FALSE, ...) {
  
  if (verbose == TRUE) {
    cat(paste("Model selection for family", family, "\n"))
  }
  
  myGamBiCopEst <- function(formula) 
    suppressMessages(gamBiCopEst(data, formula, family, tau, 
                                 method, tol.rel, n.iters))
  
  n <- dim(data)[1]
  d <- dim(data)[2]-2
  ## Create a list with formulas of smooth terms corresponding to the covariates  
  nn <- names(data)[-which( (names(data) == "u1") | (names(data) == "u2"))]
  
  basis <- rep(5,d)
  formula.expr <- mapply(get.smooth,nn,basis)
  formula.lin <- NULL
  
  ## Update the list by removing unsignificant predictors 
  if (verbose == TRUE) {
    cat("Remove unsignificant covariates.......\n")
  }
  
  sel.smooth <- smooth2lin <- FALSE
  sel.lin <- NULL
  while((!all(c(sel.lin,sel.smooth)) && length(basis) > 0) ||
          any(smooth2lin)) {
    
    smooth2lin <- FALSE
    sel.lin <- NULL
    formula.tmp <- get.formula(c(formula.lin,formula.expr))
    
    if (verbose == TRUE) {
      cat("Model formula:\n")
      print(formula.tmp)
    }
    
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
      basis <- rep(5,length(nn))
      formula.expr <- mapply(get.smooth,nn,basis)
      
      ## Check whether some smooth need to be set to linear
      ## (and adjust the model if required)
      smooth2lin <- sel.smooth & summary(tmp$res@model)$edf < edf
      if (any(smooth2lin)) {
        sel.lin <- which(sel.smooth) 
        sel.lin <- summary(tmp$res@model)$edf[sel.lin] < edf
        sel.smooth <- !sel.lin
        formula.lin <- c(formula.lin,nn[sel.lin])
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
  }  else {
    formula.tmp <- get.formula(c(formula.lin,formula.expr))
    tmp <- myGamBiCopEst(formula.tmp)
  }
  
  ## Increasing the basis size appropriately
  
  tmp2 <- get.pval(tmp$res@model)
  
  sel <- tmp2[,1]-tmp2[,2] < 1 & tmp2[,4] < level
  #sel <- summary(tmp$res@model)$edf > (basis-1)/2
  if (verbose == TRUE && any(sel)) {
    cat(paste("Select the basis sizes .......\n"))
  }
  #kk <- round(exp(seq(log(8), log(n/30),length.out=20)))
  #kcount <- 1
  #while (any(sel) && kcount < 20) {
  while (any(sel) && all(basis < n/30)) {
    
    ## Increase basis sizes
    #basis[sel] <- kk[kcount]
    #kcount <- kcount + 1
    #basis[sel] <- round(basis[sel]*1.5)
    basis[sel] <- 2*basis[sel]
    
    #     ## Extract and fit model to residuals for each smooth components
    #     data$y <- residuals(tmp$res@model)
    #     data$w <- tmp$res@model$weights
    #     residuals.expr <- mapply(function(x,y) 
    #       get.smooth(x,y,bs="cs"), nn[sel], basis[sel])
    #     residuals.formula <- sapply(residuals.expr,function(expr) 
    #       get.formula(expr,TRUE))
    #     residuals.fit <- lapply(residuals.formula, function(f) 
    #       gam(f, data = data, gamma = 1.4, weights = w)) 
    # 
    #     ## Suspect residuals
    #     residuals.edf <- sapply(residuals.fit,function(x) sum(x$edf)-1) 
    #     sel[sel] <- residuals.edf > (basis[sel]-1)/2
    
    ## Update the smooth terms, formula and fit the new model
    if (any(sel)) {
      formula.expr[sel] <- mapply(get.smooth,nn[sel],basis[sel])
      formula.tmp <- get.formula(c(formula.lin,formula.expr))
      if (verbose == TRUE) {
        cat("Updated model formula:\n")
        print(formula.tmp)
      }      
      tmp <- myGamBiCopEst(formula.tmp)
      #sel[sel] <- summary(tmp$res@model)$edf[sel] > (basis[sel]-1)/2
      tmp2 <- get.pval(tmp$res@model)
      sel <- sel & tmp2[,1]-tmp2[,2] < 1 & tmp2[,4] < level
    }
    
    ## Check that the basis size is smaller than 1/2 the size of 
    ## the corresponding covariate's support
    if (any(sel)) {
      if (sum(sel) == 1) {
        sel[sel] <- 2*basis[sel] < length(unique(data[,nn[sel]]))
      } else {
        sel[sel] <- 2*basis[sel] < apply(data[,nn[sel]], 2, function(x) 
          length(unique(x)))
      } 
    }
  }
  
  return(tmp)
}

get.smooth <- function(x,k,bs="cr") {
  paste("s(",x,", k=", k,", bs='", bs, "')",sep = "")
} 
get.formula <- function(expr,y=FALSE) {
  if (!y) {
    as.formula(paste("~",paste(expr,collapse = " + ")))
  } else {
    as.formula(paste("y ~",paste(expr,collapse = " + ")))
  }
} 

get.linear <- function(x,nn){
  names(unlist(sapply(nn,function(z)grep(z,x))))
}

valid.gamBiCopSel <- function(data, rotations, familyset, selcrit, level, edf, 
                              tau, method, tol.rel, n.iters, parallel, 
                              verbose) {
  
  tmp <- valid.gamBiCopEst(data, n.iters, tau, tol.rel, method, verbose, 1)
  if (tmp != TRUE) {
    return(tmp)
  }
  
  if (!valid.familyset(familyset)) {
    return(return(msg.familyset(var2char(familyset))))
  }
  
  if (!valid.logical(rotations)) {
    return(msg.logical(var2char(rotations)))
  }
  
  if(is.null(selcrit) || length(selcrit) != 1 || 
       (selcrit != "AIC" && selcrit != "BIC")) {
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
  
  return(TRUE)
}

## Internal code borrowed from mgcv
get.pval <- function (b, subsample = 5000, n.rep = 400) 
{
  m <- length(b$smooth)
  if (m == 0) 
    return(NULL)
  rsd <- residuals(b)
  ve <- rep(0, n.rep)
  p.val <- v.obs <- kc <- edf <- rep(0, m)
  snames <- rep("", m)
  n <- nrow(b$model)
  if (n > subsample) {
    ind <- sample(1:n, subsample)
    modf <- b$model[ind, ]
    rsd <- rsd[ind]
  }
  else modf <- b$model
  nr <- length(rsd)
  for (k in 1:m) {
    ok <- TRUE
    b$smooth[[k]]$by <- "NA"
    dat <- ExtractData(b$smooth[[k]], modf, NULL)$data
    if (!is.null(attr(dat, "index")) || 
          !is.null(attr(dat[[1]], "matrix")) || is.matrix(dat[[1]])) 
      ok <- FALSE
    if (ok) 
      dat <- as.data.frame(dat)
    snames[k] <- b$smooth[[k]]$label
    ind <- b$smooth[[k]]$first.para:b$smooth[[k]]$last.para
    kc[k] <- length(ind)
    edf[k] <- sum(b$edf[ind])
    nc <- b$smooth[[k]]$dim
    if (ok && ncol(dat) > nc) 
      dat <- dat[, 1:nc, drop = FALSE]
    for (j in 1:nc) if (is.factor(dat[[j]])) 
      ok <- FALSE
    if (!ok) {
      p.val[k] <- v.obs[k] <- NA
    }
    else {
      if (nc == 1) {
        e <- diff(rsd[order(dat[, 1])])
        v.obs[k] <- mean(e^2)/2
        for (i in 1:n.rep) {
          e <- diff(rsd[sample(1:nr, nr)])
          ve[i] <- mean(e^2)/2
        }
        p.val[k] <- mean(ve < v.obs[k])
        v.obs[k] <- v.obs[k]/mean(rsd^2)
      }
      else {
        if (!is.null(b$smooth[[k]]$margin)) {
          beta <- coef(b)[ind]
          f0 <- PredictMat(b$smooth[[k]], dat) %*% beta
          gr.f <- rep(0, ncol(dat))
          for (i in 1:nc) {
            datp <- dat
            dx <- diff(range(dat[, i]))/1000
            datp[, i] <- datp[, i] + dx
            fp <- PredictMat(b$smooth[[k]], datp) %*% 
              beta
            gr.f[i] <- mean(abs(fp - f0))/dx
          }
          for (i in 1:nc) {
            dat[, i] <- dat[, i] - min(dat[, i])
            dat[, i] <- gr.f[i] * dat[, i]/max(dat[, 
                                                   i])
          }
        }
        nn <- 3
        ni <- nearest(nn, as.matrix(dat))$ni
        e <- rsd - rsd[ni[, 1]]
        for (j in 2:nn) e <- c(e, rsd - rsd[ni[, j]])
        v.obs[k] <- mean(e^2)/2
        for (i in 1:n.rep) {
          rsdr <- rsd[sample(1:nr, nr)]
          e <- rsdr - rsdr[ni[, 1]]
          for (j in 2:nn) e <- c(e, rsdr - rsdr[ni[, 
                                                   j]])
          ve[i] <- mean(e^2)/2
        }
        p.val[k] <- mean(ve < v.obs[k])
        v.obs[k] <- v.obs[k]/mean(rsd^2)
      }
    }
  }
  k.table <- cbind(kc, edf, v.obs, p.val)
  dimnames(k.table) <- list(snames, c("k'", "edf", "k-index", 
                                      "p-value"))
  k.table
}
## Internal code borrowed from mgcv
ExtractData <- function (object, data, knots) {
  knt <- dat <- list()
  for (i in 1:length(object$term)) {
    dat[[object$term[i]]] <- get.var(object$term[i], data)
    knt[[object$term[i]]] <- get.var(object$term[i], knots)
  }
  names(dat) <- object$term
  m <- length(object$term)
  if (!is.null(attr(dat[[1]], "matrix"))) {
    n <- length(dat[[1]])
    X <- matrix(unlist(dat), n, m)
    if (is.numeric(X)) {
      X <- uniquecombs(X)
      if (nrow(X) < n * 0.9) {
        for (i in 1:m) dat[[i]] <- X[, i]
        attr(dat, "index") <- attr(X, "index")
      }
    }
  }
  if (object$by != "NA") {
    by <- get.var(object$by, data)
    if (!is.null(by)) {
      dat[[m + 1]] <- by
      names(dat)[m + 1] <- object$by
    }
  }
  return(list(data = dat, knots = knt))
}