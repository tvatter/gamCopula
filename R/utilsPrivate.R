## Update parameters with respect to fited gam model
"pars.update" <- function(mm, family, data, tau = TRUE) {
  
  out <- list()
  fitted <- predict(mm, data)
  
  out$partrans <- as.numeric(fitted)
  
  if (tau == TRUE) {
    # From transformed parameters to Kendall's tau
    if (is.element(family, c(1, 2, 5))) {
      out$tau <- tanh(out$partrans/2)
      out$tau <- pmax(-1 + 1e-08, out$tau)
      out$tau <- pmin(1 - 1e-08, out$tau)
      out$partrans <- 2 * atanh(out$tau)
    } else {
      out$tau <- (1 + tanh(out$partrans/2))/2
      out$tau <- pmax(1e-08, out$tau)
      out$tau <- pmin(1 - 1e-08, out$tau)
      out$partrans <- 2 * atanh(2 * out$tau - 1)
    }
    
    # From Kendall's tau to copula parameter
    if (is.element(family, c(1, 2, 5))) {
      if (is.element(family, c(1, 2))) {
        out$par <- sapply(out$tau, function(x) BiCopTau2Par(1, x))
      } else {
        out$par <- sapply(out$tau, function(x) BiCopTau2Par(5, x))
      }
      sel1 <- out$par == 1
      sel2 <- out$par == -1
      out$par[sel1] <- 1 - 1e-08
      out$par[sel2] <- -1 + 1e-08
    } else {
      out$par <- sapply(out$tau, function(x) BiCopTau2Par(family, x))
      if (family == 3) {
        out$par <- pmax(1e-08, out$par)
      }
      if (family == 4) {
        out$par <- pmax(1 + 1e-08, out$par)
      }
    }
  } else {
    # From transformed parameters to copula parameter
    if (is.element(family, c(1, 2))) {
      out$par <- tanh(out$partrans/2)
      sel1 <- out$par == 1
      sel2 <- out$par == -1
      out$par[sel1] <- 1 - 1e-08
      out$par[sel2] <- -1 + 1e-08
      out$partrans <- 2 * atanh(out$par)
    } else if (family == 3) {
      out$par <- exp(out$partrans)
      sel1 <- out$par == Inf
      sel2 <- out$par == 0
      out$par[sel1] <- 1e+05
      out$par[sel2] <- 1e-08
      out$partrans <- log(out$par)
    } else if (family == 4) {
      out$par <- 1 + exp(out$partrans)
      sel1 <- out$par == Inf
      sel2 <- out$par == 1
      out$par[sel1] <- 1e+05
      out$par[sel2] <- 1 + 1e-08
      out$partrans <- log(out$par - 1)
    } else {
      out$par <- out$partrans
      sel1 <- out$par == Inf
      sel2 <- out$par == 0
      sel3 <- out$par == -Inf
      out$par[sel1] <- 1e+05
      out$par[sel2] <- 1e-08
      out$par[sel3] <- -1e+05
      out$partrans <- out$par
    }
  }
  
  return(out)
}

## Update trace
"trace.update" <- function(old.par, new.par) {
  
  traces <- max(abs((old.par - new.par)/old.par))
  eps <- traces
  
  out <- list()
  out$traces <- traces
  out$eps <- eps
  
  return(out)
}

## Newton-Raphson/Fisher-scoring step
"wz.update" <- function(dd, new.pars, family, method, tau = TRUE) {
  
  u <- dd$d1 * dd$dpar
  if (tau == TRUE) {
    u <- u * dd$dtau
  }
  if (method == "NR") {
    if (tau == TRUE) {
      w <- dd$dtau^2 * (dd$dpar^2 * (dd$d2/dd$p - dd$d1^2) + dd$dpar2 * dd$d1) + 
        dd$dtau2 * dd$d1 * dd$dpar * dd$dtau2
    } else {
      w <- dd$dpar^2 * (dd$d2/dd$p - dd$d1^2) + dd$dpar2 * dd$d1
    }
    w <- -w
    w[w<0] <- median(w[w>0])
  } else {
    if (family == 2) {
      if (tau == TRUE) {
        w <- (dd$dpar * dd$dtau)^2 * FisherBiCop(family, new.pars$par, new.pars$par2)
      } else {
        w <- dd$dpar^2 * FisherBiCop(family, new.pars$par, new.pars$par2)
      }
    } else {
      if (tau == TRUE) {
        w <- (dd$dpar * dd$dtau)^2 * FisherBiCop(family, new.pars$par)
      } else {
        w <- dd$dpar^2 * FisherBiCop(family, new.pars$par)
      }
    }
  }
  
  z <- new.pars$partrans + u/w
  
  out <- list()
  out$w <- w
  out$z <- z
  
  return(out)
}

## Compute first derivatives of the copula, copula parameters transformations and
## dependence measure transformations
"derivatives.par" <- function(data, new.pars, family, method, tau = TRUE) {
  
  # Derivatives of the copula with respect to its own parameter
  if (family == 2) {
    dp <- apply(data, 1, function(x) c(BiCopPDF(x[1], x[2], family = 2, x[3], 
      x[4]), BiCopDeriv(x[1], x[2], family = 2, x[3], x[4], deriv = "par", 
      log = TRUE)))
    if (method == "NR") {
      temp <- apply(data, 1, function(x) BiCopDeriv2(x[1], x[2], family = 2, 
        x[3], x[4], deriv = "par"))
      dp <- rbind(dp, temp)
    }
  } else {
    dp <- apply(data, 1, function(x) c(BiCopPDF(x[1], x[2], family = family, 
      x[3]), BiCopDeriv(x[1], x[2], family = family, x[3], deriv = "par", log = TRUE)))
    if (method == "NR") {
      temp <- apply(data, 1, function(x) BiCopDeriv2(x[1], x[2], family = family, 
        x[3], deriv = "par"))
      dp <- rbind(dp, temp)
    }
  }
  
  dp[1, which(dp[1, ] == 0)] <- 1e-16
  
  if (family == 4) {
    sel <- is.nan(dp[2, ])
    if (sum(sel) > 1) {
      dp[2, sel] <- apply(data[sel, ], 1, function(x) derivGumbel(x[3], x[1], 
        x[2]))/dp[1, sel]
    } else if (sum(sel) == 1) {
      dp[2, sel] <- derivGumbel(data[sel, 3], data[sel, 1], data[sel, 2])/dp[1, 
        sel]
    }
  }
  
  out <- list()
  out$p <- dp[1, ]
  out$d1 <- dp[2, ]
  
  if (method == "NR") {
    out$d2 <- dp[3, ]
  }
  
  if (tau == TRUE) {
    # Derivatives of the dependence measure with respect to the model parameters
    dtau <- 1/(1 + cosh(new.pars$partrans))
    if (method == "NR") {
      dtau2 <- -4 * sinh(new.pars$partrans/2)^4 * (1/sinh(new.pars$partrans))^3
    }
    
    # Derivatives of the copula parameter with respect to the dependence measure
    if (is.element(family, c(1, 2))) {
      dpar <- cos(new.pars$tau * pi/2) * pi/2
      if (method == "NR") {
        dpar2 <- -sin(new.pars$tau * pi/2) * (pi/2)^2
      }
    } else if (family == 3) {
      dpar <- 2/(new.pars$tau - 1)^2
      if (method == "NR") {
        dpar2 <- -4/(new.pars$tau - 1)^3
      }
    } else if (family == 4) {
      dpar <- 1/(new.pars$tau - 1)^2
      if (method == "NR") {
        dpar2 <- -2/(new.pars$tau - 1)^3
      }
    }
    
    out$dpar <- dpar
    out$dtau <- dtau
    if (method == "NR") {
      out$dtau2 <- dtau2
      out$dpar2 <- dpar2
    }
  } else {
    # Derivatives of the copula parameter with respect to the model parameters
    if (is.element(family, c(1, 2))) {
      dpar <- 1/(1 + cosh(new.pars$partrans))
      if (method == "NR") {
        dpar2 <- -4 * sinh(new.pars$partrans/2)^4 * (1/sinh(new.pars$partrans))^3
      }
    } else if (is.element(family, c(3, 4))) {
      dpar <- exp(new.pars$partrans)
      if (method == "NR") {
        dpar2 <- exp(new.pars$partrans)
      }
    } else {
      dpar <- rep(1, length(new.pars$partrans))
      if (method == "NR") {
        dpar2 <- rep(0, length(new.pars$partrans))
      }
    }
    
    out$dpar <- dpar
    if (method == "NR") {
      out$dpar2 <- dpar2
    }
  }
  
  return(out)
}

## Derivative of the Gumbel copula density with respect to the parameter (Higher
## precision than in the VineCopula package)
derivGumbel <- function(s, u, v) {
  a <- -log(u)
  b <- -log(v)
  A <- a^s + b^s
  x <- -1 + s
  y <- (-3 + 1/s)
  
  temp <- 1/(u * v * s^2) * exp(-A^(1/s))
  temp <- (a^(x/y) * A * b^(x/y))^y * temp
  
  temp2 <- s * (-a^s * ((-1 + s)^2 + (-3 + 2 * s) * A^(1/s) + A^(2/s)) + s * (-1 + 
    s + A^(1/s)) * b^s) * log(a)
  temp2 <- c(temp2, (1 - s + (-3 + s) * A^(1/s) + A^(2/s)) * A * log(A))
  temp2 <- c(temp2, s * (s * a^s * (1 + (-1 + s + A^(1/s)) * log(b)) + b^s * (s - 
    ((-1 + s)^2 + (-3 + 2 * s) * A^(1/s) + A^(2/s)) * log(b))))
  
  return(temp * sum(temp2))
}

# Check parameter(s) for a given copula family
family.check <- function(family, par, par2 = 0) {
  if (!(family %in% c(1, 2, 3, 4, 13, 14, 23, 24, 33, 34))) {
    return("Copula family not yet implemented.")
  } else if ((family == 1 || family == 2) && abs(par) >= 1) {
    return("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
  } else if (family == 2 && par2 <= 2) {
    return("The degrees of freedom parameter of the t-copula has to be larger than 2.")
  } else if ((family == 3 || family == 13) && par <= 0) {
    return("The parameter of the Clayton copula has to be positive.")
  } else if ((family == 4 || family == 14) && par < 1) {
    return("The parameter of the Gumbel copula has to be in the interval [1,oo).")
  } else if ((family == 23 || family == 33) && par >= 0) {
    return("The parameter of the rotated Clayton copula has to be negative.")
  } else if ((family == 24 || family == 34) && par > -1) {
    return("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
  } else {
    return(TRUE)
  }
}

createMaxMat <- function(Matrix){
  
  if(dim(Matrix)[1]!=dim(Matrix)[2]) 
    stop("Structure matrix has to be quadratic.")
  
  MaxMat <- reorderRVineMatrix(Matrix)  
  n  <- nrow(MaxMat)  
  for(j in 1:(n-1)){
    for(i in (n-1):j){
      MaxMat[i,j] <- max(MaxMat[i:(i+1),j])
    }
  }
  
  tMaxMat <- MaxMat
  tMaxMat[is.na(tMaxMat)] <- 0  
  oldSort <- diag(Matrix)
  oldSort <- oldSort[n:1]  
  for(i in 1:n){
    MaxMat[tMaxMat == i] <- oldSort[i]
  }
  
  return(MaxMat)
}

neededCondDistr <- function(Vine){
  if(dim(Vine)[1]!=dim(Vine)[2]) stop("Structure matrix has to be quadratic.")
  
  Vine <- reorderRVineMatrix(Vine)
  MaxMat <- createMaxMat(Vine)
  d <- nrow(Vine)
  
  M <- list()
  M$direct <- matrix(FALSE,d,d)
  M$indirect <- matrix(FALSE,d,d)  
  M$direct[2:d,1] <- TRUE
  
  for(i in 2:(d-1)){
    v <- d-i+1 
    bw <- as.matrix(MaxMat[i:d,1:(i-1)]) == v  
    direct <- Vine[i:d,1:(i-1)] == v  
    M$indirect[i:d,i] <- apply(as.matrix(bw & (!direct)),1,any)   
    M$direct[i:d,i] <- TRUE  
    M$direct[i,i] <- any(as.matrix(bw)[1,] & as.matrix(direct)[1,])
  }
  
  return(M)
}

reorderRVineMatrix <- function(Matrix){
  oldOrder <- diag(Matrix) 
  O <- apply(t(1:nrow(Matrix)),2,"==", Matrix) 
  for(i in 1:nrow(Matrix)){
    Matrix[O[,oldOrder[i]]] <- nrow(Matrix)-i+1
  }  
  return(Matrix)
}

# 
# ## mclapply.hack Nathan VanHoudnos nathanvan AT northwestern FULL STOP edu July
# ## 14, 2014 A script to implement a hackish version of parallel:mclapply() on
# ## Windows machines.  On Linux or Mac, the script has no effect.
# 
# ## Define the hack
# mclapply.hack <- function(...) {
#   ## Create a cluster
#   size.of.list <- length(list(...)[[1]])
#   cl <- parallel::makeCluster(min(size.of.list, parallel::detectCores()))
#   
#   ## Find out the names of the loaded packages
#   loaded.package.names <- c(
#     ## Base packages
#     sessionInfo()$basePkgs,
#     ## Additional packages
#     names( sessionInfo()$otherPkgs ))
#   
#   tryCatch({
#     
#     ## Copy over all of the objects within scope to all clusters.
#     this.env <- environment()
#     while (identical(this.env, globalenv()) == FALSE) {
#       parallel::clusterExport(cl, ls(all.names = TRUE, env = this.env), envir = this.env)
#       this.env <- parent.env(environment())
#     }
#     parallel::clusterExport(cl, ls(all.names = TRUE, env = globalenv()), envir = globalenv())
#     
#     ## Load the libraries on all the clusters N.B. length(cl) returns the number of
#     ## clusters
#     parallel::parLapply(cl, 1:length(cl), function(xx) {
#       lapply(loaded.package.names, function(yy) {
#         require(yy, character.only = TRUE)
#       })
#     })
#     
#     ## Run the lapply in parallel
#     return(parallel::parLapply(cl, ...))
#   }, finally = {
#     ## Stop the cluster
#     parallel::stopCluster(cl)
#   })
# }
# 
# 
# ## If the OS is Windows, set mclapply to the hackish version. Otherwise, leave
# ## the definition alone.
# mclapply <- switch(Sys.info()[["sysname"]], Windows = {
#   mclapply.hack
# }, Linux = {
#   parallel::mclapply
# }, Darwin = {
#   parallel::mclapply
# })