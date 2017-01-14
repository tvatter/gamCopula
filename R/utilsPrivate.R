fasttau <- function(x, y, weights = NA) {
  m <- length(x)
  n <- length(y)
  if (m == 0 || n == 0) 
    stop("both 'x' and 'y' must be non-empty")
  if (m != n) 
    stop("'x' and 'y' must have the same length")
  ktau <- TauMatrix(matrix(c(x, y), length(x), 2), weights)[2, 1]
  return(ktau)
}

get.modelCount <- function(d) {
  model.count <- rep(0, d^2)
  temp <- 1:(d * (d - 1)/2)
  t1 <- 1
  sel <- seq(d, d^2 - d, by = d)
  for (i in 1:(d - 1)) {
    t2 <- t1 + d - i - 1
    model.count[sel[1:(d - i)] - i + 1] <- temp[t1:t2]
    t1 <- t2 + 1
  }
  model.count <- matrix(model.count, d, d)
}


prepare.data <- function(data, covariates, trunclevel = NA, 
                         familyset = c(1,2,301,401), rotations = TRUE) {
  
  data <- data.frame(data)
  covariates <- as.character(covariates)
  if (!(length(covariates) == 1 && is.na(covariates))) {
    l <- length(covariates)
  } else {
    l <- 0
  }
  n <- dim(data)[1]
  d <- dim(data)[2] - l
  
  if (is.null(colnames(data))) {
    nn <- paste("V",1:d,sep="") 
    if (l != 0) {
      nn <- c(nn,covariates)
    }
    colnames(data) <- nn
  } else {
    nn <- colnames(data)
  }  
  
  if (is.na(trunclevel)) {
    trunclevel <- d
  }
  if (trunclevel == 0) {
    familyset <- 0
  }
  if (length(familyset) == 1 && is.na(familyset)) {
    familyset <- get.familyset()
  }
  if (rotations) {
    familyset <- withRotations(familyset)
  }
  
  return(list(nn = nn, covariates = covariates, data = data,
              l = l, n = n, d = d, 
              trunclevel = trunclevel, familyset = familyset))
}

prepare.data2 <- function(udata, lin.covs, smooth.covs, trunclevel = NA, 
                          familyset = c(1,2,301,401), rotations = TRUE) {
  udata <- as.data.frame(udata)
  if (is.null(colnames(udata)))
    colnames(udata) <- c(paste0("u", 1:ncol(udata)))
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
  
  if (ncol(data) > ncol(udata)) {
    l <- ncol(data) - ncol(udata)
  } else {
    l <- 0
  }
  n <- nrow(udata)
  d <- ncol(udata)
  
  nn <- colnames(data)
  
  if (is.na(trunclevel)) {
    trunclevel <- d
  }
  if (trunclevel == 0) {
    familyset <- 0
  }
  if (length(familyset) == 1 && is.na(familyset)) {
    familyset <- get.familyset()
  }
  if (rotations) {
    familyset <- withRotations(familyset)
  }
  
  return(list(nn = nn, covariates = colnames(data)[-(1:d)], data = data,
              l = l, n = n, d = d, 
              trunclevel = trunclevel, familyset = familyset))
}
