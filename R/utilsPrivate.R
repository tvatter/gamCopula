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