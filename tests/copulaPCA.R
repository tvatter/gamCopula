#' library(copula)
#' library(cds)
#' library(fastICA)
#' library(MASS)
#' 
#' n <- 1e2
#' d <- 11
#' f <- 10
#' 
#' Posdef <- function (n, ev = runif(n, 0, 10)) 
#' {
#'   Z <- matrix(ncol=n, rnorm(n^2))
#'   decomp <- qr(Z)
#'   Q <- qr.Q(decomp) 
#'   R <- qr.R(decomp)
#'   d <- diag(R)
#'   ph <- d / abs(d)
#'   O <- Q %*% diag(ph)
#'   Z <- t(O) %*% diag(ev) %*% O
#'   return(Z)
#' }
#' 
#' #pdmat <- Posdef(n=5, ev=5:1)
#' #eigen(pdmat)$val
#' 
#' vc <- Posdef(f)
#' dd <- solve(sqrt(diag(diag(vc))))
#' cc <- dd %*% vc %*% dd
#' diag(cc) <- 1
#' cop <- normalCopula(param = P2p(cc), dim = f, dispstr = "un")
#' u <- cbind(rCopula(n,cop), matrix(runif(n*(d-f)),n,(d-f)))
#' 
#' # cc <- rcormat(d, r = f, sparse.prop = 0.6)$R
#' # # sc <- svd(cc)
#' # # sc$d[abs(sc$d) < 1e-12] <- 1e-12
#' # # eigen(cc)$val
#' # # eigen(sc$u %*% diag(sc$d) %*% t(sc$v))$val
#' # diag(cc) <- 1
#' # cop <- normalCopula(param = P2p(cc), dim = d, dispstr = "un")
#' # u <- rCopula(n,cop)
#'  
#' logit <- function(u) log(u/(1-u))
#' ilogit <- function(x) exp(x)/(1+exp(x))
#' # uu <- seq(0,1,1e-3)
#' # xx <- seq(-10,10,1e-2)
#' # par(mfrow=c(1,2))
#' # curve(logit,0,1,1e3)
#' # curve(qnorm,0,1,1e3,add=T, col = "red")
#' # curve(ilogit,-10,10,1e3)
#' # curve(pnorm,-10,10,1e3,add=T, col = "red")
#' 
#' methods <- c("pearson", "spearman", "kendall")
#' links <- list(probit = list(qnorm,pnorm), logit = list(logit,ilogit))
#' 
#' nn <- apply(expand.grid(names(links), methods), 1, function(x) 
#'   paste(x[2], x[1], sep = "/"))
#' 
#' doOne <- function(link, method, u) {
#'   x <- link[[1]](u)
#'   pca <- eigen(cor(x, method = method),T)
#'   p <- pca$values/sum(pca$values)
#'   xt <- x %*% pca$vectors
#'   ut <- link[[2]](xt)
#'   return(list(x = x, loadings = pca$vectors, xt = xt, ut = ut, p = p))
#' }
#' 
#' res <- apply(expand.grid(1:length(links),1:length(methods)), 1, function(j) 
#'   doOne(links[[j[1]]], methods[j[2]],u))
#' 
#' p <- sapply(res, function(x) x$p)
#' ut <- sapply(res, function(x) head(x$ut[,1]))
#' xt <- sapply(res, function(x) head(x$xt[,1]))
#' colnames(p) <- colnames(ut) <- colnames(xt) <- nn
#' 
#' p
#' apply(p,2,cumsum)
#' ut
#' xt
#' 
