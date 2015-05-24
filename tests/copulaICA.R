S <- cbind(sin((1:1000)/20), rep((((1:200)-100)/100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A

a1 <- icafast(X, 1)
a2 <- icafast(X, 2)

# generate noisy data (p!=r)
set.seed(123)
nobs=1000
d <- 10
Amat=cbind(icasamp("a","rnd",nobs),icasamp("b","rnd",nobs))
Bmat=matrix(2*runif(2*nobs),d,2)
Emat=matrix(rnorm(d*nobs),nobs,d)
Xmat=tcrossprod(Amat,Bmat)+Emat

kmax <- d
vafs <- matrix(0,kmax-1,kmax)
for (j in 1:(kmax-1)) {
 vafs[j,1:(j+1)] <- icafast(Xmat,j+1)$vafs
}

kmax <- d
vafs2 <- matrix(0,kmax-1,kmax)
for (j in 1:(kmax-1)) {
  vafs2[j,1:(j+1)] <- sort(rowSums(fastICA(Xmat,j+1)$A^2)*nobs/sum(Xmat^2),
                           decreasing = TRUE)
}

library(microbenchmark)
bm <- microbenchmark(fastICA = fastICA(Xmat,10), icafast = icafast(Xmat,10))

library(copula)
library(ica)
library(MASS)

n <- 1e2
d <- 11
f <- 10

Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

vc <- Posdef(f)
dd <- solve(sqrt(diag(diag(vc))))
cc <- dd %*% vc %*% dd
diag(cc) <- 1
cop <- normalCopula(param = P2p(cc), dim = f, dispstr = "un")
u <- cbind(rCopula(n,cop), matrix(runif(n*(d-f)),n,(d-f)))

logit <- function(u) log(u/(1-u))
ilogit <- function(x) exp(x)/(1+exp(x))
links <- list(probit = list(qnorm,pnorm), logit = list(logit,ilogit))

doOne <- function(link, u) {
  x <- link[[1]](u)
  ica <- icafast(x,dim(x)[2])
  p <- pca$values/sum(pca$values)
  xt <- x %*% pca$vectors
  ut <- link[[2]](xt)
  return(list(x = x, loadings = pca$vectors, xt = xt, ut = ut, p = p))
}

res <- apply(expand.grid(1:length(links),1:length(methods)), 1, function(j) 
  doOne(links[[j[1]]], methods[j[2]],u))

p <- sapply(res, function(x) x$p)
ut <- sapply(res, function(x) head(x$ut[,1]))
xt <- sapply(res, function(x) head(x$xt[,1]))
colnames(p) <- colnames(ut) <- colnames(xt) <- nn

p
apply(p,2,cumsum)
ut
xt


