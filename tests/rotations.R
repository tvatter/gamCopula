## A vine's pair copulas must be evaluated as they were fitted. With its default
## rotations = TRUE, gamVineCopSelect stored a double Clayton or Gumbel family
## with its 90 and 270 degrees rotations exchanged, the family of the pair with
## its arguments swapped, and gamVineStructureSelect fitted each next tree on
## the h-functions of that swapped family.
library(gamCopula)
library(VineCopula)

## A D-vine on 1 - 2 - 3 whose pairs are negatively dependent and asymmetric:
## (1,2) and (1,3 | 2) are 90 degrees Clayton, (2,3) a 270 degrees Gumbel.
Matrix <- matrix(c(1, 3, 2, 0, 2, 3, 0, 0, 3), 3, 3)
family <- matrix(c(0, 23, 23, 0, 0, 34, 0, 0, 0), 3, 3)
par <- matrix(c(0, -1.5, -3, 0, 0, -2, 0, 0, 0), 3, 3)
RVM <- RVineMatrix(Matrix, family, par)

set.seed(1)
u <- RVineSim(300, RVM)
colnames(u) <- paste0("v", 1:3)
x <- data.frame(x = runif(300))
u.test <- RVineSim(500, RVM)
colnames(u.test) <- paste0("v", 1:3)
x.test <- data.frame(x = runif(500))

## The Kullback-Leibler divergence of a fit from the true copula, on the test
## sample.
kl <- function(fit) {
  mean(log(RVinePDF(u.test, RVM)) -
    log(gamVinePDF(fit, cbind(u.test, x.test))))
}
select <- function(Matrix = NULL, level = 0.05) {
  if (is.null(Matrix)) {
    gamVineStructureSelect(u,
      smooth.covs = x, simplified = TRUE,
      familyset = c(301, 401), level = level
    )
  } else {
    gamVineCopSelect(u, Matrix,
      smooth.covs = x, simplified = TRUE,
      familyset = c(301, 401), level = level
    )
  }
}

## Both fits are close to the truth (an exchanged rotation costs about 1 in
## gamVineCopSelect, and h-functions of the exchanged family about 0.07).
stopifnot(kl(select(Matrix)) < 0.03)
stopifnot(kl(select()) < 0.03)

## With every pair constant (no covariate passes the test at this level), each
## fit is VineCopula's selection from the same families.
same <- function(fit, RVM) {
  isTRUE(all.equal(
    gamVinePDF(fit, cbind(u.test, x.test)), RVinePDF(u.test, RVM),
    tolerance = 1e-6
  ))
}
families <- c(3, 4, 13, 14, 23, 24, 33, 34)
stopifnot(same(
  select(Matrix, level = 1e-6),
  RVineCopSelect(u, families, Matrix)
))
stopifnot(same(select(level = 1e-6), RVineStructureSelect(u, families)))
