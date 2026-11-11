## In a non-simplified vine, the conditioning variables of a pair beyond the
## first tree are covariates of its model, named after the data's columns, and
## gamBiCopFit reads them from a data frame that names the pair's arguments u1
## and u2. A variable named u1 or u2 was read as the pair's argument, without
## error: such names are refused with a clear error. A name starting with xt
## was partially matched to mgcv's xt and broke the fit: it works.
library(gamCopula)
library(VineCopula)

## A D-vine on 1 - 2 - 3 whose pair (1,3 | 2) is a Clayton copula with
## Kendall's tau 0.05 + 0.5 * u2, the others Gaussian with correlation 0.85.
set.seed(1)
n <- 600
w <- matrix(runif(3 * n), n, 3)
b <- BiCopHinv1(w[, 1], w[, 3], 3, BiCopTau2Par(3, 0.05 + 0.5 * w[, 2]))
u <- cbind(
  BiCopHinv2(w[, 1], w[, 2], 1, 0.85), w[, 2],
  BiCopHinv2(b, w[, 2], 1, 0.85)
)
Matrix <- matrix(c(1, 3, 2, 0, 2, 3, 0, 0, 3), 3, 3)
fit <- function(names, ...) {
  colnames(u) <- names
  gamVineCopSelect(u, Matrix, simplified = FALSE, familyset = c(1, 301), ...)
}
error <- function(expr) {
  tryCatch({
    expr
    ""
  }, error = conditionMessage)
}

## Variable 2, the conditioning variable, cannot be named u2; variable 1, which
## conditions no pair, can be named u1.
stopifnot(grepl("rename: u2.$", error(fit(c("u1", "u2", "u3")))))

## Nor can a covariate be named z, mgcv's working response.
stopifnot(grepl("rename: z.$", error(fit(c("a", "b", "c"),
  smooth.covs = data.frame(z = runif(n))
))))

## A variable named xt2 is the covariate of the pair (1,3 | 2).
m <- fit(c("a", "xt2", "c"))@model[[3]]
stopifnot(identical(all.vars(m@model$pred.formula), "xt2"))
tau <- gamBiCopPredict(m, data.frame(xt2 = c(0.1, 0.9)), target = "tau")$tau
stopifnot(abs(tau - c(0.1, 0.5)) < 0.1)
