## gamVineStructureSelect builds each tree's candidate graph from the pairs of
## edges of the previous tree, and skips the pairs that share no variable. The
## selection must not fail when the first candidate pair is skipped (igraph >= 2
## refused the edge attributes it was then given).
library(gamCopula)
library(VineCopula)

## A Gaussian D-vine on the path 1 - 3 - 4 - 2: the first tree's edges come in
## the order (1,3), (2,4), (3,4), so the second tree's first candidate pair,
## ((1,3), (2,4)), shares no variable.
set.seed(1)
d <- 4
Matrix <- matrix(c(1, 2, 4, 3, 0, 2, 3, 4, 0, 0, 3, 4, 0, 0, 0, 4), d, d)
family <- par <- matrix(0, d, d)
family[lower.tri(family)] <- 1
par[d, 1:(d - 1)] <- 0.7
u <- RVineSim(200, RVineMatrix(Matrix, family, par))

fit <- gamVineStructureSelect(u,
  simplified = FALSE, familyset = 1,
  trunclevel = 1
)

## The first tree is the path, so the case above was the one tested.
M <- fit@Matrix
first <- sort(apply(cbind(M[d, 1:(d - 1)], diag(M)[1:(d - 1)]), 1, function(e)
  paste(sort(e), collapse = "-")))
stopifnot(identical(first, c("1-3", "2-4", "3-4")))
stopifnot(RVineMatrixCheck(M) == 1)
