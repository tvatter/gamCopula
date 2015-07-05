#' require(VineCopula)
#' set.seed(0)
#' 
#' ## A first example with a 3-dimensional GAM-Vine
#' 
#' # Define a R-vine tree structure matrix
#' d <- 3
#' Matrix <- c(2,3,1,0,3,1,0,0,1)
#' Matrix <- matrix(Matrix,d,d)
#' nnames <- paste("x", 1:d, sep = "")
#' 
#' # Copula families for each edge
#' fam <- c(3,4,1)
#' 
#' # Parameters for the first tree (two unconditional copulas)
#' par <- c(1,2)
#' 
#' # A link for the second tree (a unique conditional copula)
#' g <- function(x){
#'   tanh(x/2)
#' }
#' 
#' # Pre-allocate the GAM-Vine model list
#' count <- 1
#' model <- vector(mode = "list", length = d*(d-1)/2)
#' 
#' # The first tree contains only the two unconditional copulas
#' for(i in 1:(d-1))
#' {
#'   model[[count]] <- list(family = fam[count], par = par[count], par2 = 0)
#'   count <- count + 1
#' }
#' 
#' # The second tree contains a unique conditional copula
#' # In this first example, we take a linear calibration function (10*x-5)
#' tmp <- sapply(seq(0,1,1e-3), function(x) BiCopSim(1, fam[count], g(10*x-5))) 
#' data <- data.frame(u1 = tmp[1,], u2 = tmp[2,], x1 = seq(0,1,1e-3))
#' model[[count]] <- gamBiCopEst(data, ~ x1, fam[count])$res
#' 
#' # Define gamVine object
#' GVC <- gamVine(Matrix = Matrix, model = model, names = nnames)
#' summary(GVC)
#' 
#' # Simulate new data
#' N <- 1e3
#' simData <- data.frame(gamVineSim(N, GVC))
#' colnames(simData) <- nnames
#' 
#' # Fit data
#' summary(fitGVC <- gamVineSeqEst(simData, GVC))
#' 
#' # The second tree contains a unique conditional copula
#' # In this first example, we take a quadratic calibration function
#' quad <- function(t, Ti = 0, Tf = 1, b = 8) {
#'   Tm <- (Tf - Ti)/2
#'   a <- -(b/3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
#'   return(a + b * (t - Tm)^2)}
#' tmp <- sapply(seq(0,1,1e-3), function(x) BiCopSim(1, fam[count], g(quad(x)))) 
#' data <- data.frame(u1 = tmp[1,], u2 = tmp[2,], x1 = seq(0,1,1e-3))
#' model[[count]] <- gamBiCopEst(data, ~ s(x1, k = 5, fx = T), fam[count])$res
#' 
#' # Update the gamVine object
#' GVC <- gamVine(Matrix = Matrix, model = model, names = nnames)
#' summary(GVC)
#' 
#' # Simulate new data
#' N <- 1e3
#' simData <- data.frame(gamVineSim(N, GVC))
#' colnames(simData) <- nnames
#' 
#' # Fit data
#' summary(fitGVC <- gamVineSeqEst(simData, GVC))
