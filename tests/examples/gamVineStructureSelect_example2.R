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
#' # Pre-allocate the GAM-Vine model list
#' count <- 1
#' model <- vector(mode = "list", length = d*(d-1)/2)
#' 
#' # The first tree contains only the two unconditional copulas
#' for (i in 1:(d-1)) {
#'   model[[count]] <- list(family = fam[count], par = par[count], par2 = 0)
#'   count <- count + 1
#' }
#' 
#' # The second tree contains a unique conditional copula
#' # In this first example, we take a linear calibration function (10*x-5)
#' 
#' # Set-up a dummy dataset
#' tmp <- data.frame(u1 = runif(1e2), u2 = runif(1e2), x1 = runif(1e2))
#' 
#' # Set-up an arbitrary linear model for the calibration function
#' model[[count]] <- gamBiCopEst(tmp, ~ x1, fam[count])$res
#' 
#' # Update the coefficients of the model
#' attr(model[[count]],"model")$coefficients <- c(-5, 10)
#' 
#' # Define gamVine object
#' GVC <- gamVine(Matrix = Matrix, model = model, names = nnames)
#' GVC
#' 
#' # Simulate new data
#' N <- 1e3
#' simData <- data.frame(gamVineSim(N, GVC))
#' colnames(simData) <- nnames
#' 
#' # Fit data
#' summary(fitGVC <- gamVineSeqEst(simData, GVC))
#' summary(fitGVC2 <- gamVineStructureSelect(simData, GVC))