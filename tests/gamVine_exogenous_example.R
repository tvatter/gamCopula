rm(list=ls())
require(mgcv)
set.seed(0)

##  Simulation parameters
# Sample size and the covariate
n <- 1e3
t <- seq(0,1,length.out=n)
# Copula families
familyset <- c(1:2,301:304,401:404)
# Define a 3-dimensional R-vine tree structure matrix
d <- 3
Matrix <- c(2,3,1,0,3,1,0,0,1)
Matrix <- matrix(Matrix,d,d)
nnames <- paste("X", 1:d, sep = "")

## A function factory
eta0 <- 1
calib.surf <- list(
  calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
    Tm <- (Tf - Ti)/2
    a <- -(b/3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
    return(a + b * (t - Tm)^2)},
  calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1) {
    a <- b * (1 - 2 * Tf * pi/(f * Tf * pi +
                                 cos(2 * f * pi * (Tf - Ti))
                               - cos(2 * f * pi * Ti)))
    return((a + b)/2 + (b - a) * sin(2 * f * pi * (t - Ti))/2)},
  calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf/8) {
    Tm <- (Tf - Ti)/2
    a <- (b * s * sqrt(2 * pi)/Tf) * (pnorm(0, Tm, s) - pnorm(Tf, Tm, s))
    return(a + b * exp(-(t - Tm)^2/(2 * s^2)))})

##  Create the model
# Define gam-vine model list
count <- 1
model <- vector(mode = "list", length = d*(d-1)/2)
sel <- seq(d,d^2-d, by = d)

# A dummy dataset
data <- data.frame(u1 = runif(1e2), u2 = runif(1e2), t = t)
tmpform <- "~s(t,k=10,bs='cr')"

# Trees 1 to (d-1)
for(j in 1:(d-1)){
  for(i in 1:(d-j)){ 
    # Select a copula family and a deterministric function
    family <- sample(familyset, 1)  
    tmp <- sample(3, 1)
    fct <- function(x) {
      eta0 + calib.surf[[tmp]](x)
    }
    
    # Spline approximation of the true function
    y <- fct(t)
    form <- as.formula(paste("y", tmpform, collapse = ""))
    b <- gam(form)
    #plot(x[,1],(y-fitted(b))/y)
    
    # Create a dummy gamBiCop object
    tmp <- gamBiCopEst(data = data, formula = form, family = 1, n.iters = 1)$res
    
    # Update the copula family and the model coefficients
    attr(tmp, "model") <- b
    attr(tmp, "family") <- family
    if (family == 2) {
      attr(tmp, "par2") <- 2+exp(rnorm(1))
    }
    model[[count]] <- tmp
    count <- count+1  
  } 
}

# Create the gamVineCopula object
GVC <- gamVine(Matrix=Matrix,model = model,names=nnames,covariates = "t")
print(GVC)

## Simulate and fit the model
sim <- gamVineSim(n, GVC, newdata = data.frame(t=t))
fitGVC <- gamVineSeqEst(cbind(sim,t), GVC, "t", verbose = TRUE)
fitGVC2 <- gamVineCopSelect(cbind(sim,t), Matrix, "t", verbose = TRUE)
#fitGVC3 <- gamVineStructureSelect(sim, verbose = TRUE)

## Plot the results
dev.off()
par(mfrow=c(3,3))
plot(GVC, ylim = c(-2.5,2.5))

plot(fitGVC, ylim = c(-2.5,2.5))

plot(fitGVC2, ylim = c(-2.5,2.5))

#plot(fitGVC3, ylim = c(-2.5,2.5))




