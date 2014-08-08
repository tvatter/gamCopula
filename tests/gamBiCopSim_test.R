#' ##       											
#' ##	1) Simulate some data	
#' ##  a) Kendall's tau as a sum of three smooth components
#' ##  b) Distribution of the covariates as a Gaussian copula 
#' ##	2) Estimate the model and simulate new data							
#' 
#' require(gamVineCopula)
#' 
#' set.seed(0)
#' 
#' ##' Simulation parameters
#' # Sample size
#' n <- 2e2  
#' # Correlation between the covariates
#' rho <- 0.5
#' # Copula family (Clayton here)
#' fam <- 3  
#' # Degrees of freedom (for the t-copula when fam <- 2)
#' par2 <- 4 
#' # Should the model be specified in terms of Kendall's tau (TRUE) or copula parameter
#' tau <- FALSE
#' # Newton-Raphse ("NR) or Fisher-Scoring ("FS") algorithm
#' met <- "FS"
#' # Relative tolerance for "NR"/"FS"
#' tol <- 1e-6
#' # Max number of iterations for "NR"/"FS"
#' itermax <- 25
#' 
#' ##
#' ## 1) Simulate some data
#' ##
#' 
#' ## Quadratic calibration
#' Tf <- 1
#' b.quad <- 8*Tf
#' t0.quad <- Tf/2
#' a.quad <- -(b.quad/3)*(Tf^2-3*Tf*t0.quad+3*t0.quad^2)
#' calib.quadratic <- function(t, t0, a, b){return(a + b*(t-t0)^2)}
#' 
#' ## Sinusoidal calibration
#' b.sin <- 1
#' f.sin <- 1
#' t0.sin <- 0
#' a.sin <- b.sin*(1-2*f.sin*Tf*pi/(f.sin*Tf*pi+cos(2*f.sin*pi*(Tf-t0.sin))-cos(2*f.sin*pi*t0.sin)))
#' calib.sinusoidal <- function(t, f, t0, a, b){return((a+b)/2 + (b-a)*sin(2*pi*f*(t-t0))/2)}
#' 
#' ## Exponential calibration 
#' t0.exp  <- Tf/2
#' s.exp <- Tf/8
#' b.exp <- 2
#' a.exp <- (b.exp*s.exp*sqrt(2*pi)/Tf)*(pnorm(0,t0.exp,s.exp)-pnorm(Tf,t0.exp,s.exp))
#' calib.exponential <- function(t, a, b, t0, s){return(a+b*exp(-(t-t0)^2/(2*s^2)))}
#' 
#' ## Base level
#' eta0 <- 1
#' 
#' ## Additive calibration function
#' calib.fct <- function(x1,x2,x3){
#'   return(eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#'          calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin)+
#'          calib.exponential(x3, a.exp, b.exp, t0.exp, s.exp))}
#' 
#' 
#' ## A single dataset
#' covariates.distr <- copula::mvdc(copula::normalCopula(rho, dim = 3), 
#'    c("unif"), list(list(min = 0, max = Tf)), marginsIdentical = TRUE) 
#' X <- copula::rMvdc(n,covariates.distr)
#' temp <- CondBiCopSim(fam, calib.fct, X, par2=par2, return.par=TRUE)
#' calib <- temp$eta
#' param <- temp$par
#' pseudo <- temp$data
#' dataset <- data.frame("u1" = pseudo[,1], "u2" = pseudo[,2], 
#'    "x1" = X[,1], "x2" = X[,2], "x3" = X[,3])
# 
#' ##
#' ## 2) Estimate and simulate model
#' ##
#' 
#' ## Model fit with penalized cubic splines
#' pen <- TRUE
#' basis <- c(3,10,10)
#' formula <- ~s(x1, k=basis[1], bs = "cr", fx= !pen)+
#'   s(x2, k=basis[2], bs = "cr", fx= !pen)+
#'   s(x3, k=basis[3], bs = "cr", fx= !pen)
#' system.time(fit <- gamBiCopEst(dataset, family = fam, parFrhs = formula))
#' 
#' X <- as.data.frame(copula::rMvdc(n,covariates.distr))
#' names(X) <- c("x1", "x2", "x3")
#' sim <- gamBiCopSim(fit$res, X)