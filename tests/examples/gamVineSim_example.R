#' set.seed(0)
#' 
#' ##  Simulation parameters
#' # Sample size
#' n <- 2e2  
#' # Copula family
#' fam <- 1:4
#' # Degrees of freedom (for the t-copula when fam <- 2)
#' par2 <- 4 
#' # Should the model be specified in terms of Kendall's tau (TRUE) or copula parameter
#' tau <- TRUE
#' # Newton-Raphse ("NR) or Fisher-Scoring ("FS") algorithm
#' met <- "FS"
#' # Relative tolerance for "NR"/"FS"
#' tol.rel <- 1e-6
#' # Max number of iterations for "NR"/"FS"
#' n.iters <- 25
#' # Define a 5-dimensional R-vine tree structure matrix
#' d <- 5
#' Matrix <- c(5,2,3,1,4,0,2,3,4,1,0,0,3,4,1,0,0,0,4,1,0,0,0,0,1)
#' Matrix <- matrix(Matrix,d,d)
#' nnames <- paste("X", 1:5, sep = "")
#' 
#' ## Integration grid
#' step <- 1e-3
#' ngrid <- 1/step
#' xx <- seq(0,1,length.out = ngrid)
#' Xx <- data.frame(cbind(xx,xx,xx,xx,xx))
#' names(Xx) <- nnames
#' true <- true.approx <- matrix(NA, ngrid, dim(Xx)[2])
#' 
#' ## Quadratic calibration
#' Tf <- 1
#' b.quad <- 8*Tf
#' t0.quad <- Tf/2
#' a.quad <- -(b.quad/3)*(Tf^2-3*Tf*t0.quad+3*t0.quad^2)
#' calib.quadratic <- function(t, t0, a, b){return(a + b*(t-t0)^2)}
#' true[,1] <- calib.quadratic(xx, t0.quad, a.quad, b.quad)
#' 
#' ## Sinusoidal calibration
#' b.sin <- 1
#' f.sin <- 1
#' t0.sin <- 0
#' a.sin <- b.sin*(1-2*f.sin*Tf*pi/(f.sin*Tf*pi+cos(2*f.sin*pi*(Tf-t0.sin))-cos(2*f.sin*pi*t0.sin)))
#' calib.sinusoidal <- function(t, f, t0, a, b){return((a+b)/2 + (b-a)*sin(2*pi*f*(t-t0))/2)}
#' true[,2] <- calib.sinusoidal(xx, f.sin, t0.sin, a.sin, b.sin)
#' 
#' ## Exponential calibration 
#' t0.exp  <- Tf/2
#' s.exp <- Tf/8
#' b.exp <- 2
#' a.exp <- (b.exp*s.exp*sqrt(2*pi)/Tf)*(pnorm(0,t0.exp,s.exp)-pnorm(Tf,t0.exp,s.exp))
#' calib.exponential <- function(t, a, b, t0, s){return(a+b*exp(-(t-t0)^2/(2*s^2)))}
#' true[,3] <- calib.exponential(xx, a.exp, b.exp, t0.exp, s.exp)
#' 
#' ## Base level
#' eta0 <- 1
#' 
#' ## Additive calibration function
#' calib.fct <- function(x1,x2,x3){
#'   return(eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#'            calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin)+
#'            calib.exponential(x3, a.exp, b.exp, t0.exp, s.exp))}
#' calib.fct12 <- function(x1,x2){eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#'                                  calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin)}
#' calib.fct13 <- function(x1,x2){eta0+calib.quadratic(x1, t0.quad, a.quad, b.quad)+ 
#'                                  calib.exponential(x2, a.exp, b.exp, t0.exp, s.exp)}
#' calib.fct23 <- function(x1,x2){eta0+calib.sinusoidal(x2, f.sin, t0.sin, a.sin, b.sin) +
#'                                  calib.exponential(x1, a.exp, b.exp, t0.exp, s.exp)}
#' 
#' ## Spline approximation of each true smooth function for a basis size of 10
#' pen <- TRUE
#' basis <- c(5,10,15)
#' for(i in 1:length(basis)){
#'   temp <- rep(0,3)
#'   temp[i] <- 1
#'   temp <- data.frame("y" = calib.fct(xx*temp[1],xx*temp[2],xx*temp[3]), "x" = xx)  
#'   form <- ~s(x, k=basis[i], bs = "cr", fx= !pen)
#'   temp <- gam(update(form,y~.), data = temp)
#'   true.approx[,i] <- predict.gam(temp, type = "terms")
#' }
#' 
#' 
#' # define gam-vine model list
#' count <- 1
#' model <- vector(mode = "list", length = d*(d-1)/2)
#' # First tree
#' for(i in 1:(d-1))
#' {
#'   # Select copula family
#'   family <- sample(fam, 1)
#'   
#'   # Simulate and fit some data
#'   temp <- BiCopSim(n, family, par = BiCopEta2Par(family, eta0), par2 = par2)
#'   dataset <- data.frame("u1" = temp[,1], "u2" = temp[,2])      
#'   model[[count]] <- gamBiCopEst(family = family, parFrhs = ~1, dataset, tau = tau,
#'                                 method = met, tol.rel = tol.rel, n.iters = n.iters)$res
#'   count <- count + 1
#' }
#' 
#' # Trees 2 to (d-1)
#' for(j in 2:(d-1)){
#'   for(i in 1:(d-j)){ 
#'     # Select a copula family
#'     family <- sample(fam, 1)  
#' 
#'     # Select a true underlying smooth component for each conditioning variable  
#'     cond <- nnames[sort(Matrix[(d-j+2):d,i])]
#'     l <- length(cond)
#'     temp <- rep(0,3)
#'     temp[sample(3, l)] <- 1
#'     
#'     # Simulate a dataset
#'     if(l != 1){
#'       rho <- 0
#'       covariates.distr <- copula:::mvdc(copula:::normalCopula(rho, dim = l), c("unif"), list(list(min = 0, max = Tf)), marginsIdentical = TRUE) 
#'       X <- copula:::rMvdc(n,covariates.distr)
#'       formula.expr <- c()
#'       for(o in 1:l){
#'         formula.expr <- c(formula.expr,paste("s(",cond[o], ", k=",basis[which(temp != 0)[o]], ", bs = 'cr')", sep = ""))
#'       }
#'       formula.expr <- paste(formula.expr, collapse = " + ")
#'       formula.expr <-  paste("~", formula.expr, collapse = "", sep = "")     
#'     }else{
#'       X <- runif(n)
#'       formula.expr <- paste("~s(",cond, ", k=", basis[which(temp != 0)], ", bs = 'cr')", sep = "")
#'     }
#'     formula.expr <- as.formula(formula.expr)
#'     
#'     if(l == 1){
#'       switch(which(temp == 1),
#'              calib.temp <- function(x){eta0 + calib.quadratic(x, t0.quad, a.quad, b.quad)},
#'              calib.temp <- function(x){eta0 + calib.sinusoidal(x, f.sin, t0.sin, a.sin, b.sin)},
#'              calib.temp <- function(x){eta0 + calib.exponential(x, a.exp, b.exp, t0.exp, s.exp)})
#'     }else if(l == 2){
#'       switch(paste(which(temp != 0), collapse = ""),
#'              "12"={calib.temp <- function(x1,x2){calib.fct12(x1,x2)}},
#'              "13"={calib.temp <- function(x1,x2){calib.fct13(x1,x2)}},
#'              "23"={calib.temp <- function(x1,x2){calib.fct23(x1,x2)}})
#'     }else{
#'       calib.temp <- function(x1,x2,x3) calib.fct(x1,x2,x3)
#'     }
#'     temp <- CondBiCopSim(family, calib.temp, X, par2=par2, return.par=FALSE, tau = tau)  
#'     dataset <- data.frame("u1" = temp[,1], "u2" = temp[,2])
#'     dataset <- cbind(dataset, X)
#'     names(dataset)[3:(2+l)] <- cond   
#'     
#'     # Estimate the gam model
#'     model[[count]] <- gamBiCopEst(family = family, parFrhs = formula.expr, dataset, tau = tau,
#'                                   method = met, tol.rel = tol.rel, n.iters = n.iters)$res   
#'     count <- count+1  
#'   } 
#' }
#' 
#' # define gamVineMatrix object
#' GVC <- gamVineCop(Matrix=Matrix,model = model,names=nnames)
#' print(GVC)
#' 
#' N <- 5e2
#' system.time(simData <- gamVineSim(N, GVC))
