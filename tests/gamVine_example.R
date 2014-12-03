#'  set.seed(1)
#'  
#'  ##  Simulation parameters
#'  # Sample size
#'  n <- 2e2  
#'  # Copula family
#'  fam <- 1:4
#'  # Degrees of freedom (for the t-copula when fam <- 2)
#'  par2 <- 4 
#'  # Should the model be specified in terms of Kendall's tau (TRUE) or copula parameter
#'  tau <- FALSE
#'  # Newton-Raphse ("NR) or Fisher-Scoring ("FS") algorithm
#'  met <- "NR"
#'  # Relative tolerance for "NR"/"FS"
#'  tol.rel <- 1e-6
#'  # Max number of iterations for "NR"/"FS"
#'  n.iters <- 25
#'  # Define a 5-dimensional R-vine tree structure matrix
#'  d <- 5
#'  Matrix <- c(5,2,3,1,4,0,2,3,4,1,0,0,3,4,1,0,0,0,4,1,0,0,0,0,1)
#'  Matrix <- matrix(Matrix,d,d)
#'  nnames <- paste("X", 1:5, sep = "")
#'  
#'  ## Covariates distribution
#'  rho <- 0.5
#'  
#'  ## Calibration functions
#'  eta0 <- 1
#'  
#'  calib.quad <- function(t, Ti = 0, Tf = 1, b = 8){
#'    Tm <- (Tf-Ti)/2
#'    a <- -(b/3)*(Tf^2-3*Tf*Tm+3*Tm^2)
#'    return(a + b*(t-Tm)^2)}
#'  
#'  calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1){
#'    a <- b*(1-2*Tf*pi/(f*Tf*pi+cos(2*f*pi*(Tf-Ti))-cos(2*f*pi*Ti)))
#'    return((a+b)/2 + (b-a)*sin(2*f*pi*(t-Ti))/2)}
#'  
#'  calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf/8){
#'    Tm <- (Tf-Ti)/2
#'    a <- (b*s*sqrt(2*pi)/Tf)*(pnorm(0,Tm,s)-pnorm(Tf,Tm,s))
#'    return(a+b*exp(-(t-Tm)^2/(2*s^2)))}
#'  
#'  calib.fct <- c(calib.quad, calib.sin, calib.exp)
#'  
#'  # define gam-vine model list
#'  count <- 1
#'  model <- vector(mode = "list", length = d*(d-1)/2)
#'  par <- matrix(0,d,d)
#'  sel <- seq(d,d^2-d, by = d)
#'  
#'  # First tree
#'  for(i in 1:(d-1))
#'  {
#'    # Select copula family
#'    family <- sample(fam, 1)
#'    model[[count]]$family <- family
#'    
#'    # Simulate and fit some data
#'    temp <- VineCopula::BiCopSim(n, family, par = BiCopEta2Par(family, eta0), par2 = par2)
#'    temp <- BiCopEst(temp[,1],temp[,2],family)    
#'    model[[count]]$par <- temp$par
#'    model[[count]]$par2 <- temp$par2
#'    count <- count + 1
#'  }
#'  
#'  # Trees 2 to (d-1)
#'  for(j in 2:(d-1)){
#'    for(i in 1:(d-j)){ 
#'      # Select a copula family
#'      family <- sample(fam, 1)  
#'      
#'      # Select a true underlying smooth component for each conditioning variable  
#'      cond <- nnames[sort(Matrix[(d-j+2):d,i])]
#'      l <- length(cond)
#'      temp <- sample(3, l, replace = TRUE)
#'      
#'      # Simulate a dataset
#'      if(l != 1){
#'        covariates.distr <- copula::mvdc(copula::normalCopula(rho, dim = l), c("unif"), list(list(min = 0, max = 1)), marginsIdentical = TRUE) 
#'        X <- copula::rMvdc(n,covariates.distr)         
#'        tmp.fct <- paste("function(",paste(cond,collapse=","),"){eta0+",
#'                   paste(sapply(1:l, function(x) paste("calib.fct[[",temp[x],"]](",cond[x],")",
#'                   sep="")), collapse="+"),"}",sep="")
#'        tmp.fct <- eval(parse(text = tmp.fct))
#'      }else{
#'        X <- matrix(runif(n),ncol=1)
#'        tmp.fct <- function(x) eta0+calib.fct[[temp]](x)    
#'      }
#'      colnames(X) <- cond
#'      U <- CondBiCopSim(family, tmp.fct, X, par2 = 4, return.par = FALSE, tau = tau) 
#'      dataset <- data.frame("u1" = U[,1], "u2" = U[,2])
#'      dataset <- cbind(dataset, X)
#'      names(dataset)[3:(2+l)] <- cond   
#'      
#'      # Estimate the gam model
#'      formula.expr <- paste("~",paste(paste("s(",cond,", k=10, bs='cr')", sep = ""),collapse=" + "))
#'      formula.expr <- as.formula(formula.expr)
#'      model[[count]] <- gamBiCopEst(family = family, parFrhs = formula.expr, dataset, tau = tau,
#'                                    method = met, tol.rel = tol.rel, n.iters = n.iters)$res 
#'      count <- count+1  
#'    } 
#'  }
#'  # define gamVineMatrix object
#'  GVC <- gamVine(Matrix=Matrix,model = model,names=nnames)
#'  print(GVC)
#' 
#'  #  # simulate the gamVineMatrix
#'  #  N <- 1e3
#'  #  sim <- gamVineSim(N, GVC)
#'  #  GVC3 <- gamVineSeqEst(sim, GVC, verbose = TRUE)
#'  
