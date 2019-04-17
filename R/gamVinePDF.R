#' Conditional density function of a gamVine
#'
#' This function returns the density of a conditional pair-copula constructions,
#' where either the copula parameters or the Kendall's taus are modeled as a function
#' of the covariates.
#'
#' @param object \code{\link{gamVine-class}} object.
#' @param data (Same as in \code{\link{predict.gam}} from the
#' \code{\link[mgcv:mgcv-package]{mgcv}} package) A matrix or data frame
#' containing the values of the model covariates at which predictions are
#' required, along with a number of additional columns corresponding to the
#' variables in the pair copula decomposition.
#' @return The conditional density.
#' @examples
#' require(mgcv)
#' set.seed(0)
#'
#' ##  Simulation parameters
#' # Sample size
#' n <- 1e3
#' # Copula families
#' familyset <- c(1:2,301:304,401:404)
#' # Define a 4-dimensional R-vine tree structure matrix
#' d <- 4
#' Matrix <- c(2,3,4,1,0,3,4,1,0,0,4,1,0,0,0,1)
#' Matrix <- matrix(Matrix,d,d)
#' nnames <- paste("X", 1:d, sep = "")
#'
#' ## A function factory
#' eta0 <- 1
#' calib.surf <- list(
#'   calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
#'     Tm <- (Tf - Ti)/2
#'     a <- -(b/3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
#'     return(a + b * (t - Tm)^2)},
#'   calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1) {
#'     a <- b * (1 - 2 * Tf * pi/(f * Tf * pi +
#'                                  cos(2 * f * pi * (Tf - Ti))
#'                                - cos(2 * f * pi * Ti)))
#'     return((a + b)/2 + (b - a) * sin(2 * f * pi * (t - Ti))/2)},
#'   calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf/8) {
#'     Tm <- (Tf - Ti)/2
#'     a <- (b * s * sqrt(2 * pi)/Tf) * (pnorm(0, Tm, s) - pnorm(Tf, Tm, s))
#'     return(a + b * exp(-(t - Tm)^2/(2 * s^2)))})
#'
#' ##  Create the model
#' # Define gam-vine model list
#' count <- 1
#' model <- vector(mode = "list", length = d*(d-1)/2)
#' sel <- seq(d,d^2-d, by = d)
#'
#' # First tree
#' for (i in 1:(d-1)) {
#'   # Select a copula family
#'   family <- sample(familyset, 1)
#'   model[[count]]$family <- family
#'
#'   # Use the canonical link and a randomly generated parameter
#'   if (is.element(family,c(1,2))) {
#'     model[[count]]$par <- tanh(rnorm(1)/2)
#'     if (family == 2) {
#'       model[[count]]$par2 <- 2+exp(rnorm(1))
#'     }
#'   } else {
#'     if (is.element(family,c(401:404))) {
#'       rr <- rnorm(1)
#'       model[[count]]$par <- sign(rr)*(1+abs(rr))
#'     } else {
#'       model[[count]]$par <- rnorm(1)
#'     }
#'     model[[count]]$par2 <- 0
#'   }
#'   count <- count + 1
#' }
#'
#' # A dummy dataset
#' data <- data.frame(u1 = runif(1e2), u2 = runif(1e2), matrix(runif(1e2*d),1e2,d))
#'
#' # Trees 2 to (d-1)
#' for(j in 2:(d-1)){
#'   for(i in 1:(d-j)){
#'     # Select a copula family
#'     family <- sample(familyset, 1)
#'
#'     # Select the conditiong set and create a model formula
#'     cond <- nnames[sort(Matrix[(d-j+2):d,i])]
#'     tmpform <- paste("~",paste(paste("s(", cond, ", k=10, bs='cr')",
#'                                      sep = ""), collapse=" + "))
#'     l <- length(cond)
#'     temp <- sample(3, l, replace = TRUE)
#'
#'     # Spline approximation of the true function
#'     m <- 1e2
#'     x <- matrix(seq(0,1,length.out=m), nrow = m, ncol = 1)
#'     if(l != 1){
#'       tmp.fct <- paste("function(x){eta0+",
#'                        paste(sapply(1:l, function(x)
#'                          paste("calib.surf[[",temp[x],"]](x[",x,"])",
#'                                sep="")), collapse="+"),"}",sep="")
#'       tmp.fct <- eval(parse(text = tmp.fct))
#'       x <- eval(parse(text = paste0("expand.grid(",
#'                                    paste0(rep("x",l), collapse = ","),")",
#'                                    collapse = "")))
#'       y <- apply(x,1,tmp.fct)
#'     }else{
#'       tmp.fct <- function(x) eta0+calib.surf[[temp]](x)
#'       colnames(x) <- cond
#'       y <- tmp.fct(x)
#'     }
#'
#'     # Estimate the gam model
#'     form <- as.formula(paste0("y", tmpform))
#'     dd <- data.frame(y, x)
#'     names(dd) <- c("y", cond)
#'     b <- gam(form, data = dd)
#'     #plot(x[,1],(y-fitted(b))/y)
#'
#'     # Create a dummy gamBiCop object
#'     tmp <- gamBiCopFit(data = data, formula = form, family = 1, n.iters = 1)$res
#'
#'     # Update the copula family and the model coefficients
#'     attr(tmp, "model")$coefficients <- coefficients(b)
#'     attr(tmp, "model")$smooth <- b$smooth
#'     attr(tmp, "family") <- family
#'     if (family == 2) {
#'       attr(tmp, "par2") <- 2+exp(rnorm(1))
#'     }
#'     model[[count]] <- tmp
#'     count <- count+1
#'   }
#' }
#'
#' # Create the gamVineCopula object
#' GVC <- gamVine(Matrix=Matrix,model = model,names=nnames)
#' print(GVC)
#'
#' \dontrun{
#' ## Simulate and fit the model
#' sim <- gamVineSimulate(n, GVC)
#' fitGVC <- gamVineSeqFit(sim, GVC, verbose = TRUE)
#' fitGVC2 <- gamVineCopSelect(sim, Matrix, verbose = TRUE)
#' (gamVinePDF(GVC, sim[1:10, ]))
#'
#' ## Plot the results
#' dev.off()
#' par(mfrow=c(3,4))
#' plot(GVC, ylim = c(-2.5,2.5))
#'
#' plot(fitGVC, ylim = c(-2.5,2.5))
#'
#' plot(fitGVC2, ylim = c(-2.5,2.5))}
#'
#' @seealso \code{\link{gamVine}}, \code{\link{gamVineCopSelect}},
#' \code{\link{gamVineStructureSelect}}, \code{\link{gamVine-class}},
#' \code{\link{gamVineSimulate}} and \code{\link{gamBiCopFit}}.
gamVinePDF <- function(object, data) {
  tmp <- valid.gamVine(object)
  if (tmp != "TRUE") {
    return(tmp)
  }
  covariates <- object@covariates

  ## Transform to dataframe, get dimensions, etc (see in utilsPrivate)
  tmp <- prepare.data(data, covariates)
  n <- tmp$n
  d <- tmp$d
  l <- tmp$l
  nn <- tmp$nn
  data <- tmp$data
  covariates <- tmp$data
  if (l > 0)
    covariates <- cbind(covariates, tmp$covariates)

  oldobject <- object
  oldMat <- object@Matrix
  o <- diag(oldMat)
  oo <- o[length(o):1]
  if (any(o != length(o):1)) {
    object <- gamVineNormalize(object)
    data[, 1:d] <- data[, oo]
  }

  Mat <- object@Matrix
  fam <- gamVineFamily(object)
  MaxMat <- createMaxMat(Mat)
  CondDistr <- neededCondDistr(Mat)

  V <- list()
  V$dens <- array(1, dim = c(d, d, n))
  V$direct <- array(NA, dim = c(d, d, n))
  V$indirect <- array(NA, dim = c(d, d, n))
  V$direct[d, , ] <- t(data[, d:1])

  model.count <- get.modelCount(d)
  for (i in (d - 1):1) {
    for (k in d:(i + 1)) {
      # print(model.count[k, i])
      m <- MaxMat[k, i]
      zr1 <- V$direct[k, i, ]

      if (m == Mat[k, i]) {
        zr2 <- V$direct[k, (d - m + 1), ]
      } else {
        zr2 <- V$indirect[k, (d - m + 1), ]
      }

      mki <- model.count[k, i]
      mm <- object@model[[mki]]
      if (valid.gamBiCop(mm) != TRUE) {
        par <- rep(mm$par, n)
        par2 <- mm$par2
      } else {
        par <- gamBiCopPredict(mm, newdata = covariates, target = "par")$par
        par2 <- mm@par2
      }
      fams <- vapply(
        1:length(par),
        function(j) famTrans(fam[k, i], inv = FALSE, par = par[j]),
        numeric(1)
      )
      
      V$dens[k, i, ] <- BiCopPDF(zr2, zr1, fams, par, par2, check.pars = FALSE)

      if (CondDistr$direct[k - 1, i]) {
        V$direct[k - 1, i, ] <- BiCopHfunc(zr2, zr1,
          fams, par, par2,
          check.pars = FALSE
        )$hfunc1
      }

      if (CondDistr$indirect[k - 1, i]) {
        V$indirect[k - 1, i, ] <- BiCopHfunc(zr2, zr1,
          fams, par, par2,
          check.pars = FALSE
        )$hfunc2
      }
    }
  }

  c(apply(V$dens, 3, prod))
}

# npars.gamVine <- function(object, ...) {
#   sum(sapply(object@model, pair_npar))
# }
# 
# pair_npar <- function(x) {
#   if (gamCopula:::valid.gamBiCop(x) != TRUE) {
#     return(VineCopula::BiCop(x$family, x$par, x$par2)$npars)
#   }
#   l <- gamCopula:::logLik.gamBiCop(x)
#   attributes(l)$df
# }
# 
# AIC.gamVine <- function(object, data) {
#   -2 * sum(log(gamVinePDF(data, object))) + 2 * npars.gamVine(object)
# }
