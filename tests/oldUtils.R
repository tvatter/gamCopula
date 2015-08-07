# ## Fisher information with respect to Kendall's tau for a bivariate copula
# "FisherBiCop2" <- function(family, tau, par2 = NULL) {
#   
#   if (family == 1) {
#     par <- sin(pi/2 * tau)
#     out <- FisherGaussian(par) * (pi * cos(pi * tau/2)/2)^2
#   } else if (family == 2) {
#     par <- sin(pi/2 * tau)
#     out <- FisherStudentRho(par, par2) * (pi * cos(pi * tau/2)/2)^2
#   } else if (family == 3) {
#     par <- 2 * tau/(1 - tau)
#     out <- FisherClayton(par) * (2/(tau - 1)^2)^2
#   } else if (family == 4) {
#     par <- 1/(1 - tau)
#     out <- FisherGumbel(par) * (1/(tau - 1)^2)^2
#   }
#   
#   return(out)
# }
# 
# 
# ## Derivative of the Gumbel copula density with respect to the parameter (Higher
# ## precision than in the VineCopula package)
# derivGumbel <- function(s, u, v) {
#   a <- -log(u)
#   b <- -log(v)
#   A <- a^s + b^s
#   x <- -1 + s
#   y <- (-3 + 1/s)
#   
#   temp <- 1/(u * v * s^2) * exp(-A^(1/s))
#   temp <- (a^(x/y) * A * b^(x/y))^y * temp
#   
#   temp2 <- s * (-a^s * ((-1 + s)^2 + (-3 + 2 * s) * A^(1/s) + A^(2/s)) + 
#                   s * (-1 +  s + A^(1/s)) * b^s) * log(a)
#   temp2 <- c(temp2, (1 - s + (-3 + s) * A^(1/s) + A^(2/s)) * A * log(A))
#   temp2 <- c(temp2, s * (s * a^s * (1 + (-1 + s + A^(1/s)) * log(b)) + 
#                            b^s * (s -  ((-1 + s)^2 + (-3 + 2 * s) * A^(1/s) +
#                                           A^(2/s)) * log(b))))
#   
#   return(temp * sum(temp2))
# }
# 
# 
# getRotations <- function (i) {
#   out <- i
#   if (i %in% c(3, 13, 23, 33)) 
#     out <- c(3, 13, 23, 33)
#   if (i %in% c(4, 14, 24, 34)) 
#     out <- c(4, 14, 24, 34)
#   if (i %in% c(6, 16, 26, 36)) 
#     out <- c(6, 16, 26, 36)
#   if (i %in% c(7, 17, 27, 37)) 
#     out <- c(7, 17, 27, 37)
#   if (i %in% c(8, 18, 28, 38)) 
#     out <- c(8, 18, 28, 38)
#   if (i %in% c(9, 19, 29, 39)) 
#     out <- c(9, 19, 29, 39)
#   if (i %in% c(10, 20, 30, 40)) 
#     out <- c(10, 20, 30, 40)
#   if (i %in% c(104, 114, 124, 134)) 
#     out <- c(104, 114, 124, 134)
#   if (i %in% c(204, 214, 224, 234)) 
#     out <- c(204, 214, 224, 234)
#   out
# }
# 
# 
# 
# ## Rotate data
# rotate.data <- function(data, angle = 0) {
#   
#   valid <- c(0, 90, 180, 270)
#   if (!is.element(angle, valid)) {
#     print("Angle must be 0, 90, 180 or 270.")
#     return(data)
#   }
#   
#   n <- dim(data)
#   if (n[2] != 2) {
#     print("Data must be bivariate.")
#     return(data)
#   }
#   
#   if (angle == 90) {
#     return(cbind(1 - data[, 1], data[, 2]))
#   } else if (angle == 180) {
#     return(cbind(1 - data[, 1], 1 - data[, 2]))
#   } else if (angle == 270) {
#     return(cbind(data[, 1], 1 - data[, 2]))
#   } else {
#     return(data)
#   }
# }
# 
# ## Rotate copula parameter
# rotate.param <- function(param, angle = 0) {
#   
#   #valid <- c(0, 90, 180, 270)
#   #if (!is.element(angle, valid)) {
#   #  print("Angle must be 0, 90, 180 or 270.")
#   #  return(param)
#   #}
#   
#   if (angle == 90) {
#     return(-param)
#   } else if (angle == 180) {
#     return(param)
#   } else if (angle == 270) {
#     return(-param)
#   } else {
#     return(param)
#   }
# }
# 
# 
# withRotations <- function(nums) {
#   unique(unlist(lapply(nums, getRotations)))
# }