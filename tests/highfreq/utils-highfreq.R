# get.medRV <- function (r, N) {
#   
#   q = abs(as.numeric(r))
#   q = as.numeric(rollmedian(q, k = 3, align = "center"))        
#   fact <- (pi/(6 - 4 * sqrt(3) + pi)) * (N/(N - 2))    
#   medrv <- fact*rollsum(q^2, k = N-2, align = "center")
#   
#   return(medrv)
#   
# }
# 
# get.FFF <- function(r, K, m, w = NA){
#     
#   X <- matrix(0, ncol = 1+2*K, nrow = length(r))  
#   tt <- mod((1:length(r))/m,1)
#   X[,1] <- 1
#   for(j in 1:K){
#     X[,1+j] <- cos(2*pi*j*tt)
#     X[,1+j+K] <- sin(2*pi*j*tt)
#   }  
#   
#   vol <- 2*log(abs(as.numeric(r)))
#   data <- cbind(vol,X)
#   
#   if(is.na(w)){
#     model <- lm.fit(y=data[,1],x=data[,-1])
#     s <- exp(model$fitted.values/2)
#     normalization <- mean(s[1:m])
#     return(s/normalization)
#   } else {
#     rollreg <- function(data) {
#       model <- lm.fit(y=data[,1],x=data[,-1])
#       s <- exp(model$fitted.values/2)
#       normalization <- mean(s[1:m])
#       return(s[w]/normalization)
#     }   
#     return(rollapply(data, width = w, FUN = rollreg, by.column = FALSE))
#   }
# }