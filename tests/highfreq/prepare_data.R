# library(zoo)
# 
# pairs <- c("CHFUSD", "EURUSD", "GBPUSD", "JPYUSD")
# ret <- list()
# for(i in 1:length(pairs)){
#   ret[[i]] <- read.csv(paste("data/ret_",pairs[i],".txt", sep = ""))
# }
# 
# ret <- do.call(cbind.data.frame, ret)[,c(1,seq(2,(length(pairs)*2),2))]
# 
# 
# ret <- zoo(ret[,2:(1+length(pairs))],as.POSIXct(ret[,1], tz = "UTC", "%d-%b-%Y %H:%M:%S"))
# names(ret) <- pairs
# # plot(ret)
# # plot(ret[(which.min(ret[,4])-30):(which.min(ret[,4])+288*2),])
# # plot(ret[(which.min(ret[,1])-30):(which.min(ret[,1])+288*2),])
# save(ret,file = "data/FXreturns.Rda")
# 
