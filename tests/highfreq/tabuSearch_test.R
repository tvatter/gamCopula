# library(mgcv)
# library(tabuSearch)
# # http://archive.ics.uci.edu/ml/datasets/Housing
# housing=read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/housing/housing.data")
# colnames(housing)=c("CRIM", "ZN", "INDUS", "CHAS", "NOX", "RM", "AGE", "DIS", "RAD", "TAX", "PTRATIO", "B", "LSTAT", "MEDV")
# 
# housing$CHAS=as.factor(housing$CHAS)
# housing$RAD=as.factor(housing$RAD) # Changed to factor bc only 9 unique values
# summary(housing)
# 
# set.seed(20120823)
# cv=sample(nrow(housing))
# train=housing[cv[1:300],]
# valid=housing[cv[301:400],]
# test=housing[cv[401:nrow(housing)],]
# 
# ssto=sum((valid$MEDV-mean(valid$MEDV))^2)
# evaluate <- function(th){ 
#   num.cols=sum(th)
#   if (num.cols == 0) return(0)
#   fo.str="MEDV ~"
#   cum.cols=0
#   for (i in 1:length(th)) {
#     if (th[i]>0) {
#       if (is.factor(train[,i])) {
#         fo.str=paste(fo.str,colnames(train)[i],sep=" ")
#       } else {
#         fo.str=paste(fo.str," s(",colnames(train)[i],")",sep="")
#       }
#       cum.cols=cum.cols+1
#       if (cum.cols<num.cols) {
#         fo.str=paste(fo.str,"+")
#       }
#     }
#   }
#   #   colnames(train)[c(th,0)==1]
#   fo=as.formula(fo.str)
#   gam1 <- gam(fo,data=train)
#   pred1 <- predict(gam1,valid,se=FALSE)
#   sse <- sum((pred1-valid$MEDV)^2,na.rm=TRUE)
#   return(1-sse/ssto)
# }
# 
# res <- tabuSearch(size = ncol(train)-1, iters = 20, objFunc = evaluate, listSize = 5,
#                   config = rbinom(ncol(train)-1,1,.5), nRestarts = 4,verbose=TRUE)