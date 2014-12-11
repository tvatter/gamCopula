## Required functions borrowed from the VineCopula package because
## they are not exported

createMaxMat <- function(Matrix){
  
  if(dim(Matrix)[1]!=dim(Matrix)[2]) 
    stop("Structure matrix has to be quadratic.")
  
  MaxMat <- reorderRVineMatrix(Matrix)  
  n  <- nrow(MaxMat)  
  for(j in 1:(n-1)){
    for(i in (n-1):j){
      MaxMat[i,j] <- max(MaxMat[i:(i+1),j])
    }
  }
  
  tMaxMat <- MaxMat
  tMaxMat[is.na(tMaxMat)] <- 0  
  oldSort <- diag(Matrix)
  oldSort <- oldSort[n:1]  
  for(i in 1:n){
    MaxMat[tMaxMat == i] <- oldSort[i]
  }
  
  return(MaxMat)
}

neededCondDistr <- function(Vine){
  if(dim(Vine)[1]!=dim(Vine)[2]) stop("Structure matrix has to be quadratic.")
  
  Vine <- reorderRVineMatrix(Vine)
  MaxMat <- createMaxMat(Vine)
  d <- nrow(Vine)
  
  M <- list()
  M$direct <- matrix(FALSE,d,d)
  M$indirect <- matrix(FALSE,d,d)  
  M$direct[2:d,1] <- TRUE
  
  for(i in 2:(d-1)){
    v <- d-i+1 
    bw <- as.matrix(MaxMat[i:d,1:(i-1)]) == v  
    direct <- Vine[i:d,1:(i-1)] == v  
    M$indirect[i:d,i] <- apply(as.matrix(bw & (!direct)),1,any)   
    M$direct[i:d,i] <- TRUE  
    M$direct[i,i] <- any(as.matrix(bw)[1,] & as.matrix(direct)[1,])
  }
  
  return(M)
}

reorderRVineMatrix <- function(Matrix){
  oldOrder <- diag(Matrix) 
  O <- apply(t(1:nrow(Matrix)),2,"==", Matrix) 
  for(i in 1:nrow(Matrix)){
    Matrix[O[,oldOrder[i]]] <- nrow(Matrix)-i+1
  }  
  return(Matrix)
}
