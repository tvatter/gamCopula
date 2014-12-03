# gamVineSeqEst <- function(dataset, GVC, 
#                           method = "NR", tol.rel = 0.001, n.iters = 10, 
#                           verbose = FALSE) {
#   
#   if (!is.matrix(dataset) && !is.data.frame(dataset)) {
#     stop("Dataset has to be either a matrix or a data frame")
#   } else {
#     dataset <- data.frame(dataset)
#     nn <- names(dataset)
#     n <- dim(dataset)[1]
#     d <- dim(dataset)[2]
#   }
#   
#   if (n < 2) 
#     stop("Number of observations has to be at least 2.")
#   if (any(dataset > 1) || any(dataset < 0)) 
#     stop("Data has be in the interval [0,1].")
#   
#   # if(!valid.gamVine(GVC)){ stop('gamBiVineSeqEst can only be used to estimate
#   # from gamVine objects') }
#   
#   # count <- d # First tree for(j in 2:(d-1)){ for(i in 1:(d-j)){ mm <-
#   # model[[count]] if(valid.gamBiCop(mm) == TRUE){ return(paste('Element', count,
#   # 'of the model list, (tree', j, ') should be a valid gamBiCop object or a list
#   # containing three items (family, par, par2).'))} count <- count+1 } }
#   
#   o <- diag(GVC@Matrix)
#   if (length(o) != d) {
#     stop("The dimension of the gamVine object is incorrect.")
#   }
#   
#   oldGVC <- GVC
#   if (any(o != length(o):1)) {
#     GVC <- gamVineNormalize(GVC)
#     dataset <- dataset[, o[length(o):1]]
#   }
#   
#   options(warn = -1)
#   if (!is.na(as.integer(n.iters))) {
#     if ((as.integer(n.iters) < 1) || (as.integer(n.iters) != as.numeric(n.iters))) {
#       stop("N.iters should be a positive integer!")
#     } else {
#       n.iters <- as.integer(n.iters)
#     }
#   } else {
#     stop("N.iters should be a positive integer!")
#   }
#   
#   if (!is.na(as.numeric(tol.rel))) {
#     if ((as.numeric(tol.rel) < 0) || (as.numeric(tol.rel) > 1)) {
#       stop("Tol.rel should be a real number in [0,1]!")
#     } else {
#       tol.rel <- as.numeric(tol.rel)
#     }
#   } else {
#     stop("Tol.rel should be a real number in [0,1]!")
#   }
#   
#   if (!is.element(method, c("FS", "NR"))) {
#     stop("Method should be a string, either FS (Fisher-scoring) or NR (Newton-Raphson)!")
#   }
#   
#   if (!(is.logical(verbose) || (verbose == 0) || (verbose == 1))) {
#     stop("Verbose should takes 0/1 or FALSE/TRUE.")
#   } else {
#     verbose <- as.logical(verbose)
#   }
#   options(warn = 0)
#   
#   Mat <- GVC@Matrix
#   oldMat <- oldGVC@Matrix
#   fam <- gamVineFamily(GVC)
#   MaxMat <- VineCopula:::createMaxMat(Mat)
#   CondDistr <- VineCopula:::neededCondDistr(Mat)
#   
#   V <- list()
#   V$direct <- array(NA, dim = c(d, d, n))
#   V$indirect <- array(NA, dim = c(d, d, n))
#   V$direct[d, , ] <- t(dataset[, d:1])
#   
#   model.count <- rep(0, d^2)
#   temp <- 1:(d * (d - 1)/2)
#   t1 <- 1
#   sel <- seq(d, d^2 - d, by = d)
#   for (i in 1:(d - 1)) {
#     t2 <- t1 + d - i - 1
#     model.count[sel[1:(d - i)] - i + 1] <- temp[t1:t2]
#     t1 <- t2 + 1
#   }
#   model.count <- matrix(model.count, d, d)
#   
#   for (i in (d - 1):1) {
#     for (k in d:(i + 1)) {
#       print(model.count[k, i])
#       m <- MaxMat[k, i]
#       zr1 <- V$direct[k, i, ]
#       
#       if (m == Mat[k, i]) {
#         zr2 <- V$direct[k, (d - m + 1), ]
#       } else {
#         zr2 <- V$indirect[k, (d - m + 1), ]
#       }
#       
#       if (verbose == TRUE) {
#         if (k == d) 
#           message(oldMat[i, i], ",", oldMat[k, i]) else message(oldMat[i, i], ",", oldMat[k, i], "|", paste(oldMat[(k + 
#           1):d, i], collapse = ","))
#       }
#       
#       
#       if (k == d) {
#         temp <- BiCopEst(zr2, zr1, fam[k, i])
#         GVC@model[[model.count[k, i]]]$par <- temp$par
#         GVC@model[[model.count[k, i]]]$par2 <- temp$par2
#         par <- rep(temp$par, n)
#         par2 <- temp$par2
#       } else {
#         cond <- oldMat[(k + 1):d, i]
#         data <- data.frame(u1 = zr2, u2 = zr1)
#         data <- cbind(data, dataset[, cond])
#         names(data)[3:(2 + length(cond))] <- nn[cond]
#         GVC@model[[model.count[k, i]]] <- gamBiCopEst(fam[k, i], GVC@model[[model.count[k, 
#           i]]]@model$formula, data, GVC@model[[model.count[k, i]]]@tau, method, 
#           tol.rel, n.iters)$res
#         par <- gamBiCopPred(GVC@model[[model.count[k, i]]], target = "par")$par
#         par2 <- GVC@model[[model.count[k, i]]]@par2
#       }
#       
#       if (CondDistr$direct[k - 1, i]) {
#         tmp <- rep(0, n)
#         tmp <- sapply(1:n, function(x) .C("Hfunc1", as.integer(fam[k, i]), 
#           as.integer(1), as.double(zr1[x]), as.double(zr2[x]), as.double(par[x]), 
#           as.double(par2), as.double(tmp[x]), PACKAGE = "VineCopula")[[7]])
#         V$direct[k - 1, i, ] <- tmp
#       }
#       
#       if (CondDistr$indirect[k - 1, i]) {
#         tmp <- rep(0, n)
#         tmp <- sapply(1:n, function(x) .C("Hfunc2", as.integer(fam[k, i]), 
#           as.integer(1), as.double(zr2[x]), as.double(zr1[x]), as.double(par[x]), 
#           as.double(par2), as.double(tmp[x]), PACKAGE = "VineCopula")[[7]])
#         V$indirect[k - 1, i, ] <- tmp
#       }
#     }
#   }
#   
#   oldGVC@model <- GVC@model
#   return(oldGVC)
# } 
