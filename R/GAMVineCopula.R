# #########################
# ##  gam vine copulas  ##
# #########################

# validGAMVineCopula = function(object) {
  # dim <- object@dimension
  # if( dim <= 2)
    # return("Number of dimension too small (>2).")
  # if(length(object@copulas)!=(dim*(dim-1)/2))
    # return("Number of provided copulas does not match given dimension.")
  # if(!any(unlist(lapply(object@copulas,function(x) is(x,"copula")))))
    # return("Not all provided copulas in your list are indeed copulas.")
  # return (TRUE)
# }

# setOldClass("GAMVineMatrix")

# setClass("GAMVineCopula",
         # representation = representation(copulas="list", dimension="integer", 
                                         # GVM="GAMVineMatrix"),
         # prototype = prototype(GVM=structure(list(),class="GAMVineMatrix")),
         # validity = validGAMVineCopula,
         # contains = list("copula")
# )

# # constructor
# GAMVineCopula <- function (GVM) { # GVM <- 4L
  
  # stopifnot(class(GVM)=="GAMVineMatrix") 
  
  # ltr <- lower.tri(GVM$Matrix)
  # copDef <- cbind(GVM$family[ltr], GVM$par[ltr], GVM$par2[ltr])
  # copulas <- rev(apply(copDef,1, function(x) { 
                                   # copulaFromFamilyIndex(x[1],x[2],x[3])
                                 # }))
  
  # new("GAMVineCopula", copulas=copulas, dimension = as.integer(nrow(GVM$Matrix)),
      # GVM=GVM, parameters = numeric(),
      # param.names = character(), param.lowbnd = numeric(), 
      # param.upbnd = numeric(), fullname = paste("GAMVine copula family."))
# }

# showGAMVineCopula <- function(object) {
  # dim <- object@dimension
  # cat(object@fullname, "\n")
  # cat("Dimension: ", dim, "\n")
  # cat("Represented by the following",dim*(dim-1)/2, "copulas:\n")
  # for (i in 1:length(object@copulas)) {
    # cat("  ", class(object@copulas[[i]]), "with parameter(s)", 
        # object@copulas[[i]]@parameters, "\n")
  # }
# }

# setMethod("show", signature("GAMVineCopula"), showGAMVineCopula)

# ## density ##

# dGAMVine <- function(u, copula, log=FALSE) {
  # GVM <- copula@GVM
  # vineLoglik <- GAMVineLogLik(u, GVM, separate=TRUE)$loglik
  # if(log)
    # return(vineLoglik)
  # else
    # return(exp(vineLoglik))
# }

# setMethod("dCopula", signature("numeric","GAMVineCopula"), 
          # function(u, copula, log, ...) {
            # dGAMVine(matrix(u, ncol=copula@dimension), copula, log, ...)
          # })
# setMethod("dCopula", signature("matrix","GAMVineCopula"), dGAMVine)
# setMethod("dCopula", signature("data.frame","GAMVineCopula"), 
          # function(u, copula, log, ...) {
            # dGAMVine(as.matrix(u), copula, log, ...)
          # })

# ## simulation

# GAMVine <- function(n, copula) {
  # GVM <- copula@GVM
  # GAMVineSim(n, GVM)
# }

# setMethod("rCopula", signature("numeric","GAMVineCopula"), GAMVine)

# # fitting using GAMVine
# fitGAMVineCop <- function(copula, data, 
                       # method=list(StructureSelect=FALSE, indeptest=FALSE)) {
  # stopifnot(copula@dimension==ncol(data))
  # if("familyset" %in% names(method))
    # familyset <- method[["familyset"]]
  # else
    # familyset <- NA
  # if("indeptest" %in% names(method))
    # indept <- method[["indeptest"]]
  # else
    # indept <- FALSE
  # if("StructureSelect" %in% names(method)) {
    # if(method[["StructureSelect"]])
      # GAMVineCop <- GAMVineCopula(GAMVineStructureSelect(data, familyset, indeptest=indept))
    # else 
      # GAMVineCop <- GAMVineCopula(GAMVineCopSelect(data, familyset, copula@GVM$Matrix, indeptest=indept))
  # } else {
    # GAMVineCop <- GAMVineCopula(GAMVineCopSelect(data, familyset, copula@GVM$Matrix, indeptest=indept))
  # }
  
  # return(GAMVineCop)
# }

# setMethod("fitCopula", signature=signature("GAMVineCopula"), fitGAMVineCop) 