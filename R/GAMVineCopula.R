# #########################
# ##  gam vine copulas  ##
# #########################

# validgamVineCopula = function(object) {
  # dim <- object@dimension
  # if( dim <= 2)
    # return("Number of dimension too small (>2).")
  # if(length(object@copulas)!=(dim*(dim-1)/2))
    # return("Number of provided copulas does not match given dimension.")
  # if(!any(unlist(lapply(object@copulas,function(x) is(x,"copula")))))
    # return("Not all provided copulas in your list are indeed copulas.")
  # return (TRUE)
# }

# setOldClass("gamVineMatrix")

# setClass("gamVineCopula",
         # representation = representation(copulas="list", dimension="integer", 
                                         # GVM="gamVineMatrix"),
         # prototype = prototype(GVM=structure(list(),class="gamVineMatrix")),
         # validity = validgamVineCopula,
         # contains = list("copula")
# )

# # constructor
# gamVineCopula <- function (GVM) { # GVM <- 4L
  
  # stopifnot(class(GVM)=="gamVineMatrix") 
  
  # ltr <- lower.tri(GVM$Matrix)
  # copDef <- cbind(GVM$family[ltr], GVM$par[ltr], GVM$par2[ltr])
  # copulas <- rev(apply(copDef,1, function(x) { 
                                   # copulaFromFamilyIndex(x[1],x[2],x[3])
                                 # }))
  
  # new("gamVineCopula", copulas=copulas, dimension = as.integer(nrow(GVM$Matrix)),
      # GVM=GVM, parameters = numeric(),
      # param.names = character(), param.lowbnd = numeric(), 
      # param.upbnd = numeric(), fullname = paste("gamVine copula family."))
# }

# showgamVineCopula <- function(object) {
  # dim <- object@dimension
  # cat(object@fullname, "\n")
  # cat("Dimension: ", dim, "\n")
  # cat("Represented by the following",dim*(dim-1)/2, "copulas:\n")
  # for (i in 1:length(object@copulas)) {
    # cat("  ", class(object@copulas[[i]]), "with parameter(s)", 
        # object@copulas[[i]]@parameters, "\n")
  # }
# }

# setMethod("show", signature("gamVineCopula"), showgamVineCopula)

# ## density ##

# dgamVine <- function(u, copula, log=FALSE) {
  # GVM <- copula@GVM
  # vineLoglik <- gamVineLogLik(u, GVM, separate=TRUE)$loglik
  # if(log)
    # return(vineLoglik)
  # else
    # return(exp(vineLoglik))
# }

# setMethod("dCopula", signature("numeric","gamVineCopula"), 
          # function(u, copula, log, ...) {
            # dgamVine(matrix(u, ncol=copula@dimension), copula, log, ...)
          # })
# setMethod("dCopula", signature("matrix","gamVineCopula"), dgamVine)
# setMethod("dCopula", signature("data.frame","gamVineCopula"), 
          # function(u, copula, log, ...) {
            # dgamVine(as.matrix(u), copula, log, ...)
          # })

# ## simulation

# gamVine <- function(n, copula) {
  # GVM <- copula@GVM
  # gamVineSim(n, GVM)
# }

# setMethod("rCopula", signature("numeric","gamVineCopula"), gamVine)

# # fitting using gamVine
# fitgamVineCop <- function(copula, data, 
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
      # gamVineCop <- gamVineCopula(gamVineStructureSelect(data, familyset, indeptest=indept))
    # else 
      # gamVineCop <- gamVineCopula(gamVineCopSelect(data, familyset, copula@GVM$Matrix, indeptest=indept))
  # } else {
    # gamVineCop <- gamVineCopula(gamVineCopSelect(data, familyset, copula@GVM$Matrix, indeptest=indept))
  # }
  
  # return(gamVineCop)
# }

# setMethod("fitCopula", signature=signature("gamVineCopula"), fitgamVineCop) 