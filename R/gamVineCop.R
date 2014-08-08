#############################
##  gam vine copulas  ##
#############################
valid.gamVineCop = function(object) { 
  
  d <- length(attributes(object))
  if((d < 3) || names(attributes(object))[1:3] != c("Matrix", "model",  "names")){
    return("A gamVineCop contains at least 1) a R-Vine matrix, 2) a list of lists 
            with three numeric items (family, par and par2) and/or gamBiCop objects
            and 3) a vector of names.")
  }
  Matrix <- object@Matrix  
  names <- object@names
  Matrix[is.na(Matrix)]=0  
  d <- dim(Matrix)[1]
  if(dim(Matrix)[2] != d) return("Structure matrix has to be quadratic.")
  if(max(Matrix)>d) return("Error in the structure matrix.")
  if(RVineMatrixCheck(Matrix)!=1) return("'Matrix' is not a valid R-vine matrix")
  if(length(names)>1 & length(names)!= d){ 
    return("Length of the vector 'names' is not correct.")
  }else if(length(names) == 0){
    names <- paste("x", 1:d, sep = "")
  } 
  
  model <- object@model
  if(length(model)!= d*(d-1)/2){ 
    return("Length of the list 'model' is not correct.") 
  }
  count <- 1
  # First tree
  for(i in 1:(d-1)){
    mm <- model[[count]]
    
    if(!is.list(mm) || any(names(mm) != c("family", "par", "par2")) || !is.numeric(unlist(mm))){
      return(paste("Element", count, "of the model list (first tree) should 
          be a list with three numeric items (family, par and par2)."))
    }
    
    chk <- family.check(mm$family, mm$par, mm$par2)
    if(chk != TRUE){
      return(paste("In element", count, "of the model list (first tree):", chk))
    }
    count <- count + 1
  }
 
  # Trees 2 to (d-1)
  for(j in 2:(d-1)){
    for(i in 1:(d-j)){ 
      mm <- model[[count]]
      if(valid.gamBiCop(mm) == TRUE){
        cond <- sort(all.vars(model$pred.formula))
        cond2 <- names[sort(Matrix[(d-j+2):d,i])]
        if(!all(cond == cond2)){
          return(paste("The formula of element", count, "of the model list, (tree", j, ") does not
                     not contain the correct conditioning variables."))}
      }else{
        if(!is.list(mm) || any(names(mm) != c("family", "par", "par2")) || !is.numeric(unlist(mm))){
          return(paste("Element", count, "of the model list, (tree", j, ") should 
                     be a valid gamBiCop object 
                     or a list containing three items (family, par, par2)."))}
        
        chk <- family.check(mm$family, mm$par, mm$par2)
        if(chk != TRUE){
          return(paste("In element", count, "of the model list, (tree", j, "):", chk))}
      } 
      count <- count+1  
    } 
  } 
  return(TRUE)       	
}

#'  The \code{\link{gamVineCop-class}}
#'
#'  \code{\link{gamVineCop-class}} is an S4 class to store 
#'  a generalized additive model on a vine copula.
#'
#' @slot Matrix lower triangular d x d matrix that defines the R-vine tree structure.
#' @slot model list containing d x (d-1)/2 lists with three numeric items (family, 
#' par and par2) and/or \code{\link{gamBiCop-class}} objects.
#' @slot names vector of d names.
#' @seealso \code{\link{gamVineCop}}, \code{\link{RVineMatrix}} and \code{\link{gamBiCop-class}}.
#' @export
setClass("gamVineCop",
         slots = c(Matrix="matrix", model = "list", names = "character"),
         validity = valid.gamVineCop
)

#' Constructor of the \code{\link{gamVineCop-class}}
#'
#'  A constructor for objects of the \code{\link{gamVineCop-class}}.
#'
#' @param Matrix lower triangular d x d matrix that defines the R-vine tree structure.
#' @param model list containing d x (d-1)/2 lists with three numeric items (family, 
#' par and par2) and/or \code{\link{gamBiCop-class}} objects.
#' @param names vector of d names.
#' @return A \code{\link{gamVineCop-class}} object.
#' @seealso \code{\link{gamVineCop-class}}, \code{\link{RVineMatrix}} and \code{\link{gamBiCop-class}}.
gamVineCop  <- 
  function(Matrix,model,names=NA){
	#MaxMat=VineCopula:::createMaxMat(Matrix)
	#CondDistr=VineCopula:::neededCondDistr(Matrix)
  
  new("gamVineCop", Matrix = Matrix, model = model, names = as.character(names))
}
