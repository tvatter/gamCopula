var2char <- function(var) {
  deparse(substitute(var))
}

valid.logical <- function(x) {
  !is.null(x) &&  length(x) == 1 && !is.na(x) &&
    (is.logical(x) || (x == 0) || (x == 1))
}

msg.logical <- function(x) {
  paste("'", x, "' should take 0/1 or FALSE/TRUE.", sep = "")
}

valid.unif <- function(x) {
  !is.null(x) &&  length(x) == 1 && !is.na(x) && 
    is.numeric(x) && x >= 0 && x <= 1
}

msg.unif <- function(x) {
  paste("'", x, "' should be a real number in [0,1].", sep = "")
}

valid.real <- function(x) {
  !is.null(x) && is.numeric(x) &&  length(x) == 1 && !is.na(x)
}

msg.real <- function(x) {
  paste("'", x, "' should be a real number.", sep = "")
}

valid.posint <- function(x) {
  !is.null(x) &&  length(x) == 1 && !is.na(x) && 
    is.numeric(x) &&  as.integer(x) ==  x && x > 0
}

msg.posint <- function(x) {
  paste("'", x, "' should be a positive integer.", sep = "")
}

valid.family <- function(x) { 
  valid.posint(x) && is.element(x, get.familyset())
}

msg.family <- function(x) {
  paste("Copula family not implemented. '", x, 
        "' should be in {", 
        paste(get.familyset(), collapse = ","), "}.", sep = "")
}

valid.familypos <- function(x, tau) {
  tau > 0 && !is.element(x, c(23, 24, 33, 34))
}

msg.familypos <- function(x) {
  "This copula family cannot be used for positively dependent data."
}

valid.familyneg <- function(x, tau) {
  tau < 0 && !is.element(x, c(3, 4, 13, 14))
}

msg.familyneg <- function(x) {
  "This copula family cannot be used for negatively dependent data."
}

valid.familyset <- function(x) {
  !is.null(x) && ((length(x) == 1 && (is.na(x) || valid.family(x))) ||
                    all(sapply(x,valid.family)))
}

msg.familyset <- function(x) {
  paste("'", x, "' should be either NA or a vector with elements in {", 
        paste(get.familyset(), collapse = ","), "}.", sep = "")
}

valid.familysetpos <- function(x, tau) {
  any(sapply(x, function(y) valid.familypos(y,tau)))
}

msg.familysetpos <- function(x) {
  paste("'", x, "' needs at least ",
        "one bivariate copula family for positive dependence.", sep = "")
}

valid.familysetneg <- function(x, tau) {
  any(sapply(x, function(y) valid.familyneg(y,tau)))
}

msg.familysetneg <- function(x) {
  paste("'", x, "' needs at least ",
        "one bivariate copula family for negative dependence.", sep = "")
}

valid.covariates <- function(x, msg) {
  if (!is.vector(x) || any(class(x) != "character")) {
    return(msg)
  }
  if (!(length(x) == 1 && is.na(x))) {
    l <- length(x)
  } else {
    l <- 0
  }
  return(l)
}