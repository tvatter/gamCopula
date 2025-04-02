#' Structure selection and estimation of a GAM-Vine model.
#'
#' This function select the structure and estimates the parameter(s) of a
#' Generalized Additive model
#' (GAM) Vine model, where GAMs for individual edges are specified either for
#' the copula parameter or Kendall's tau.
#' It solves the maximum penalized likelihood estimation for the copula families
#' supported in this package by reformulating each Newton-Raphson iteration as
#' a generalized ridge regression, which is solved using
#' the \code{\link[mgcv:mgcv-package]{mgcv}} package.
#'
#' @param udata A matrix or data frame containing the data in [0,1]^d.
#' @param lin.covs A matrix or data frame containing the parametric (i.e.,
#' linear) covariates (default: \code{lin.covs = NULL}).
#' @param smooth.covs A matrix or data frame containing the non-parametric
#' (i.e., smooth) covariates (default: \code{smooth.covs = NULL}).
#' @param simplified If \code{TRUE} (default), then a simplified vine is fitted
#' (which is possible only if there are exogenous covariates). If \code{FALSE},
#' then a non-simplified vine is fitted.
#' @param type \code{type = 0} (default) for a R-Vine and \code{type = 1} for a
#' C-Vine.
#' @param familyset An integer vector of pair-copula families to select from
#' (the independence copula MUST NOT be specified in this vector unless one
#' wants to fit an independence vine!).
#' Not listed copula families might be included to better handle
#' limit cases. If \code{familyset = NA} (default), selection among all
#' possible families is performed. Coding of pair-copula families:
#' \code{1} Gaussian,
#' \code{2} Student t,
#' \code{5} Frank,
#' \code{301} Double Clayton type I (standard and rotated 90 degrees),
#' \code{302} Double Clayton type II (standard and rotated 270 degrees),
#' \code{303} Double Clayton type III (survival and rotated 90 degrees),
#' \code{304} Double Clayton type IV (survival and rotated 270 degrees),
#' \code{401} Double Gumbel type I (standard and rotated 90 degrees),
#' \code{402} Double Gumbel type II (standard and rotated 270 degrees),
#' \code{403} Double Gumbel type III (survival and rotated 90 degrees),
#' \code{404} Double Gumbel type IV (survival and rotated 270 degrees).
#' @param rotations If \code{TRUE}, all rotations of the families in familyset
#' are included.
#' @param familycrit Character indicating the criterion for bivariate copula
#' selection. Possible choices: \code{familycrit = 'AIC'} (default) or
#' \code{'BIC'}, as in \code{\link[VineCopula]{BiCopSelect}} from the
#' \code{\link[VineCopula:VineCopula-package]{VineCopula}} package.
#' @param treecrit Character indicating how pairs are selected in each tree.
#' \code{treecrit = "tau"} uses the maximum spanning tree of the Kendall's tau
#' (i.e., the tree of maximal overall dependence), \code{treecrit = "rho"}
#' uses the Spearman's rho.
#' @param level Numerical; Passed to \code{\link{gamBiCopSelect}}, it is the
#' significance level of the test for removing individual
#' predictors (default: \code{level = 0.05}) for each conditional pair-copula.
#' @param trunclevel Integer; level of truncation.
#' @param tau \code{TRUE} (default) for a calibration function specified for
#' Kendall's tau or \code{FALSE} for a calibration function specified
#' for the Copula parameter.
#' @param method \code{'NR'} for Newton-Raphson
#' and  \code{'FS'} for Fisher-scoring (default).
#' @param tol.rel Relative tolerance for \code{'FS'}/\code{'NR'} algorithm.
#' @param n.iters Maximal number of iterations for
#' \code{'FS'}/\code{'NR'} algorithm.
#' @param parallel \code{TRUE} (default) for parallel selection of copula
#' family at each edge or \code{FALSE} for the sequential version.
#' for the Copula parameter.
#' @param verbose \code{TRUE} if informations should be printed during the
#' estimation and \code{FALSE} (default) for a silent version.
#' from \code{\link[mgcv:mgcv-package]{mgcv}}.
#' @param select.once if \code{TRUE} the GAM structure is only selected once,
#'   for the family that appears first in \code{familyset}.
#' @return \code{gamVineSeqFit} returns a \code{\link{gamVine-class}} object.
#' @examples
#' require(VineCopula)
#' set.seed(0)
#'
#'
#' ## An example with a 3-dimensional GAM-Vine
#'
#' # Sample size
#' n <- 1e3
#'
#' # Define a R-vine tree structure matrix
#' d <- 3
#' Matrix <- c(2, 3, 1, 0, 3, 1, 0, 0, 1)
#' Matrix <- matrix(Matrix, d, d)
#' nnames <- paste("x", 1:d, sep = "")
#'
#' # Copula families for each edge
#' fam <- c(301, 401, 1)
#'
#' # Parameters for the first tree (two unconditional copulas)
#' par <- c(1, 2)
#'
#' # Pre-allocate the GAM-Vine model list
#' count <- 1
#' model <- vector(mode = "list", length = d * (d - 1) / 2)
#'
#' # The first tree contains only the two unconditional copulas
#' for (i in 1:(d - 1)) {
#'   model[[count]] <- list(family = fam[count], par = par[count], par2 = 0)
#'   count <- count + 1
#' }
#'
#' # The second tree contains a unique conditional copula
#' # In this first example, we take a linear calibration function (10*x-5)
#'
#' # Set-up a dummy dataset
#' tmp <- data.frame(u1 = runif(1e2), u2 = runif(1e2), x1 = runif(1e2))
#'
#' # Set-up an arbitrary linear model for the calibration function
#' model[[count]] <- gamBiCopFit(tmp, ~x1, fam[count])$res
#'
#' # Update the coefficients of the model
#' attr(model[[count]], "model")$coefficients <- c(-5, 10)
#'
#' # Define gamVine object
#' GVC <- gamVine(Matrix = Matrix, model = model, names = nnames)
#' GVC
#' \dontrun{
#' # Simulate new data
#' simData <- data.frame(gamVineSimulate(n, GVC))
#' colnames(simData) <- nnames
#'
#' # Fit data using sequential estimation assuming true model known
#' summary(fitGVC <- gamVineSeqFit(simData, GVC))
#'
#' # Fit data using structure selection and sequential estimation
#' summary(fitGVC2 <- gamVineStructureSelect(simData, simplified = FALSE))
#' }
#'
#' @seealso  \code{\link{gamVineSeqFit}},\code{\link{gamVineCopSelect}},
#'  \code{\link{gamVine-class}}, \code{\link{gamVineSimulate}} and
#'  \code{\link{gamBiCopSelect}}.
gamVineStructureSelect <- function(udata, lin.covs = NULL, smooth.covs = NULL,
                                   simplified = TRUE, type = 0, familyset = NA,
                                   rotations = TRUE, familycrit = "AIC",
                                   treecrit = "tau", level = 0.05,
                                   trunclevel = NA, tau = TRUE, method = "FS",
                                   tol.rel = 0.001, n.iters = 10,
                                   parallel = FALSE, verbose = FALSE,
                                   select.once = TRUE) {
  tmp <- valid.gamVineStructureSelect(
    udata, lin.covs, smooth.covs, simplified,
    type, familyset, rotations, familycrit,
    treecrit, level, trunclevel,
    tau, method, tol.rel, n.iters,
    parallel, verbose, select.once
  )
  if (tmp != TRUE) {
    stop(tmp)
  }

  ## Transform to dataframe, get dimensions, etc (see in utilsPrivate)
  tmp <- prepare.data2(
    udata, lin.covs, smooth.covs,
    trunclevel, familyset, rotations
  )
  n <- tmp$n
  d <- tmp$d
  l <- tmp$l
  nn <- tmp$nn
  data <- tmp$data
  covariates <- tmp$covariates
  trunclevel <- tmp$trunclevel
  familyset <- tmp$familyset

  if (type == 0) {
    type <- "RVine"
  } else {
    type <- "CVine"
  }

  out <- list(Tree = vector("list", d - 1), Graph = vector("list", d - 1))
  # browser()
  if (is.null(colnames(udata))) {
    colnames(udata) <- paste0("V", 1:ncol(udata))
  }
  graph <- initFirstGraph(udata, treecrit)
  mst <- findMST(graph, type)
  tree <- fitFirstTree(
    mst, udata, lin.covs, smooth.covs, l, covariates,
    familyset, familycrit, treecrit, level,
    tau, method, tol.rel, n.iters, parallel, verbose
  )
  out$Tree[[1]] <- tree
  out$Graph[[1]] <- graph

  oldtree <- tree
  for (i in 2:(d - 1)) {
    if (trunclevel == i - 1) {
      familyset <- 0
    }
    graph <- buildNextGraph(
      tree, udata, l, lin.covs, smooth.covs,
      simplified, treecrit
    )
    mst <- findMST(graph, type)
    tree <- fitTree(
      mst, tree, udata, l, lin.covs, smooth.covs, simplified,
      familyset, familycrit, treecrit, level,
      tau, method, tol.rel, n.iters, parallel, verbose
    )

    out$Tree[[i]] <- tree
    out$Graph[[i]] <- graph
  }

  return(as.GVC(out, covariates))
}

initFirstGraph <- function(data, treecrit) {
  if (treecrit == "tau") {
    C <- TauMatrix(data)
  }
  if (treecrit == "rho") {
    C <- cor(data, method = "spearman")
  }

  rownames(C) <- colnames(C) <- colnames(data)

  g <- graph.adjacency(C, mode = "lower", weighted = TRUE, diag = FALSE)
  E(g)$tau <- E(g)$weight
  E(g)$weight <- 1 - abs(E(g)$weight)
  E(g)$name <- paste(get.edgelist(g)[, 1], get.edgelist(g)[, 2], sep = ",")
  for (i in 1:ecount(g)) {
    E(g)$conditionedSet[[i]] <- ends(g, i, FALSE)
  }

  return(g)
}

findMST <- function(g, mode = "RVine") {
  if (mode == "RVine") {
    return(minimum.spanning.tree(g, weights = E(g)$weight))
  } else {
    M <- abs(get.adjacency(g, attr = "weight", sparse = 0))
    root <- which.min(rowSums(M))
    Ecken <- ends(g, 1:ecount(g), FALSE)
    pos <- Ecken[, 2] == root | Ecken[, 1] == root
    return(delete.edges(g, E(g)[!pos]))
  }
}

fitFirstTree <- function(mst, udata, lin.covs, smooth.covs, l, covariates,
                         familyset, familycrit, treecrit, level,
                         tau, method, tol.rel, n.iters, parallel, verbose) {
  d <- ecount(mst)

  parameterForACopula <- vector("list", d)

  for (i in 1:d) {
    a <- ends(mst, i, FALSE)

    if (is.null(V(mst)[a[1]]$name)) {
      E(mst)[i]$CondName1 <- a[1]
    } else {
      E(mst)[i]$CondName1 <- V(mst)[a[1]]$name
    }

    if (is.null(V(mst)[a[2]]$name)) {
      E(mst)[i]$CondName2 <- a[2]
    } else {
      E(mst)[i]$CondName2 <- V(mst)[a[2]]$name
    }

    if (is.null(V(mst)[a[1]]$name) || is.null(V(mst)[a[2]]$name)) {
      E(mst)[i]$Name <- paste(a[1], a[2], sep = " , ")
    } else {
      E(mst)[i]$Name <- paste(V(mst)[a[1]]$name,
        V(mst)[a[2]]$name,
        sep = " , "
      )
    }

    if (l == 0) {
      parameterForACopula[[i]] <- list()
      parameterForACopula[[i]]$zr1 <- udata[, a[1]]
      parameterForACopula[[i]]$zr2 <- udata[, a[2]]
    } else {
      parameterForACopula[[i]] <- list(
        cbind(
          as.numeric(udata[, a[1]]),
          as.numeric(udata[, a[2]])
        ),
        lin.covs,
        smooth.covs
      )
    }
  }

  if (l == 0) {
    outForACopula <- lapply(
      parameterForACopula, wrapper.fitACopula,
      familyset, familycrit, level
    )
  } else {
    outForACopula <- lapply(
      parameterForACopula, wrapper.fitACopula,
      familyset, familycrit, level,
      tau, method, tol.rel, n.iters, parallel
    )
  }
  for (i in 1:d) {
    E(mst)[i]$model <- list(outForACopula[[i]]$model)
    E(mst)[i]$CondData1 <- list(outForACopula[[i]]$CondOn1)
    E(mst)[i]$CondData2 <- list(outForACopula[[i]]$CondOn2)
  }

  return(mst)
}

fitTree <- function(mst, oldVineGraph, udata, l, lin.covs, smooth.covs,
                    simplified, familyset, familycrit, treecrit, level,
                    tau, method, tol.rel, n.iters, parallel, verbose) {
  d <- ecount(mst)

  parameterForACopula <- list()

  for (i in 1:d) {
    con <- ends(mst, i, FALSE)
    tmp <- ends(oldVineGraph, con, FALSE)

    if ((tmp[1, 1] == tmp[2, 1]) || (tmp[1, 2] == tmp[2, 1])) {
      same <- tmp[2, 1]
    } else {
      if ((tmp[1, 1] == tmp[2, 2]) || (tmp[1, 2] == tmp[2, 2])) {
        same <- tmp[2, 2]
      }
    }

    other1 <- tmp[1, tmp[1, ] != same]
    other2 <- tmp[2, tmp[2, ] != same]

    if (tmp[1, 1] == same) {
      zr1 <- E(oldVineGraph)[con[1]]$CondData2
      n1 <- E(oldVineGraph)[con[1]]$CondName2
    } else {
      zr1 <- E(oldVineGraph)[con[1]]$CondData1
      n1 <- E(oldVineGraph)[con[1]]$CondName1
    }

    if (tmp[2, 1] == same) {
      zr2 <- E(oldVineGraph)[con[2]]$CondData2
      n2 <- E(oldVineGraph)[con[2]]$CondName2
    } else {
      zr2 <- E(oldVineGraph)[con[2]]$CondData1
      n2 <- E(oldVineGraph)[con[2]]$CondName1
    }

    if (is.list(zr1)) {
      zr1a <- as.vector(zr1[[1]])
      zr2a <- as.vector(zr2[[1]])
      n1a <- as.vector(n1[[1]])
      n2a <- as.vector(n2[[1]])
    }
    else {
      zr1a <- zr1
      zr2a <- zr2
      n1a <- n1
      n2a <- n2
    }

    if (verbose == TRUE) message(n1a, " + ", n2a, " --> ", E(mst)[i]$name)

    E(mst)[i]$CondName2 <- n1a
    E(mst)[i]$CondName1 <- n2a

    tmp <- E(mst)$conditioningSet[[i]]
    if (simplified) {
      if (l == 0) {
        parameterForACopula[[i]] <- list()
        parameterForACopula[[i]]$zr1 <- as.numeric(zr1a)
        parameterForACopula[[i]]$zr2 <- as.numeric(zr2a)
      } else {
        parameterForACopula[[i]] <- list(
          cbind(
            as.numeric(zr1a),
            as.numeric(zr2a)
          ),
          lin.covs,
          smooth.covs
        )
      }
    } else {
      temp <- as.data.frame(udata[, tmp])
      colnames(temp) <- colnames(udata)[tmp]
      if (!is.null(smooth.covs)) {
        temp <- cbind(temp, smooth.covs)
      }
      parameterForACopula[[i]] <- list(
        cbind(
          as.numeric(zr1a),
          as.numeric(zr2a)
        ),
        lin.covs, temp
      )
    }
  }

  outForACopula <- lapply(
    parameterForACopula, wrapper.fitACopula,
    familyset, familycrit, level,
    tau, method, tol.rel, n.iters, parallel
  )

  for (i in 1:d) {
    E(mst)[i]$model <- list(outForACopula[[i]]$model)
    E(mst)[i]$CondData2 <- list(outForACopula[[i]]$CondOn2)
    E(mst)[i]$CondData1 <- list(outForACopula[[i]]$CondOn1)
  }

  return(mst)
}

buildNextGraph <- function(graph, udata, l, lin.covs, smooth.covs, simplified,
                           treecrit) {
  EL <- get.edgelist(graph)
  d <- ecount(graph)

  g <- graph.full(d)
  V(g)$name <- E(graph)$name
  V(g)$conditionedSet <- E(graph)$conditionedSet

  if (!is.null(E(graph)$conditioningSet)) {
    V(g)$conditioningSet <- E(graph)$conditioningSet
  }

  for (i in 1:ecount(g)) {
    con <- ends(g, i, FALSE)
    tmp <- ends(graph, con, FALSE)

    ok <- FALSE
    if ((tmp[1, 1] == tmp[2, 1]) || (tmp[1, 2] == tmp[2, 1])) {
      ok <- TRUE
      same <- tmp[2, 1]
    } else {
      if ((tmp[1, 1] == tmp[2, 2]) || (tmp[1, 2] == tmp[2, 2])) {
        ok <- TRUE
        same <- tmp[2, 2]
      }
    }

    if (ok) {
      other1 <- tmp[1, tmp[1, ] != same]
      other2 <- tmp[2, tmp[2, ] != same]

      if (tmp[1, 1] == same) {
        zr1 <- E(graph)[con[1]]$CondData2
      } else {
        zr1 <- E(graph)[con[1]]$CondData1
      }

      if (tmp[2, 1] == same) {
        zr2 <- E(graph)[con[2]]$CondData2
      } else {
        zr2 <- E(graph)[con[2]]$CondData1
      }

      if (is.list(zr1)) {
        zr1a <- as.vector(zr1[[1]])
        zr2a <- as.vector(zr2[[1]])
      } else {
        zr1a <- zr1
        zr2a <- zr2
      }
      noNAs <- !(is.na(zr1a) | is.na(zr2a))

      name.node1 <- strsplit(V(g)[con[1]]$name, split = " *[,|] *")[[1]]
      name.node2 <- strsplit(V(g)[con[2]]$name, split = " *[,|] *")[[1]]
      intersection <- c()

      for (j in 1:length(name.node1)) {
        for (k in 1:length(name.node2)) {
          if (name.node1[j] == name.node2[k]) {
            intersection <- c(intersection, name.node1[j])
            name.node1[j] <- ""
            name.node2[k] <- ""
            break
          }
        }
      }

      difference <- c()
      for (j in 1:length(name.node1)) {
        if (name.node1[j] != "") {
          difference <- c(difference, name.node1[j])
        }
      }

      for (j in 1:length(name.node2)) {
        if (name.node2[j] != "") {
          difference <- c(difference, name.node2[j])
        }
      }

      E(g)[i]$name <- paste(paste(difference, collapse = ","),
        paste(intersection, collapse = ","),
        sep = " | "
      )

      if (is.list(V(g)[con[1]]$conditionedSet)) {
        l1 <- c(
          as.vector(V(g)[con[1]]$conditionedSet[[1]]),
          as.vector(V(g)[con[1]]$conditioningSet[[1]])
        )
        l2 <- c(
          as.vector(V(g)[con[2]]$conditionedSet[[1]]),
          as.vector(V(g)[con[2]]$conditioningSet[[1]])
        )
      } else {
        l1 <- c(V(g)[con[1]]$conditionedSet, V(g)[con[1]]$conditioningSet)
        l2 <- c(V(g)[con[2]]$conditionedSet, V(g)[con[2]]$conditioningSet)
      }
      out <- intersectDifference(l1, l2)

      suppressWarnings({
        E(g)$conditionedSet[i] <- list(out$difference)
      })
      suppressWarnings({
        E(g)$conditioningSet[i] <- list(out$intersection)
      })

      if (treecrit == "tau") {
        E(g)[i]$weight <- 1 - abs(fasttau(zr1a[noNAs], zr2a[noNAs]))
      } else if (treecrit == "rho") {
        E(g)[i]$weight <- 1 - abs(cor(zr1a[noNAs], zr2a[noNAs], method = "spearman"))
      }
    }

    E(g)[i]$todel <- !ok
  }

  E(g)$tau <- E(g)$weight
  g <- delete.edges(g, E(g)[E(g)$todel])

  return(g)
}

wrapper.fitACopula <- function(parameterForACopula, ...) {
  if (length(parameterForACopula) == 2) {
    return(fitACopula(parameterForACopula$zr1, parameterForACopula$zr2, ...))
  } else {
    return(fitAGAMCopula(parameterForACopula, ...))
  }
}

intersectDifference <- function(liste1, liste2) {
  out <- list()
  out$intersection <- c()
  out$difference <- c()

  for (j in 1:length(liste1)) {
    for (k in 1:length(liste2)) {
      if (!is.na(liste2[k]) && liste1[j] == liste2[k]) {
        out$intersection <- c(out$intersection, liste1[j])
        liste1[j] <- NA
        liste2[k] <- NA
        break
      }
    }
  }

  for (j in 1:length(liste1)) {
    if (!is.na(liste1[j])) {
      out$difference <- c(out$difference, liste1[j])
    }
  }
  for (j in 1:length(liste2)) {
    if (!is.na(liste2[j])) {
      out$difference <- c(out$difference, liste2[j])
    }
  }

  return(out)
}

fitACopula <- function(u1, u2, familyset = NA, familycrit = "AIC",
                       level = 0.05, rotation = TRUE) {

  ## transform the familyset to codes for the VineCopula package
  fams <- famTrans(familyset, inv = FALSE, set = TRUE)

  ## select family and estimate parameter(s) for the pair copula
  out <- BiCopSelect(u1, u2,
    fams,
    familycrit,
    indeptest = TRUE,
    level = level,
    rotations = FALSE
  )
  fam <- famTrans(out$family,
    inv = TRUE,
    par = cor(u1, u2), familyset = familyset
  )
  if (rotation == TRUE) {
    if (fam %in% c(301, 303, 401, 403)) {
      fam <- fam + 1
    } else if (fam %in% c(302, 304, 402, 404)) {
      fam <- fam - 1
    }
  }
  out$family <- fam

  ## store pseudo-observations for estimation in next tree
  tmp <- bicoppd1d2(cbind(u1, u2, out$par, out$par2), out$family, p = FALSE, h = TRUE)

  ## save the model and pseudo-observations
  model <- list(family = out$family, par = out$par, par2 = out$par2)
  out <- list(model = model, CondOn1 = tmp[2, ], CondOn2 = tmp[1, ])

  return(out)
}

fitAGAMCopula <- function(data, familyset, familycrit, level, tau,
                          method, tol.rel, n.iters, parallel, rotation = TRUE,
                          select.once = TRUE, ...) {
  out <- list()
  u1 <- data[[1]][, 1]
  u2 <- data[[1]][, 2]
  udata <- cbind(u1, u2)
  lin.covs <- data[[2]]
  smooth.covs <- data[[3]]
  data <- do.call(cbind, data[!(sapply(data, is.null))])

  if (length(familyset) == 1 && familyset == 0) {
    ## independence copula

    out$model <- list(family = 0, par = 0, par2 = 0)
    par <- rep(0, length(u1))
    fam <- 0
    par2 <- 0
  } else {
    ## transform the familyset to codes for the VineCopula package
    fams <- famTrans(familyset, inv = FALSE, set = TRUE)

    ## conditional copula
    tmp <- gamBiCopSelect(udata, lin.covs, smooth.covs,
      familyset, FALSE, familycrit, level, 1.5, tau,
      method, tol.rel, n.iters, parallel,
      select.once = select.once, FALSE, ...
    )
    if (!is.character(tmp)) {
      nvar <- unique(all.vars(tmp$res@model$pred.formula))
    }

    if (!is.character(tmp) && length(nvar) > 0) {
      out$model <- tmp$res
      par <- gamBiCopPredict(out$model, target = "par")$par
      fam <- out$model@family
      par2 <- out$model@par2
    } else {
      tmp <- BiCopSelect(u1, u2, fams,
        selectioncrit = familycrit,
        indeptest = TRUE, level = level, rotations = FALSE
      )
      par <- rep(tmp$par, length(u1))
      out$model$par <- tmp$par
      out$model$family <- fam <- tmp$family
      out$model$par2 <- par2 <- tmp$par2
    }
    fam <- famTrans(fam, inv = TRUE, par = cor(u1, u2), familyset = familyset)
    if (rotation == TRUE) {
      if (fam %in% c(301, 303, 401, 403)) {
        fam <- fam + 1
      } else if (fam %in% c(302, 304, 402, 404)) {
        fam <- fam - 1
      }
    }
    if (isS4(out$model)) {
      attr(out$model, "family") <- fam
    } else {
      out$model$family <- fam
    }
  }

  ## store pseudo-observations for estimation in next tree
  fams <- vapply(1:length(par), function(j)
    famTrans(fam, inv = FALSE, par = par[j]), numeric(1))
  out$CondOn1 <- BiCopHfunc(u1, u2, fams, par, par2, check.pars = FALSE)$hfunc2
  out$CondOn2 <- BiCopHfunc(u1, u2, fams, par, par2, check.pars = FALSE)$hfunc1

  return(out)
}

as.GVC <- function(GVC, covariates) {
  n <- length(GVC$Tree) + 1
  con <- list()
  nam <- V(GVC$Tree[[1]])$name

  conditionedSets <- NULL
  correspondingModel <- NULL

  if (is.list(E(GVC$Tree[[n - 1]])$conditionedSet)) {
    conditionedSets[[n - 1]][[1]] <- (E(GVC$Tree[[n - 1]])$conditionedSet[[1]])
  }
  else {
    conditionedSets[[n - 1]][[1]] <- (E(GVC$Tree[[n - 1]])$conditionedSet)
  }

  for (k in 1:(n - 2)) {
    conditionedSets[[k]] <- E(GVC$Tree[[k]])$conditionedSet
    correspondingModel[[k]] <- as.list(E(GVC$Tree[[k]])$model)
  }

  correspondingModel[[n - 1]] <- E(GVC$Tree[[n - 1]])$model

  model.count <- get.modelCount(n)
  model <- vector("list", n * (n - 1) / 2)
  M <- matrix(NA, n, n)

  for (k in 1:(n - 1)) {
    w <- conditionedSets[[n - k]][[1]][1]

    M[k, k] <- w
    M[(k + 1), k] <- conditionedSets[[n - k]][[1]][2]

    model[[model.count[(k + 1), k]]] <- correspondingModel[[n - k]][[1]]
    if (k == (n - 1)) {
      M[(k + 1), (k + 1)] <- conditionedSets[[n - k]][[1]][2]
    } else {
      for (i in (k + 2):n) {
        for (j in 1:length(conditionedSets[[n - i + 1]])) {
          cs <- conditionedSets[[n - i + 1]][[j]]
          tmp <- correspondingModel[[n - i + 1]][[j]]
          if (cs[1] == w) {
            M[i, k] <- cs[2]
            model[[model.count[i, k]]] <- tmp
            break
          } else if (cs[2] == w) {
            M[i, k] <- cs[1]
            if (isS4(tmp)) {
              fam <- attr(tmp, "family")
            } else {
              fam <- tmp$family
            }
            if (fam %in% c(301, 303, 401, 403)) {
              fam <- fam + 1
            } else if (fam %in% c(302, 304, 402, 404)) {
              fam <- fam - 1
            }

            if (isS4(tmp)) {
              attr(tmp, "family") <- fam
            } else {
              tmp$family <- fam
            }
            model[[model.count[i, k]]] <- tmp
            break
          }
        }

        conditionedSets[[n - i + 1]][[j]] <- NULL
        correspondingModel[[n - i + 1]][[j]] <- NULL
      }
    }
  }
  M[is.na(M)] <- 0
  return(gamVine(M, model, nam, covariates))
}

valid.gamVineStructureSelect <- function(data, lin.covs, smooth.covs,
                                         simplified, type,
                                         familyset, rotations, familycrit,
                                         treecrit, level, trunclevel,
                                         tau, method, tol.rel, n.iters,
                                         parallel, verbose, select.once) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    return("data has to be either a matrix or a data frame")
  }
  tmp <- valid.gamBiCopSelect(
    data[, 1:2], lin.covs, smooth.covs, rotations,
    familyset, familycrit, level, 2, tau, method,
    tol.rel, n.iters, parallel, verbose, select.once
  )
  if (tmp != TRUE) {
    return(tmp)
  }

  n <- dim(data)[1]
  d <- dim(data)[2]

  if (d < 2) {
    return("Number of dimensions has to be at least 2.")
  }
  if (n < 2) {
    return("Number of observations has to be at least 2.")
  }
  if (any(data[, 1:d] > 1) || any(data[, 1:d] < 0)) {
    return("Data has be in the interval [0,1].")
  }

  if (!valid.logical(simplified)) {
    return(msg.logical(var2char(simplified)))
  }

  if (is.null(lin.covs) && is.null(smooth.covs) && simplified == TRUE) {
    return("When there are no covariates, the vine can't be simplified.")
  }

  if (is.null(type) || length(type) != 1 || !is.element(type, c(0, 1))) {
    return("Vine model not implemented.")
  }

  if (length(treecrit) != 1 || (treecrit != "tau" && treecrit != "rho")) {
    return("Selection criterion for the pair selection not implemented.")
  }

  if (!is.na(trunclevel) && !valid.posint(trunclevel)) {
    return("'trunclevel' should be a positive integer or NA.")
  }

  if (!is.logical(select.once)) {
    return("'select.once' must be logical")
  }

  return(TRUE)
}
