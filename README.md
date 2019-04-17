gamCopula
=========

[![Linux and OSX build status](https://travis-ci.org/tvatter/gamCopula.svg?branch=master)](https://travis-ci.org/tvatter/gamCopula)
[![Windows build status](http://ci.appveyor.com/api/projects/status/github/tvatter/gamCopula?svg=true)](https://ci.appveyor.com/project/tvatter/gamCopula)
[![CRAN version](http://www.r-pkg.org/badges/version/gamCopula)](https://cran.r-project.org/package=gamCopula)
[![CRAN downloads](http://cranlogs.r-pkg.org/badges/gamCopula)](https://cran.r-project.org/package=gamCopula)

This R package implements the generalized additive modeling framework for bivariate
copulas introduced by Vatter and Chavez-Demoulin (2015) and its extension to 
Pair-Copula Constructions (or Vine Copulas) by Vatter and Nagler (2017). 
It includes functions for parameter estimation, model selection, simulation, 
and visualization. The package is still under development. Please see the [API documentation](https://tvatter.github.io/gamCopula/)
for a detailed description of all functions.

Table of contents
-----------------

- [How to install](#how-to-install)
- [Package overview](#package-overview)
	- [Bivariate copula modeling: the gamBiCop-class](#bivariate-copula-modeling-the-gambicop-class)
	- [Vine copula modeling: the gamVine-class](#vine-copula-modeling-the-gamvine-class)
	- [Bivariate copula families](#bivariate-copula-families)
- [References](#references)

------------------------------------------------------------------------

How to install
--------------


You can install:

-   the stable release on CRAN:

    ``` r
    install.packages("gamCopula")
    ```

-   the latest development version:

    ``` r
    devtools::install_github("tvatter/gamCopula")
    ```

------------------------------------------------------------------------

Package overview
----------------

Below, we list most functions and features you should know about. As usual in 
copula models, data are assumed to be serially independent and lie in the unit
hypercube. 

### Bivariate copula modeling: the gamBiCop-class

  * `gamBiCop`: Creates a GAM bivariate copula by specifying a family and model,
  namely a `gamObject` as return by the `gam` function from the `mgcv` package.
  Returns an object of class `gamBiCop`. The class has the following methods:
     
     * `print`, `summary`: a brief or comprehensive overview of the bivariate
        copula, respectively. 
            
     * `plot`: plot method based on `plot.gam` from `mgcv`.
     
     * `logLik`, `AIC`, `BIC`, `nobs`: usual fit statistics.
     
     * `EDF`: Equivalent degrees of freedom for the components of the model.
        
  * `gamBiCopSimulate`: Simulates from a bivariate GAM copula.

  * `gamBiCopFit`: Estimates parameters of a bivariate copula with a prespecified
    family. Returns an object of class `gamBiCop`.
     
  * `gamBiCopSelect`: Estimates the parameters of a bivariate copula for a set 
    of families and selects the best fitting model (using either AIC or BIC). 
    Returns an object of class `gamBiCop`.
    
  * `gamBiCopPredict`, `gamBiCopPDF`, `gamBiCopCDF`: Predict and PDF/CDF methods 
    for the GAM copula model.

### Vine copula modeling: the gamVine-class

  * `gamVine`: Creates a GAM vine copula model by specifying a tree structure and 
  list of `gamBicop` objects corresponding to each edge. Returns an object of 
  class `gamVine`. The class has the following methods:
    
    * `print`, `summary`: a brief or comprehensive overview of the bivariate
      copula, respectively. 
        
    * `plot`: plots based on `plot.gamBiCop`.

  * `gamVineSimulate`: Simulates from a GAM vine copula model.
      
  * `gamVineSeqFit`: Estimates the parameters of a GAM vine copula model with 
    prespecified structure and families.
    
  * `gamVineCopSelect`: Estimates the parameters and selects the best family for a
    GAM vine copula model with prespecified structure matrix.
    
  * `gamVineStructureSelect`: Fits a GAM vine copula model assuming no prior knowledge.
    It selects the R-vine structure using Dissmann et al. (2013)'s 
    method, estimates parameters for various families, and selects the best 
    family for each pair.

  * `gamVinePDF`: Computes the PDF for a given `gamVine` object.

  * `RVM2GVC`: converts an `RVineMatrix` object from the `VineCopula` package 
    into a `gamVine`

### Bivariate copula families

In this package several bivariate copula families are included for bivariate 
and multivariate analysis using vine copulas. It provides 
functionality of elliptical (Gaussian and Student-t) as well as Archimedean 
(Clayton, Gumbel, Frank) copulas to cover a large
range of dependence patterns. For the Clayton and Gumbel copula families,
rotated versions are included to cover negative dependence as well.

A copula family: 1 Gaussian, 2 Student t, 5 Frank, 301 Double Clayton type I (standard and rotated 90 degrees), 302 Double Clayton type II (standard and rotated 270 degrees), 303 Double Clayton type III (survival and rotated 90 degrees), 304 Double Clayton type IV (survival and rotated 270 degrees), 401 Double Gumbel type I (standard and rotated 90 degrees), 402 Double Gumbel type II (standard and rotated 270 degrees), 403 Double Gumbel type III (survival and rotated 90 degrees), 404 Double Gumbel type IV (survival and rotated 270 degrees).


The following table shows the parameter ranges of bivariate copula families with 
parameters `par` and `par2` and internal coding `family`:

| Copula family                                     | `family` | `par`         | `par2`       |
|:--------------------------------------------------|:---------|:--------------|:-------------|
| Gaussian                                          | `1`      | `(-1, 1)`     | -            |
| Student t                                         | `2`      | `(-1, 1)`     | `(2,Inf)`    |
| Double Clayton type I (standard and 90 degrees)   | `301`    | `(-Inf, Inf)` | -            |
| Double Clayton type II (standard and 270 degrees) | `302`    | `(-Inf, Inf)` | -            |
| Double Clayton type III (survival and 90 degrees) | `303`    | `(-Inf, Inf)` | -            |
| Double Clayton type IV (survival and 270 degrees) | `304`    | `(-Inf, Inf)` | -            |
| Double Gumbel type I (standard and 90 degrees)    | `401`    | `(-Inf, Inf)` | -            |
| Double Gumbel type II (standard and 270 degrees)  | `402`    | `(-Inf, Inf)` | -            |
| Double Gumbel type III (survival and 90 degrees)  | `403`    | `(-Inf, Inf)` | -            |
| Double Gumbel type IV (survival and 270 degrees)  | `404`    | `(-Inf, Inf)` | -            |
| Frank                                             | `5`      | `R \ {0}`     | -            |

------------------------------------------------------------------------

References
----------
Vatter, T., Nagler, T. (2017)  
*Generalized Additive Models for Pair-Copula Constructions*.  
Preprint available at [arXiv:1608.01593](https://arxiv.org/abs/1608.01593).

Vatter, T.,  Chavez-Demoulin, V. (2015).  
*Generalized additive models for conditional dependence structures*.  
Journal of Multivariate Analysis, 141: 147-167, http://dx.doi.org/10.1016/j.jmva.2015.07.003.
