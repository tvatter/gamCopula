gamCopula
=========

> Generalized additive models for bivariate copulas and pair-copula constructions


This R package implements the generalized additive modeling framework for copulas introduced by Vatter and  Chavez-Demoulin (2015).
An extension to pair-copula constructions (also called vine copulas) is also provied. The package is still under development.

You can install the latest development version as follows:

``` r
devtools::install_github("tvatter/gamCopula")
```
The pacakge requires the latest version of the VineCopula package (Schepsmeier et al., 2016). To install:

``` r
devtools::install_github("tnagler/VineCopula")
```


References
----------
Vatter, T.,  Chavez-Demoulin, V. (2015).  
*Generalized additive models for conditional dependence structures*.  
Journal of Multivariate Analysis, 141: 147-167, http://dx.doi.org/10.1016/j.jmva.2015.07.003.

Schepsmeier, U., Stoeber, J., Brechmann, E. C., Graeler, B., Nagler, T., Erhardt, T. (2016).   
*VineCopula: Statistical Inference of Vine Copulas*.  
R package version 2.0.0.  https://github.com/tnagler/VineCopula
