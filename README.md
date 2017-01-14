gamCopula
=========

> Generalized additive models for bivariate copulas and pair-copula constructions

[![Linux and OSX build status](https://travis-ci.org/tvatter/gamCopula.svg?branch=master)](https://travis-ci.org/tvatter/gamCopula)
[![Windows build status](http://ci.appveyor.com/api/projects/status/github/tvatter/gamCopula?svg=true)](https://ci.appveyor.com/project/tvatter/gamCopula)


This R package implements the generalized additive modeling framework for copulas introduced by Vatter and  Chavez-Demoulin (2015).
An extension to pair-copula constructions (also called vine copulas) is also provied. The package is still under development.

You can install the latest development version as follows:

``` r
devtools::install_github("tvatter/gamCopula")
```


References
----------
Vatter, T., Nagler, T. (2016)  
*Generalized Additive Models for Pair-Copula Constructions*.  
Preprint available at [arXiv:1608.01593](https://arxiv.org/abs/1608.01593).

Vatter, T.,  Chavez-Demoulin, V. (2015).  
*Generalized additive models for conditional dependence structures*.  
Journal of Multivariate Analysis, 141: 147-167, http://dx.doi.org/10.1016/j.jmva.2015.07.003.

Schepsmeier, U., Stoeber, J., Brechmann, E. C., Graeler, B., Nagler, T., Erhardt, T. (2016).   
*VineCopula: Statistical Inference of Vine Copulas*.  
R package version 2.0.7.  https://github.com/tnagler/VineCopula.
