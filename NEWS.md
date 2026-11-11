gamCopula (development version)
----------------------------------------------------------------

BUG FIXES

* `gamVineCopSelect` evaluates each pair-copula as it was fitted: with
  `rotations = TRUE` (the default), a double Clayton or Gumbel family was
  stored with its 90 and 270 degrees rotations exchanged. `rotations` only
  includes all rotations of the families in `familyset`, as documented.
* `gamVineStructureSelect` fits each tree on the h-functions of the
  pair-copulas fitted in the previous one, not on those of their rotations
  exchanged.
* `gamVineStructureSelect` no longer fails with igraph >= 2 ("Length of new
  attribute value must be 1 or ...") when the first candidate edge of a tree
  breaks the proximity condition.
* `gamVineCopSelect` and `gamVineStructureSelect` stop with an error when a
  covariate, or a conditioning variable in a non-simplified vine, is named
  u1, u2, w, z or xt, shares its name with another, or does not have a
  syntactically valid name; a conditioning variable named u1 or u2 was read
  as the pair-copula's argument.
* `gamBiCopFit` no longer takes a column whose name starts with "xt" for
  mgcv's `xt`.

gamCopula 0.0.8 (April 2, 2025)
----------------------------------------------------------------

NEW FEATURES

* Allow for the Frank copula

BUG FIXES

* Update documentation to fix CRAN warnings.

gamCopula 0.0.7 (February 5, 2020)
----------------------------------------------------------------

BUG FIXES

* Update documentation to fix CRAN warnings.

gamCopula 0.0.6 (December 4, 2019)
----------------------------------------------------------------

BUG FIXES

* Switch from tests like `class(.) == *` into `is(., *)`.
  
gamCopula 0.0.5 (April 17, 2019)
----------------------------------------------------------------

NEW FEATURES

* `gamVinePDF` function added to compute the pdf for a `gamVine` object.

BUG FIXES

* Fixed internal bug in `gamVineStructureSelect`.
  
gamCopula 0.0.4 (August 17, 2017)
----------------------------------------------------------------
  
BUG FIXES

* Fixed bug in `gamBiCopSelect` introduced by `select.once` when either
  `lin.cov = NULL` or `smooth.cov = NULL`.

* Fixed internal bug in `gamBiCopFit`.
  
gamCopula 0.0.3 (August 14, 2017)
----------------------------------------------------------------

DEPENDS

* Now depends explicitly on `R (>= 3.1.0)`. So far, this dependence was
    implicit through our dependence on the copula package.

NEW FEATURES

* Online API documentation on <https://tvatter.github.io/gamCopula/>.
  
* Option to select GAM structure in `gamBiCopSelect` and
  `gamVineStructureSelect` only once with the `select.once` option (now the default).

* Use `unique(familyset)` in `gamBiCopSelect` to ensure that each family is
  estimated only once.
  
BUG FIXES

* Fixed bug in `summary.gamVine` for unconditional copulas.

* Fixed bug in `gamBiCopPredict` with `alpha != 0` and `newdata = NULL`.
