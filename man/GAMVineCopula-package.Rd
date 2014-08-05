\name{gamVineCopula-package}
\alias{gamVineCopula-package}
\alias{gamVineCopula}
\alias{gamVineCopula.package}
\docType{package}
\title{
Generalized additive models for bivariate dependence structures and vine copulas
}
\description{
This package implements various inference and simulation tools to apply generalized additive models to bivariate dependence structures and vine copulas.

More references in Vatter and Chavez-Demoulin (2014) and Vatter (2014).
}
\details{
\tabular{ll}{
Package: \tab gamVineCopula\cr
Type: \tab Package\cr
Version: \tab 0.0-1\cr
Date: \tab 2014-08-05\cr
License: \tab GPL-3\cr
}
}
\author{
Thibault Vatter

Maintainer: Thibault Vatter <thibault.vatter@gmail.com>
}
\references{
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009)
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.

Brechmann, E. C., C. Czado, and K. Aas (2012)
Truncated regular vines in high dimensions with applications to financial data.
Canadian Journal of Statistics 40 (1), 68-85.

Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013)
Selecting and estimating regular vine copulae and application to financial returns.
Computational Statistics & Data Analysis, 59 (1), 52-69.

Joe, H. (1997)
Multivariate Models and Dependence Concepts.
London: Chapman and Hall.

Schepsmeier, U. and J. Stoeber (2012).
Derivatives and Fisher information of bivariate copulas.
\url{http://link.springer.com/article/10.1007/s00362-013-0498-x}.

Vatter, T. and V. Chavez-Demoulin (2014)
Generalized additive models for conditional dependence structures.
\url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2265534}

Vatter, T. (2014)
Generalized additive models for non-simplified pair-copula constructions.

Wood, S.N. (2004)
Stable and efficient multiple smoothing parameter estimation for generalized additive models. 
J. Amer. Statist. Ass. 99:673-686. 

Wood, S.N. (2006) 
eneralized Additive Models: an introduction with R.
CRC.

}
\keyword{ package }
\seealso{
The present package is heavly relying on the \code{\link[mgcv:mgcv-package]{mgcv}} and \code{\link{VineCopula}} packages, as it basically extends and mix both of them.
}
\examples{
print("quand je l'aurai fait, ca sera vraiment bien ;)")
}
