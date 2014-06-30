\name{mnormt-package}
\alias{mnormt-package}
\docType{package}
\title{The 'mnormt' package: summary information}

\description{
This package provides functions for computing the density and the distribution 
function of multivariate normal and multivariate Student's \emph{t} variates
and for generating random vectors sampled from these distributions.  

}
\details{
Probabilities are computed via a non-Monte Carlo method. Different routines 
are used in the three cases \code{d=1, d=2, d>2}, if \code{d} denotes the 
number of dimensions. }

\section{Licence}{
This package and its documentation are usable under the terms of 
the \dQuote{GNU General Public License} version 3 or version 2,
as you prefer; a copy of them is available from 
\url{http://www.R-project.org/Licenses/}.}
  

\author{
Adelchi Azzalini (R code and package creation) and Alan Genz (Fortran code, 
see references below; this includes routines of other authors) 
}

\references{
 Genz, A. (1992).
  Numerical Computation of Multivariate Normal Probabilities.
  \emph{J. Computational and Graphical Statist.}, \bold{1}, 141-149.
  
  Genz, A. (1993).  Comparison of methods for the computation of
  multivariate normal probabilities. 
  \emph{Computing Science and Statistics}, \bold{25}, 400-405.
  
  Genz, A.:  Fortran code available at
  \url{http://www.math.wsu.edu/math/faculty/genz/software/fort77/mvn.f}
}
\keyword{package}
\keyword{distribution} 
\keyword{multivariate}



