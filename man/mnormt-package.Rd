\name{mnormt-package}
\alias{mnormt-package}
\docType{package}
\title{The 'mnormt' package: summary information}

\description{
  Functions are provided for computing the density and the  distribution 
  function of d-dimensional normal and \emph{t} random variables,
  possibly truncated (on one side or two sides, with componentwise choice),
  and for generating random vectors sampled from these distributions,
  except sampling from the truncated \emph{t}. Moments of arbitrary order of 
  a truncated normal are computed, and converted to cumulants up to order 4. }

\details{ 
  Probabilities are computed via non-Monte Carlo methods; different routines 
  are used in the case \code{d=1, d=2, d=3, d>2}, if \code{d} denotes the
  dimensionality. }
  
\section{Licence}{
  This package and its documentation are usable under the terms of the 
  \dQuote{GNU General Public License} version 3 or version 2, as you prefer; 
  a copy of them is available from 
  \url{https://www.R-project.org/Licenses/}.}
  

\author{
Adelchi Azzalini (R code and package creation) and Alan Genz (Fortran code, 
see the references below; this code incorporates routines of other authors).
% Function \code{recintab} is the \R porting by Adelchi Azzalini of the 
% Matlab function with the same name  by Kan and Robotti.
}

\references{
  Genz, A. (1992).
  Numerical Computation of Multivariate Normal Probabilities.
  \emph{J. Computational and Graphical Statist.} \bold{1}, 141-149.
  
  Genz, A. (1993).  Comparison of methods for the computation of
  multivariate normal probabilities. 
  \emph{Computing Science and Statistics}, \bold{25}, 400-405.
  
  Genz, A.:  Fortran code downloaded in 2006 from the author web page, 
  located at \url{https://www.math.wsu.edu/faculty/genz/software/software.html},
  as of 2020-06-01.
   
}
\keyword{package}
\keyword{distribution} 
\keyword{multivariate}
\concept{multivariate normal distribution}
\concept{multivariate truncated normal distribution}
\concept{multivariate t distribution}
\concept{multivariate truncated t distribution}
\concept{moments}
\concept{cumulants}



