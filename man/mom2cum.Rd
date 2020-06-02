\name{mom2cum}
\alias{mom2cum}
\title{
Conversion of an array of moments to cumulants  
}
\description{
Given an array of moments of a multivariate distribution, the corresponding
cumulants up to the 4th order and other connected quantities are computed,
notably the Mardia's measures of multivariate skewness and kurtosis
}
\usage{mom2cum(mom)}

\arguments{
  \item{mom}{an array whose entries are assumed to represent moments of a
             multivariate distribution; see \sQuote{Details} for additional
             information.}
}
\details{The structure of the input array \code{mom} is the one of the output
produced by function \code{\link{recintab}}. For a \code{d}-dimensional
random variable, \code{mom} is an \code{k}-fold \code{d}-dimensional array,
where \code{k} is the highest order of moments being considered;
see the documentation of \code{mom} for a more detailed description. 
However, it is not necessary that \code{mom} originates from \code{recintab};
the moments can refer to any distribution, as long as  \code{mom} 
has the appropriate structure and content.

Also, it is not necessary that all entries of \code{mom} are there;
values not required for the processing can be left as \code{NA}.
For computing cumulants of order \code{k}, say, we only need  cross moments 
whose exponents add up to \code{k} or less.

Conversion from moments to cumulants is performed by using formulae (2.7)
of McCullagh (1987). See also \eqn{\rho_{23}^2} in his (2.15) and 
\eqn{\rho_4} in (2.16) for computing the Mardia's (1970, 1974) measures of 
multivariate skewness and kurtosis.
}
\value{
A list with the following elements, provided moments of the required order are
available, up to the maximal order 4. 
\item{cum1}{the \code{d}-vector of first-order cumulants, 
   AKA the expected value or the mean value; 
   this will be there if \code{mom} contains all moments of order 1.}
\item{order2}{a list with the following components: 
   \code{m2}, the \code{(d,d)} matrix of second order moments; 
   \code{cum2}, the \code{(d,d)} matrix of second order cumulants, 
     AKA the variance-covariance matrix, the variance matrix, the covariance matrix; 
   \code{conc.matrix}, the concentration matrix, that is, the inverse of \code{cum2};
   \code{log.det.cum2}, the logarithm of the determinant of \code{cum2}.
   } % end order2
\item{order3}{a list with the following components: 
   \code{m3}, array of third order moments, having dimension \code{(d,d,d)}; 
   \code{cum3}, array of third order cumulants, having dimension \code{(d,d,d)};
   \code{m3.marginal}, vector of third order marginal moments;
   \code{centr.mom3.marginal}, vector of third order marginal central moments; 
   \code{gamma1.marginal}, vector of third order marginal standardized cumulants;
   \code{gamma1.Mardia}, the Mardia measure of multivariate skewness;
   \code{beta1.Mardia}, the Mardia measure of multivariate skewness, again.
    } % end order3
\item{order4}{a list with the following components:
   \code{m4}, array of fourth order moments, with dimension \code{(d,d,d,d)}; 
   \code{cum4}, array of fourth order cumulants, with dimension \code{(d,d,d,d)};
   \code{m4.marginal}, vector of fourth order marginal moments;
   \code{centr.mom4.marginal}, vector of fourth order marginal central moments; 
   \code{gamma2.marginal}, vector of fourth order marginal standardized cumulants;
   \code{gamma2.Mardia}, the Mardia measure of multivariate kurtosis, \eqn{\gamma_{2,d}};
   \code{beta2.Mardia}, the Mardia measure of multivariate kurtosis, \eqn{\beta_{2,d}}. 
   }% end order4
}


\references{
Mardia, K. V. (1970). 
  Measures of multivariate skewness and kurtosis with applications
   \emph{Biometrika}, 57, 519-530.

Mardia, K. V. (1974).
 Applications of some measures of multivariate skewness and kurtosis 
 in testing normality and robustness studies.
 \emph{Sankhya ser.B}, 36, 115-128.

McCullagh, Peter (1987).
\emph{Tensor Methods in Statistics}.
Chapman & Hall, London.
}

\author{Adelchi Azzalini}

\note{In the case of a multivariate truncated normal distribution, 
a user does not need to call this function; \code{\link{mom.mtruncnorm}} 
provides a more convenient interface for the same computations.
The present function needs to be called only if the array \code{mom} 
has been computed for the moments of some other distribution.}

\seealso{\code{\link{recintab}} }

\examples{
mu <- c(1, -0.5)
Sigma <- toeplitz(1/(1:2))
low <- c(-Inf, -3)
hi <- c(1.5, Inf)
mom <- recintab(c(3,3), low, hi, mu, Sigma)
cum <- mom2cum(mom)
print(cum$order3$gamma1.marginal)
print(cum$order3$gamma1.Mardia)
}
%----
\keyword{distribution}
\keyword{multivariate}
\concept{moments}
\concept{cumulants}
\concept{Mardia's measures of multivariate skewness and kurtosis}
\concept{multivariate truncated normal distribution}
