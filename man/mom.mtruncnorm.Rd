\name{mom.mtruncnorm}
\alias{mom.mtruncnorm}
\title{Moments and other quantities of a (possibly) truncated multivariate 
      normal distribution
}
\description{
Moments up to the specified orders of a possibly truncated \code{d}-dimensional
normal distribution; the distribution must be non-degenerate. 
Each component variable can be truncated on one side (to the left or 
to the right) or on two sides or not truncated.
After the initial stage, cumulants up to the fourth order and other
quantities are computed, provided all moments of the required order 
had been computed in the first stage.
}
\usage{
mom.mtruncnorm(powers=4, mean, varcov, lower, upper, cum = TRUE, ...)
}
 
\arguments{
  \item{powers}{a vector of non-negative integer values representing the 
       required order of moments for each component variable, or a single
       such value, in which case this value is repeated for all \code{d}
       components.}
  \item{mean}{a vector representing the mean value of the pre-truncation normal 
       random variable.}
  \item{varcov}{a symmetric positive-definite matrix representing the variance
       matrix of the pre-truncation normal random variable.}
  \item{lower}{a vector representing the lower truncation values of the   
       component variables; \code{-Inf} values are allowed. 
       If missing, it is set equal to \code{rep(-Inf, d)}.}
  \item{upper}{a vector representing the upper truncation values of the 
       component variables; \code{Inf} values are allowed.
       If missing, it is set equal to \code{rep(Inf, d)}.}
  \item{cum}{a logical value; if  code{TRUE}  (default value), cumulants 
       are other quantities are computed up to the minimum between the fourth
       order and the maximum possible order, given the available moments.}
  \item{\dots}{additional arguments passed to \code{sadmvn}; 
       see \sQuote{Details} for a description.}
  }
  
\details{
This function makes use of two workhorses, \code{recintab} and \code{mom2cum}, 
providing a user-friendly interface to these more basic tools. 
The first function computes an array of raw moments of the truncated normal; 
the second function translates them into cumulants and other quantities
such as the Mardia's measures of skewness and kurtosis, unless \code{cum=FALSE}.
See the documentation of these two underlying functions for additional 
information about the arguments and the returned quantities. 
The argument \code{...} is passed, via \code{recintab}, to \code{sadmvn} 
for regulation of its working.

Not all \code{d} component variables need to be truncated. 
In fact, the function works also with no truncated components 
(just omit \code{lower} and \code{upper}),
although for this case there exist known formulae to do the job.
}

\value{A list with the following components:
  \item{mom}{an array with raw moments as produced by \code{recintab},
     followed by normalization; see its documentation for a description.}
  \item{cum1}{the vector of first-order cumulants, 
     \acronym{AKA} the expected value or  the mean value, 
     which will be there  provided \code{cum=TRUE} and all 
     elements of \code{powers} are not less than 1.}
  \item{order2, ...}{additional lists with higher order terms up to order 4;
     these lists only exist when the available moments are of sufficiently 
     high order. See \code{\link{mom2cum}} for a more detailed description.}
}

\author{Adelchi Azzalini}
% \note{%%  ~~further notes~~}
 
\seealso{
  \code{\link{recintab}},  \code{\link{mom2cum}},  \code{\link{sadmvn}}
}
\examples{
mu <- c(1, -0.5, 0)
Sigma <- toeplitz(1/(1:3))
lower <- c(-Inf, -3, -4)
upper <- c(1.5, Inf, 2)
m <- mom.mtruncnorm(1, mu, Sigma, lower, upper)
print(m$cum1)
# 
m <- mom.mtruncnorm(3, mu, Sigma, lower, upper)
print(m$order3$gamma1.marginal)
print(m$order3$gamma1.Mardia)
#
#--
# Example 2 of Leppard & Tallis (1989, Appl.Stat. vol.38, p.547)
truncp <- c(0, 1, 2)
U <- c(0, 0, 0)
V <- 0.5*(diag(3) + matrix(1, 3, 3))
m <- mom.mtruncnorm(2, U, V, truncp)
print(m$cum1, digits=9)
print(m$order2$cum2, digits=9)
}
%----
\keyword{distribution}
\keyword{multivariate}
\concept{moments}
\concept{multivariate truncated normal distribution}

