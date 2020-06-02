\name{recintab}
\alias{recintab}
\title{
Moments of arbitrary order of a (possibly) truncated multivariate normal variable
}
\description{
Produces an array with the moments up to specified orders of a (possibly)
truncated multivariate normal distribution. 
Each component variable can be truncated on one side (to the left or 
to the right) or on two sides or not truncated.
}
\usage{
recintab(kappa, a, b, mu, S, ...)
}

\arguments{
  \item{kappa}{a vector of non-negative integer values representing the 
       required order of moments for each component variable.}
  \item{a}{a vector representing the lower truncation values of the component 
        variables; \code{-Inf} values are allowed.}
  \item{b}{a vector representing the upper truncation values of the component  
     variables; \code{Inf} values are allowed.}
  \item{mu}{a vector representing the mean value of the pre-truncation normal 
       random variable.}
  \item{S}{a symmetric positive-definite matrix representing the variance
      matrix of the pre-truncation normal random variable.}
  \item{\dots}{parameters passed to \code{sadmvn}; see the \sQuote{Details}. }
}
\details{
This function is the \R translation of the Matlab function with the same name
belonging to the package \code{ftnorm}, 
which is associated to the paper of Kan and Robotti (2017).
The Matlab package \code{ftnorm} is available at
\url{http://www-2.rotman.utoronto.ca/~kan/research.htm}, 
accessed on 2020-04-23.

The algorithm is based on a recursion starting from the integral of the normal 
distribution over the specified hyper-rectangle. This integral is evaluated 
by \code{sadmvn}, whose tuning parameters  \code{maxpts, abseps, releps} 
can be regulated via the \code{...} argument. 
      
}
\value{an array of dimension \code{c((k1+1),...,(kd+1))} 
      if \code{kappa=c(k1,\dots, kd)}, such that the cross moment of order
      \code{(i, j, k,...)} is given at entry \code{[i+1, j+1, k+1,...]}.}
\references{
Kan, Raymond and Robotti, Cesare (2017).
On moments of folded and truncated multivariate normal distributions.
\emph{Journal of Computational and Graphical Statistics}, 26, 930-934, 
DOI: 10.1080/10618600.2017.1322092
}
\author{
Original Matlab code by Raymond Kan and Cesare Robotti,
porting to R by Adelchi Azzalini.
}
\note{This function is not intended for direct call by a user.
Function \code{\link{mom.mtruncnorm}} represents a more user-friendly tool.}

\seealso{\code{\link{mom.mtruncnorm}} for a more user-friendly function,
 \code{\link{sadmvn}} for regulating accuracy}

\examples{
mu <- c(1, -0.5, 0)
Sigma <- toeplitz(1/(1:3))
low <- c(-Inf, -3, -4)
hi <- c(1.5, Inf, 2)
mom <- recintab(c(2,3,1), low, hi, mu, Sigma)
# cross-moments up to order 2 for X1, up to the 3 for X2, up to 1 for X3
# if the components of the trivariate variable are denoted (X1,X2,X3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{distribution}
\keyword{multivariate}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{moments}
\concept{truncated normal distribution}

