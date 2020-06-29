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
The Matlab package \code{ftnorm} has been downloaded from
\url{http://www-2.rotman.utoronto.ca/~kan/research.htm}, 
on 2020-04-23.

The function returns an array, \code{M} say, whose entries represent integrals 
of type  \eqn{\int_a^b x^\kappa f(x) dx}, where \eqn{f(x)} denotes the 
\eqn{d}-dimensional normal density function. 
Typically, interest is in the transformed array \code{M/M[1]} whose entries 
represent moments of the truncated distribution.
  
The algorithm is based on a recursion starting from the integral of the normal 
distribution over the specified hyper-rectangle. This integral is evaluated 
by \code{sadmvn}, whose tuning parameters  \code{maxpts, abseps, releps} 
can be regulated via the \code{...} argument.       
}

\value{
  In the multivariate case, for an input vector \code{kappa=c(k1,\dots, kd)}, 
  the functions returns  an array of dimension \code{c((k1+1),...,(kd+1))} 
  whose entries represent integrals described in section \sQuote{Details}.
  In other words, the array element \code{M[i+1, j+1, k+1,...]} contains the 
  \emph{unnormalized}  cross moment of order \code{(i, j, k,...)}; 
  this must be divided by \code{M[1]} to obtain the regular cross moment.
  
  In the univariate case, a vector is returned with similar meaning.
  }

\section{Warning}{
  Although the underlying algorithm is exact in principle, the actual 
  computation hinges crucially on the initial integration of the multivariate 
  normal density over the truncation hyper-cube.
  This integration may result in numerical inaccuracies, whose amount depends 
  on the supplied arguments. Moreover, the recursion employed by the algorithm
  propagates the initial error to other terms.  
  
  When problematic cases have been processed by the original Matlab function,
  the same issues have occurred, up to minor variations.
  
  Instances of such errors may be detected when the array \code{M/M[1]} is 
  passed to \code{\link{mom2cum}}, but there is no guarantee that all such 
  problems are detected.
}
  
\references{
Kan, Raymond and Robotti, Cesare (2017).
On moments of folded and truncated multivariate normal distributions.
\emph{Journal of Computational and Graphical Statistics}, 26, 930-934, 
DOI: 10.1080/10618600.2017.1322092

Leppard, P. and Tallis, G. M. (1989).
Algorithm AS249:  Evaluation of the mean and covariance of the 
  truncated multinormal distribution
\emph{Applied Statistics} 38, 543-553)
}
\author{
Original Matlab code by Raymond Kan and Cesare Robotti,
porting to R by Adelchi Azzalini.
}
\note{This function is not intended for direct call by a user, at least in
commonly encountered situations.
Function \code{\link{mom.mtruncnorm}} represents a more user-friendly tool.}

\seealso{\code{\link{mom.mtruncnorm}} for a more user-friendly function,
  \code{\link{mom2cum}} for transformation to cumulants,
 \code{\link{sadmvn}} for regulating accuracy if \code{d>2}}

\examples{
mu <- c(1, -0.5, 0)
Sigma <- toeplitz(1/(1:3))
low <- c(-Inf, -3, -4)
hi <- c(1.5, Inf, 2)
M <- recintab(c(2,3,1), low, hi, mu, Sigma)
M/M[1]
# cross-moments up to order 2 for X1, up to the 3 for X2, up to 1 for X3,
# if the components of the trivariate variable are denoted (X1,X2,X3)
#--
# Example 2 of Leppard & Tallis (1989, Appl.Stat. vol.38, p.547)
truncp <- c(0, 1, 2)
U <- c(0, 0, 0)
V <- 0.5*(diag(3) + matrix(1, 3, 3))
M <- recintab(c(2,2,2), truncp, rep(Inf,3), U, V)
mom <- M/M[1]
EX <-  c(mom[2,1,1], mom[1,2,1], mom[1,1,2])
print(EX, digits=9)
EX2 <- matrix(c(
          mom[3,1,1], mom[2,2,1], mom[2,1,2],
          mom[2,2,1], mom[1,3,1], mom[1,2,2],
          mom[2,1,2], mom[1,2,2], mom[1,1,3]), 
          3, 3, byrow=TRUE)
varX <- EX2 - outer(EX ,EX)       
print(varX, digits=9)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{distribution}
\keyword{multivariate}
\concept{moments}
\concept{truncated multivariate normal distribution}

