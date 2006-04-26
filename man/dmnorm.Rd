\name{dmnorm}
\alias{dmnorm}
\alias{pmnorm}
\alias{rmnorm}
\alias{sadmvn}
\title{Multivariate normal distribution}
\description{
  The probability density function, the distribution function and random
  number generation for the multivariate normal (Gaussian)  probability
  distribution
}
\usage{
dmnorm(x, mean = rep(0, d), varcov, log = FALSE) 
pmnorm(x, mean = rep(0, length(x)), varcov, ...) 
rmnorm(n = 1, mean = rep(0, d), varcov) 
sadmvn(lower, upper, mean, varcov, maxpts = 2000 * d, abseps = 1e-06, releps = 0) 
}
\arguments{
  \item{x}{ for \code{dmnorm}, this is either a vector of length \code{d} or 
           a matrix with \code{d} columns, where \code{d=ncol(varcov)},
           giving the coordinates of the point(s) where the density must
           be evaluated;
          for \code{pmnorm}, only a vector of length \code{d} is allowed,
          and \code{d} cannot exceed 100}
  \item{mean}{a numeric vector representing the expected value of the 
           distribution; it must be of length \code{d}, as defined above}
  \item{varcov}{a positive definite matrix representing the 
        variance-covariance  matrix of the distribution; 
        a vector of length 1 is also allowed 
       (in this case, \code{d=1} is set)} 
  \item{log}{a logical value; if \code{TRUE}, 
        the logarithm of the density is computed
            }
  \item{...}{parameters passed to \code{sadmvn}, 
          among \code{maxpts}, \code{absrel}, \code{releps} }
  \item{n}{the number of random  numbers to be generated}
  \item{lower}{a numeric vector of lower integration limits of 
        the density function; must be of maximal length 100; 
        \code{+Inf} and \code{-Inf} entries are allowed }
  \item{upper}{ a numeric vector of upper integration limits 
          of the density function; must be of maximal length 100; 
          \code{+Inf} and \code{-Inf} entries are allowed }
  \item{maxpts}{the maximum number of function evaluations 
               (default value: \code{2000*d})}
  \item{abseps}{absolute error tolerance (default value: \code{1e-6})}
  \item{releps}{relative error tolerance (default value: \code{0})}
}
\details{
  Function \code{sadmvn} is an interface to a Fortran-77 routine with 
  the same name written  by Alan Genz, and available from his web page; 
  this makes uses of some auxiliary  functions whose authors are 
  documented in the Fortran code. The routine uses an adaptive 
  integration method.
}
\value{
  \code{dmnorm} returns a vector of density values (possibly log-transformed);
  \code{pmnorm} and \code{sadmvn} return a single probability with 
   attributes giving details on the achieved accuracy;
  \code{rmnorm} returns a matrix of \code{n} rows of random vectors
}
\references{
  Genz, A. (1992).
  Numerical Computation of Multivariate Normal Probabilities.
  \emph{J. Computational and Graphical Statist.}, \bold{1}, 141-149.
  
  Genz, A. (1993).  Comparison of methods for the computation of
  multivariate normal probabilities. 
  \emph{Computing Science and Statistics}, \bold{25}, 400-405.
  
  Genz, A.:  Fortran code available at
  \url{http://www.math.wsu.edu/math/faculty/genz/software/mvn.f}
}
\author{
  Fortran code of \code{SADMVN} and most auxiliary functions by Alan Genz,
  some additional auxiliary functions by people referred to within his 
  program.  Porting to R and additional R code by Adelchi Azzalini}
\note{ 
The attributes \code{error} and \code{status} of the probability
returned by \code{pmnorm} and \code{sadmvn} indicate whether the function 
had a normal termination, achieving the required accuracy. If
this is not the case, re-run the function with an higher value of
\code{maxpts}
}
\seealso{  \code{\link{dnorm}}, \code{\link{dmt}}}
\examples{
x <- seq(-2,4,length=21)
y <- 2*x+10
z <- x+cos(y) 
mu <- c(1,12,2)
Sigma <- matrix(c(1,2,0,2,5,0.5,0,0.5,3), 3, 3)
f <- dmnorm(cbind(x,y,z), mu, Sigma)
p1 <- pmnorm(c(2,11,3), mu, Sigma)
p2 <- pmnorm(c(2,11,3), mu, Sigma, maxpts=10000, abseps=1e-10)
x <- rmnorm(10, mu, Sigma)
p <- sadmvn(lower=c(2,11,3), upper=rep(Inf,3), mu, Sigma) # upper tail
#
p1 <- pnorm(0, 1, 2)
p2 <- pmnorm(0, 1, 2^2)
}
\keyword{distribution}
\keyword{multivariate}


