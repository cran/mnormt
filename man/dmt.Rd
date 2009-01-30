\name{dmt}
\alias{dmt}
\alias{pmt}
\alias{rmt}
\alias{sadmvt}
\alias{biv.nt.prob}
\title{Multivariate t distribution}
\description{
  The probability density function, the distribution function and random
  number generation for the multivariate t  probability
  distribution
}
\usage{
dmt(x, mean = rep(0, d), S, df=Inf, log = FALSE) 
pmt(x, mean = rep(0, length(x)), S, df=Inf, ...) 
rmt(n = 1, mean = rep(0, d), S, df=Inf) 
sadmvt(df, lower, upper, mean, S, maxpts = 2000 * d, abseps = 1e-06, releps = 0) 
biv.nt.prob(df, lower, upper, mean, S)
}
\arguments{
  \item{x}{ for \code{dmt}, this is either a vector of length \code{d} 
    or a matrix with \code{d} columns, where \code{d=ncol(S)}, giving
    the coordinates of the point(s) where the density must be evaluated;
    for \code{pmt}, only a vector of length \code{d} is allowed,
    and \code{d} cannot exceed 20}
  \item{mean}{a numeric vector representing the location parameter
         of the distribution (equal to the expected value when \code{df>1});
           it must be of length \code{d}, as defined above}
  \item{S}{a positive definite matrix representing the 
        scale matrix of the distribution, such that \code{S*df/(df-2)} is
        the variance-covariance matrix  when \code{df>2}; 
        a vector of length 1 is also allowed 
       (in this case, \code{d=1} is set)}
  \item{df}{degrees of freedom; it must be a positive integer for
       \code{pmt}, \code{sadmvt} and \code{biv.nt.prob}, 
        otherwise a positive number. 
        If \code{df=Inf} (default value), the corresponding \code{*mnorm}
       function is called, unless \code{d=2}; in this case
       \code{biv.nt.prob} is used, If  \code{biv.nt.prob} is called with
       \code{df=Inf}, it returns the probability of a rectangle assigned by
       a bivariate normal distribution} 
  \item{log}{a logical value; if \code{TRUE}, 
        the logarithm of the density is computed
        }    
  \item{...}{parameters passed to \code{sadmvt}, 
          among \code{maxpts}, \code{absrel}, \code{releps} }
  \item{n}{the number of random  numbers to be generated}
  \item{lower}{a numeric vector of lower integration limits of 
        the density function; must be of maximal length 20; 
        \code{+Inf} and \code{-Inf} entries are allowed }
  \item{upper}{ a numeric vector of upper integration limits 
          of the density function; must be of maximal length 20; 
          \code{+Inf} and \code{-Inf} entries are allowed }
  \item{maxpts}{the maximum number of function evaluations 
               (default value: \code{2000*d})}
  \item{abseps}{absolute error tolerance (default value: \code{1e-6})}
  \item{releps}{relative error tolerance (default value: \code{0})}
}
\details{
  The functions \code{sadmvt} and \code{biv.nt.prob} are interfaces to 
  Fortran-77  routines by Alan Genz, and available from his web page; 
  they makes uses of some auxiliary  functions whose authors are 
  documented in the Fortran code. The routine \code{sadmvt} uses an adaptive 
  integration method. The routine \code{biv.nt.prob} is specific for the
  bivariate case; if \code{df<1} or \code{df=Inf}, it computes the bivariate 
  normal distribution function using a non-iterative method described in a
  reference given below.
  If \code{pmt} is called  with \code{d>2}, this is converted into
  a suitable call to \code{sadmvt}; if \code{d=2}, a call to 
  \code{biv.nt.prob} is used; if \code{d=1}, then \code{pt} is used.
}
\value{
  \code{dmt} returns a vector of density values (possibly log-transformed);
  \code{pmt} and \code{sadmvt} return a single probability with 
   attributes giving details on the achieved accuracy;
   \code{rmt} returns a matrix of \code{n} rows of random vectors
   }
\references{
  Genz, A.:  Fortran code in files \code{mvt.f} and \code{mvtdstpack.f} 
  available at \url{http://www.math.wsu.edu/math/faculty/genz/software/} 
  
  Dunnett, C.W. and Sobel, M. (1954).
  A bivariate generalization of Student's \emph{t}-distribution with tables  
  for certain special cases. \emph{Biometrika} 41, 153--169.
}
\author{
  Fortran code of \code{SADMVT} and most auxiliary functions by Alan Genz,
  some additional auxiliary functions by people referred to within his 
  program.  Porting to R and additional R code by Adelchi Azzalini}
\note{ 
The attributes \code{error} and \code{status} of the probability
returned by \code{pmt} and \code{sadmvt} indicate whether the function 
had a normal termination, achieving the required accuracy. If
this is not the case, re-run the function with an higher value of
\code{maxpts}
}
\seealso{\code{\link{dt}}, \code{\link{dmnorm}}}
\examples{
x <- seq(-2,4,length=21)
y <- 2*x+10
z <- x+cos(y) 
mu <- c(1,12,2)
Sigma <- matrix(c(1,2,0,2,5,0.5,0,0.5,3), 3, 3)
df <- 4
f  <- dmt(cbind(x,y,z), mu, Sigma,df)
p1 <- pmt(c(2,11,3), mu, Sigma, df)
p2 <- pmt(c(2,11,3), mu, Sigma, df, maxpts=10000, abseps=1e-8)
x  <- rmt(10, mu, Sigma, df)
p  <- sadmvt(df, lower=c(2,11,3), upper=rep(Inf,3), mu, Sigma) # upper tail
#
p0 <- pmt(c(2,11), mu[1:2], Sigma[1:2,1:2], df=5)
p1 <- biv.nt.prob(5, lower=rep(-Inf,2), upper=c(2, 11), mu[1:2], Sigma[1:2,1:2])
p2 <- sadmvt(5, lower=rep(-Inf,2), upper=c(2, 11), mu[1:2], Sigma[1:2,1:2]) 
c(p0, p1, p2, p0-p1, p0-p2)

}
\keyword{distribution}
\keyword{multivariate}


