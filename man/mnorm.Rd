\name{mnorm}
\alias{dmnorm}
\alias{pmnorm}
\alias{rmnorm}
\alias{sadmvn}
\title{The multivariate normal distribution}
\description{
  The probability density function, the distribution function and random
  number generation for a \code{d}-dimensional multivariate normal (Gaussian) 
  random variable.
}
\usage{
dmnorm(x, mean = rep(0, d), varcov, log = FALSE) 
pmnorm(x, mean = rep(0, d), varcov, ...) 
rmnorm(n = 1, mean = rep(0, d), varcov, sqrt=NULL) 
sadmvn(lower, upper, mean, varcov, maxpts = 2000*d, abseps = 1e-06, releps = 0)
}
\arguments{
  \item{x}{either a vector of length \code{d} or a matrix with \code{d}
    columns representing  the coordinates of the
    point(s) where the density must be evaluated; 
    see also \sQuote{Details} for restrictions on \code{d}.}
    
  \item{mean}{either a vector of length \code{d}, representing the mean value,
    or (except for \code{rmnorm}) a matrix whose rows represent different 
    mean vectors; in the matrix case, only allowed for \code{dmnorm} and 
    \code{pmnorm}, its dimensions must match those of \code{x}.}
   
  \item{varcov}{a symmetric positive-definite matrix representing the 
    variance-covariance  matrix of the distribution; 
    a vector of length 1 is also allowed (in this case, \code{d=1} is set).} 
    
  \item{sqrt}{if not \code{NULL} (default value is \code{NULL}), 
    a square root of the  intended \code{varcov} matrix; 
    see \sQuote{Details} for a full description.}

  \item{log}{a logical value (default value is \code{FALSE});
    if \code{TRUE}, the logarithm of the density is computed.}
    
  \item{...}{arguments passed to \code{sadmvn}, 
    among \code{maxpts}, \code{abseps}, \code{releps}.}
    
  \item{n}{the number of (pseudo) random  vectors to be generated.}
  
  \item{lower}{a numeric vector of lower integration limits of 
    the density function; must be of maximal length \code{20}; 
    \code{+Inf} and \code{-Inf} entries are allowed.}
    
  \item{upper}{ a numeric vector of upper integration limits 
    of the density function; must be of maximal length \code{20}; 
    \code{+Inf} and \code{-Inf} entries are allowed.}
  
  \item{maxpts}{the maximum number of function evaluations 
    (default value: \code{2000*d}).}
  \item{abseps}{absolute error tolerance (default value: \code{1e-6}).}
  \item{releps}{relative error tolerance (default value: \code{0}).}
}


\details{
  The dimension \code{d} cannot exceed \code{20} for \code{pmnorm} and 
  \code{sadmvn}. If this threshold is exceeded, \code{NA} is returned.

  The function \code{pmnorm} works by making a suitable call to
  \code{sadmvn} if \code{d>3},  or to \code{ptriv.nt} if \code{d=3},
  or to \code{biv.nt.prob} if \code{d=2},  or to \code{pnorm} if \code{d=1}.
  The \R functions \code{sadmvn}, \code{ptriv.nt} and \code{biv.nt.prob} are,
  in essence, interfaces to underlying \acronym{Fortran 77} routines by Alan 
  Genz; see the references below.  
  These routines use adaptive numerical quadrature and other non-random 
  type techniques.

  If \code{sqrt=NULL} (default value), the working of \code{rmnorm} involves 
  computation of a square root of \code{varcov} via the Cholesky decomposition.
  If a non-\code{NULL} value of \code{sqrt} is supplied, it  is  assumed 
  that it represents a matrix, \eqn{R} say, such that \eqn{R' R}
  represents the required variance-covariance matrix of the distribution;
  in this case, the argument \code{varcov} is ignored.  
  This mechanism is intended primarily for use in a sequence of calls to
  \code{rmnorm}, all sampling from a distribution with fixed variance matrix;
  a suitable matrix \code{sqrt} can then be computed only once beforehand, 
  avoiding that the same operation is repeated multiple times along the 
  sequence of calls; see the examples below.
  Another use of \code{sqrt} is to supply a different form of square root 
  of the variance-covariance matrix, in place of the Cholesky factor.
  
  For efficiency reasons, \code{rmnorm} does not perform checks on the supplied  
  arguments.
  
  If, after setting the same seed value to \code{\link[base:Random]{set.seed}},
  two calls to \code{rmnorm} are made with the same arguments except that one
  generates \code{n1} vectors and the other \code{n2} vectors, with
  \code{n1<n2}, then the \code{n1} vectors of the first call coincide with the
  initial \code{n2} vectors of the second call.  
  }

\value{
  \code{dmnorm} returns a vector of density values (possibly log-transformed);
  \code{pmnorm} returns a vector of probabilities, possibly with attributes
  on the accuracy in case \code{x}  is a vector;
   \code{sadmvn} return a single probability with 
   attributes giving details on the achieved accuracy;
  \code{rmnorm} returns a matrix of \code{n} rows of random vectors,
  or a vector in case \code{n=1} or \code{d=1}.
}

\references{
  Genz, A. (1992).
  Numerical Computation of multivariate normal probabilities.
  \emph{J. Computational and Graphical Statist.}, \bold{1}, 141-149.
  
  Genz, A. (1993).  Comparison of methods for the computation of
  multivariate normal probabilities. 
  \emph{Computing Science and Statistics}, \bold{25}, 400-405.
  
  Genz, A.:  Fortran 77 code downloaded in 2005 and again in 2007 
  from his web-page,  whose \acronym{URL} as of 2020-04-28 was
  \url{https://www.math.wsu.edu/faculty/genz/software/software.html}
  }

\author{
  \acronym{FORTRAN 77} code of \code{SADMVN}, package \code{mvtdstpack.f}, 
  package \code{tvpack}  and most auxiliary functions by Alan Genz; 
  some additional auxiliary functions by people referred to within his programs;  
  interface to \R and additional \R code (for \code{dmnormt}, 
  \code{rmnormt}, etc.)  by Adelchi Azzalini.}

\note{ 
The attributes \code{error} and \code{status} of the probability
returned by \code{pmnorm} and \code{sadmvn} indicate whether the function 
had a normal termination, achieving the required accuracy. 
If this is not the case, re-run the function with a higher value of
\code{maxpts}
}
\seealso{\code{\link[stats:Normal]{dnorm}}, \code{\link{dmt}}, 
         \code{\link{biv.nt.prob}}, \code{\link{ptriv.nt}},
         \code{\link{plot_fxy}} for plotting examples}
\examples{
x <- seq(-2, 4, length=21)
y <- cos(2*x) + 10
z <- x + sin(3*y) 
mu <- c(1,12,2)
Sigma <- matrix(c(1,2,0,2,5,0.5,0,0.5,3), 3, 3)
f <- dmnorm(cbind(x,y,z), mu, Sigma)
f0 <- dmnorm(mu, mu, Sigma)
p1 <- pmnorm(c(2,11,3), mu, Sigma)
p2 <- pmnorm(c(2,11,3), mu, Sigma, maxpts=10000, abseps=1e-10)
p <- pmnorm(cbind(x,y,z), mu, Sigma)
#
set.seed(123)
x1 <- rmnorm(5, mu, Sigma)
set.seed(123)
x2 <- rmnorm(5, mu, sqrt=chol(Sigma)) # x1=x2
eig <- eigen(Sigma, symmetric = TRUE)
R <- t(eig$vectors \%*\% diag(sqrt(eig$values)))
for(i in 1:50) x <- rmnorm(5, mu, sqrt=R)
#
p <- sadmvn(lower=c(2,11,3), upper=rep(Inf,3), mu, Sigma) # upper tail
#
p0 <- pmnorm(c(2,11), mu[1:2], Sigma[1:2,1:2])
p1 <- biv.nt.prob(0, lower=rep(-Inf,2), upper=c(2, 11), mu[1:2], Sigma[1:2,1:2])
p2 <- sadmvn(lower=rep(-Inf,2), upper=c(2, 11), mu[1:2], Sigma[1:2,1:2]) 
c(p0, p1, p2, p0-p1, p0-p2)
#
p1 <- pnorm(0, 1, 3)
p2 <- pmnorm(0, 1, 3^2)
}
\keyword{distribution}
\keyword{multivariate}
\concept{multivariate normal distribution}

