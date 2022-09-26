\name{mt}
\alias{dmt}
\alias{pmt}
\alias{rmt}
\alias{sadmvt}
\alias{biv.nt.prob}
\alias{ptriv.nt}
\title{The multivariate Student's \emph{t} distribution}
\description{
  The probability density function, the distribution function and random number 
  generation for a \code{d}-dimensional Student's \emph{t} random variable.
}
\usage{
dmt(x, mean = rep(0, d), S, df=Inf, log = FALSE) 
pmt(x, mean = rep(0, d), S, df=Inf, ...) 
rmt(n = 1, mean = rep(0, d), S, df=Inf, sqrt=NULL) 
sadmvt(df, lower, upper, mean, S, maxpts = 2000*d, abseps = 1e-06, releps = 0) 
biv.nt.prob(df, lower, upper, mean, S)
ptriv.nt(df, x, mean, S)
}
\arguments{
  \item{x}{ either a vector of length \code{d} or (for \code{dmt} and \code{pmt}) 
     a matrix with \code{d} columns representing the coordinates of the 
     point(s) where the density must be evaluated; see also \sQuote{Details}.}
  \item{mean}{either a vector of length \code{d}, representing the location
     parameter (equal to the mean vector when \code{df>1}),  
     or (for \code{dmt} and \code{pmt}) a matrix 
     whose rows represent different mean vectors; 
     in the matrix case, its dimensions must match those of \code{x}.}
  \item{S}{a symmetric positive definite matrix with dimensions \code{(d,d)} 
     representing the  scale matrix of the distribution, 
     such that \code{S*df/(df-2)} is  the variance-covariance matrix 
      when \code{df>2};  a vector of
     length \code{1} is also allowed (in this case, \code{d=1} is set).}
  \item{df}{the degrees of freedom.  
    For \code{rmt}, it must be a positive real value or \code{Inf}. 
    For all other functions, it must be a positive integer or \code{Inf}.
    A value \code{df=Inf} is translated to a call to a suitable function
    for the the multivariate normal distribution. 
    See \sQuote{Details} for its effect for the evaluation of distribution  
    functions and other probabilities.} 
  \item{log}{a logical value(default value is \code{FALSE}); if \code{TRUE}, 
     the logarithm of the density is computed.}    
  \item{sqrt}{if not \code{NULL} (default value is \code{NULL}), 
     a square root of the intended scale matrix \code{S}; 
     see \sQuote{Details} for a full description.}   
  \item{...}{arguments passed to \code{sadmvt}, 
     among \code{maxpts}, \code{absrel}, \code{releps}.}
  \item{n}{the number of random  vectors to be generated}
  \item{lower}{a numeric vector of lower integration limits of 
     the density function; must be of maximal length \code{20}; 
     \code{+Inf} and \code{-Inf} entries are allowed.}
  \item{upper}{ a numeric vector of upper integration limits 
     of the density function; must be of maximal length \code{20}; 
     \code{+Inf} and \code{-Inf} entries are allowed }
  \item{maxpts}{the maximum number of function evaluations 
               (default value: \code{2000*d})}
  \item{abseps}{absolute error tolerance (default value: \code{1e-6}).}
  \item{releps}{relative error tolerance (default value: \code{0}).}
}


\details{
  The dimension \code{d} cannot exceed \code{20} for \code{pmt} and 
  \code{sadmvt}. If this threshold is exceeded, \code{NA} is returned.

  The functions \code{sadmvt}, \code{ptriv.mt} and \code{biv.nt.prob} are
  interfaces to Fortran 77 routines by Alan Genz, available from his web page; 
  they makes use of some auxiliary functions whose authors are indicated
  in the Fortran code itself. 
  The routine \code{sadmvt} uses an adaptive  integration method. 
  If \code{df=3}, a call to \code{pmt} activates a call to \code{ptriv.nt} 
  which  is specific for the  trivariate case, and uses Genz's  Fortran
  code \code{tvpack.f};  see Genz (2004) for the  background methodology.
  A similar fact takes place when \code{df=2} with function \code{biv.nt.prob};
  note however that the underlying Fortran code is taken from 
  \code{mvtdstpack.f}, not from \code{tvpack.f}.
  If \code{pmt} is called  with \code{d>3}, this is converted into
  a suitable call to \code{sadmvt}.
  
  If \code{sqrt=NULL} (default value), the working of \code{rmt} involves 
  computation of a square root of \code{S} via the Cholesky decomposition.
  If a non-\code{NULL} value of \code{sqrt} is supplied, it is assumed that
  it represents a square root of the scale matrix,  
  otherwise represented by \code{S}, whose value is ignored in this case.  
  This mechanism is intended primarily for use in a sequence of calls to
  \code{rmt}, all sampling from a distribution with fixed scale matrix;
  a suitable matrix \code{sqrt} can then be computed only once beforehand, 
  avoiding that the same operation is repeated multiple times along the 
  sequence of calls. For examples of use of this argument, see those in the 
  documentation of \code{\link{rmnorm}}.  
  Another use of \code{sqrt} is to supply a different form of square root 
  of the scale matrix, in place of the Cholesky factor.
  
  For efficiency reasons, \code{rmt} does not perform checks on the supplied  
  arguments.  
}

\value{
  \code{dmt} returns a vector of density values (possibly log-transformed);
  \code{pmt} and \code{sadmvt} return a single probability with 
   attributes giving details on the achieved accuracy,  provided \code{x} 
   of \code{pmnorm} is a vector;
   \code{rmt} returns a matrix of \code{n} rows of random vectors,
   or a vector in case \code{n=1} or \code{d=1}.
}

\references{
  Genz, A.:  Fortran 77 code in files \code{mvt.f}, \code{mvtdstpack.f}  
  and code{tvpack}, downloaded in 2005 and again in 2007 from his webpage, 
  whose URL as of 2020-06-01 is 
  \url{https://www.math.wsu.edu/faculty/genz/software/software.html}
  
  Genz, A. (2004). 
  Numerical computation of rectangular bivariate and trivariate normal
  and \emph{t} probabilities.
  \emph{Statistics and Computing} 14, 251-260.
  
  Dunnett, C.W. and Sobel, M. (1954).
  A bivariate generalization of Student's \emph{t}-distribution with tables  
  for certain special cases. \emph{Biometrika} 41, 153--169.
}
\author{
  \acronym{FORTRAN 77} code of \code{SADMVT}, \code{MVTDSTPACK}, \code{TVPACK}
  and many  auxiliary functions by Alan Genz;
  some additional auxiliary functions by people referred to within his 
  programs; interface to \R and additional \R code (for \code{dmt}, \code{rmt}
  etc.) by Adelchi Azzalini.}

\note{ 
The attributes \code{error} and \code{status} of the probability returned 
by \code{sadmvt} and by \code{pmt} (the latter only if \code{x} is a vector 
and \code{d>2}) indicate whether the function 
had a normal termination, achieving the required accuracy. 
If this is not the case, re-run the function with a higher value of
\code{maxpts}.
}
\seealso{\code{\link[stats:TDist]{dt}}, 
       \code{\link{rmnorm}} for use of argument \code{sqrt},
       \code{\link{plot_fxy}} for plotting examples}
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
\concept{multivariate t distribution}
