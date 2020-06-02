\name{mtruncnorm}
\alias{dmtruncnorm}
\alias{pmtruncnorm}
\alias{rmtruncnorm}

\title{The multivariate truncated normal distribution}

\description{
  The probability density function, the distribution function and random
  number generation for the multivariate truncated normal (Gaussian)
  distribution
}
\usage{
dmtruncnorm(x, mean, varcov, lower, upper, log = FALSE, ...)
pmtruncnorm(x, mean, varcov, lower, upper, ...)
rmtruncnorm(n, mean, varcov, lower, upper)
}
 
\arguments{
  \item{x}{either a vector of length \code{d} or a matrix with \code{d}
    columns, where \code{d=ncol(varcov)}, representing  the coordinates of 
    the point(s) where the density must be evaluated. }
  \item{mean}{a vector representing the mean value of the pre-truncation normal 
    distribution.}
  \item{varcov}{a symmetric positive-definite matrix representing the variance
    matrix of the pre-truncation normal  distribution.}
  \item{lower}{a vector representing the lower truncation values of the   
    component variables; \code{-Inf} values are allowed. 
    If missing, it is set equal to \code{rep(-Inf, d)}.}
  \item{upper}{a vector representing the upper truncation values of the 
    component variables; \code{Inf} values are allowed.
    If missing, it is set equal to \code{rep(Inf, d)}.}

  \item{log}{a logical value (default value is \code{FALSE});
    if \code{TRUE}, the logarithm of the density is computed.}

  \item{\dots}{arguments passed to \code{sadmvn}, 
    among \code{maxpts}, \code{abseps}, \code{releps}.}
    
  \item{n}{the number of (pseudo) random vectors to be generated.}
}

\details{For \code{dmtruncnorm} and \code{pmtruncnorm},  
the dimension \code{d} cannot exceed \code{20}.

Function \code{rmtruncnorm} is just a wrapper of the imported function
\code{\link[tmvnsim]{tmvnsim}}, set up so that the names and the pattern of 
the \code{rmtruncnorm} arguments are in agreement with the other functions 
in the package.}
  
\value{
 a numeric vector in case of \code{dmtruncnorm} and \code{pmtruncnorm};
 a matrix in case of \code{rmtruncnorm}, unless \code{n=1} in which case 
 it is a vector.
}

\author{Adelchi Azzalini}
 
\seealso{  \code{\link[tmvnsim]{tmvnsim}} }

\examples{
m2 <- c(0.5, -1)
V2 <- matrix(c(3, 3, 3, 6), 2, 2)
lower <- a <- c(-1, -2.5)
upper <- b <- c(2, 1)
set.seed(1) 
# generate a set of coordinates, pts, on the plane
pts <- matrix(runif(10, min=-1.5, max=1.5), nrow=5, ncol=2)
rownames(pts) <- LETTERS[1:nrow(pts)]
# compute PDF and CDF at the chosen coordinates, pts
pdf <- dmtruncnorm(pts, mean=m2, varcov=V2, lower, upper) 
cdf <- pmtruncnorm(pts, mean=m2, varcov=V2, lower, upper) 
print(cbind(pts, pdf, cdf))
#--
# generate a sample of random numbers
sample <- rmtruncnorm(300, mean=m2, varcov=V2, lower, upper) 
#
# then plot the sample along with the pertaining box and the earlier points, 
# just for graphical illustration 
eps <- 0.5
ra <- lower -eps
rb <- upper + eps
plot(x=c(ra[1], rb[1]), y=c(ra[2], rb[2]), xlab="", ylab="",type="n")
points(sample, pch=20, col=4, cex=0.4)
xbox <- c(a[1], b[1], b[1], a[1], a[1])
ybox <- c(a[2], a[2], b[2], b[2], a[2])
polygon(xbox, ybox, lty=2, border="gray60")
text(pts[,1], pts[,2], col=2, labels=rownames(pts),  cex=0.9)
}
%--- 
\keyword{distribution}
\keyword{multivariate}
\concept{multivariate truncated normal distribution}
 
