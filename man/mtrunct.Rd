\name{mtrunct}
\alias{dmtrunct}
\alias{pmtrunct}
\title{The multivariate truncated Student's \emph{t} distribution}
\description{
  The probability density function and the distribution function  of the
  multivariate truncated  Student's \emph{t}  distribution 
  }
\usage{
dmtrunct(x, mean, S, df, lower, upper, log = FALSE, ...)
pmtrunct(x, mean, S, df, lower, upper, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ either a vector of length \code{d} 
     or a matrix with \code{d} columns, where \code{d=ncol(S)}, giving
     the coordinates of the point(s) where the density must be evaluated. }

  \item{mean}{either a vector of length \code{d}, representing the location
     parameter (equal to the mean vector when \code{df>1}) of the pre-truncation
     distribution or a matrix whose rows represent different mean vectors; 
     in the matrix case, its dimensions must match those of \code{x}.}

  \item{S}{a symmetric positive-definite matrix representing the 
     scale matrix, such that \code{S*df/(df-2)} is the variance-covariance
     matrix of the pre-truncation distribution  when \code{df>2}.}

  \item{df}{degrees of freedom; it must be a positive integer} 

  \item{lower}{a vector representing the lower truncation values of the   
       component variables; \code{-Inf} values are allowed. 
       If missing, it is set equal to \code{rep(-Inf, d)}.}
  \item{upper}{a vector representing the upper truncation values of the 
       component variables; \code{Inf} values are allowed.
       If missing, it is set equal to \code{rep(Inf, d)}.}

  \item{log}{a logical value (default value is \code{FALSE});
    if \code{TRUE}, the logarithm of the density is computed.}

  \item{\dots}{arguments passed to \code{sadmvt}, 
     among \code{maxpts}, \code{absrel}, \code{releps}.}

}
\details{The dimension \code{d} cannot exceed \code{20}.}

\value{a numeric vector}

% \references{%% ~put references to the literature/web site here ~}

\author{Adelchi Azzalini}

% \note{%%  ~~further notes~~}
\seealso{ \code{\link{sadmvt}} for regulating accuracy}

\examples{
m2 <- c(0.5, -1)
V2 <- matrix(c(1.5, -1.75, -1.75, 3), 2, 2)
lower <- a <- c(-1, -2.5)
upper <- b <- c(2, 1)
set.seed(1) 
points <- matrix(runif(10, -3, 3), nrow=5, ncol=2)
pdf <- dmtrunct(points, mean=m2, S=V2, df=4, lower, upper) 
cdf <- pmtrunct(points, mean=m2, S=V2, df=4, lower, upper) 
}
%----
\keyword{distribution}
\keyword{multivariate}
\concept{multivariate truncated t distribution}
 
