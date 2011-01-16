\name{pd.solve}
\alias{pd.solve}
\title{Inverse of a positive definite matrix}
\description{
  The inverse of a symmetrix positive definite matrix and its log-determinant
}
\usage{
 pd.solve(x, silent = FALSE) 
}
\arguments{
  \item{x}{a symmetric positive definite matrix}
  \item{silent}{a logical value which indicates the action to take in case of an
    error. If \code{silent==TRUE} an error produces a silent \code{NULL} value
    of the function;  if \code{silent==FALSE} (default) an error generates 
    a \code{stop} with an error message}
} 

\value{
  \code{pd.solve} returns the inverse matrix of \code{x}, with an attribute
  \code{"log.det"} representing the logarithm of the determinant of \code{x}
}

\details{The function checks that \code{x} is a symmetric positive definite 
  matrix. If an error is detected, an action is taken which depends on the
  value of the argument \code{silent}. 
}
   
\author{Adelchi Azzalini}

\examples{
  x <- toeplitz(rev(1:4))
  x.inv <- pd.solve(x)
  print(x.inv \%*\% x)
  logDet <- attr(x.inv, "log.det")
  print(abs(logDet - determinant(x, logarithm=TRUE)$modulus))
}
\keyword{algebra}
\keyword{array}



