\name{pd.solve}
\alias{pd.solve}
\title{Inverse of a positive definite matrix}
\description{
  The inverse of a symmetric positive-definite matrix and its log-determinant
}
\usage{
 pd.solve(x, silent = FALSE, log.det=FALSE) 
}
\arguments{
  \item{x}{a symmetric positive-definite matrix.}

  \item{silent}{a logical value which indicates the action to take in case of 
    an error. If \code{silent==TRUE} and an error occurs, the function
    silently returns a \code{NULL} value; if \code{silent==FALSE} (default) an
    error generates a \code{stop} with an error message.}

  \item{log.det}{ a logical value to indicate whether the log-determinant
    of \code{x} is required (default is \code{FALSE}).}  
} 

\value{
  the inverse matrix of \code{x};  if \code{log.det=TRUE}, this inverse has an
  attribute which contains the logarithm of the determinant of \code{x}.
}

\details{The function checks that \code{x} is a symmetric positive-definite 
  matrix. If an error is detected, an action is taken which depends on the
  value of the argument \code{silent}. 
}
   
\author{Adelchi Azzalini}

\examples{
  x <- toeplitz(rev(1:4))
  x.inv <- pd.solve(x)
  print(x.inv \%*\% x)
  x.inv <- pd.solve(x, log.det=TRUE)
  logDet <- attr(x.inv, "log.det")
  print(abs(logDet - determinant(x, logarithm=TRUE)$modulus))
}
\keyword{algebra}
\keyword{array}



