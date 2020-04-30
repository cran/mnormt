\name{sample_Mardia_measures}
\alias{sample_Mardia_measures}
\concept{Mardia measures of multivariate skewness and kurtosis}
\title{
The Mardia measures of multivariate skewness and kurtosis for a given sample
}
\description{
Given a multivariate sample, the Mardia measures of skewness and kurtosis are
computed, along with their \emph{p}-values for testing normality
}
\usage{
sample_Mardia_measures(data, correct = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data matrix}

  \item{correct}{(logical) if \code{correct=TRUE}, the \sQuote{corrected} 
   sample variance matrix is used, otherwise the \sQuote{uncorrected} version
   is used (default)}
}
\details{
For a given a data matrix, the multivariate  measures of skewness and 
kurtosis  introduced by Mardia (1970, 1974) are computed, 
along with some associated quantities. 
We follow the notation of the 1974 paper.

If \code{n} denotes the number of complete cases, the condition \code{n>3}
is required for numerical computation. Clearly, a much larger \code{n} is
required for meaningful statistical work.

The  sample variance matrix \eqn{S} appearing  in (2.2) and (2.4) 
is computed here (in the dafault setting) with the \eqn{n} denominator, 
at variance from the commonly employed \code{n-1} denominator. 
With this definition of \eqn{S}, one can reproduce the numerical outcomes 
of the example on p.127 of Mardia (1974).

The approximate observed significance levels for testing normality, 
\code{p.b1} and \code{p.b2}, are computed using expressions (5.5) and
(5.6) in Section 5 of Mardia (1974).
For  \code{p.b2}, the condition \code{(n-d-1)>0} is required, where
\code{d} denotes the number of variables.
}
\value{
A named vector with the following components:
 \item{b1}{the measure of asymmetry as given in (2.2)}
 \item{b2}{the measure of kurtosis as given in (2.4)}
 \item{g1}{the measure of asymmetry as given in (2.10)}
 \item{g2}{the measure of kurtosis as given in (2.11)}
 \item{p.b1}{observed significance level of \code{b1}}
 \item{p.b2}{observed significance level of \code{b2}}
 \item{n}{The number of complete cases in the input data matrix}
 where the quoted formulae are those of Mardia (1974).
}
\references{

Mardia, K. V. (1970). 
  Measures of multivariate skewness and kurtosis with applications
   \emph{Biometrika}, 57, 519-530.

Mardia, K. V. (1974).
 Applications of some measures of multivariate skewness and kurtosis 
 in testing normality and robustness studies.
 \emph{Sankhya ser.B}, 36, 115-128.
}
\author{
Adelchi Azzalini
}

%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}

\examples{
x <- rmnorm(100, mean=1:3, varcov=toeplitz(1/(1:3)))
sample_Mardia_measures(x)
}
\keyword{distribution}
\keyword{multivariate}
 
