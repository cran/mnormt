\name{plot_fxy}
\alias{plot_fxy}
\title{Plotting a function of two variables}
\description{
  Plot a real-valued function \code{f} evaluated on a grid of points 
  of the Cartesian plane, possibly with parameters specified by \code{...}. 
  The type of graphical display can be regulated by selecting the plotting
  function among a set of available options.
}
\usage{
plot_fxy(f, xlim, ylim, ..., npt=51, grf, grpar) 
}
\arguments{
  \item{f}{either a function or a character string with the name of a 
     real-valued function  whose first argument  represents the 
     coordinates of points where \code{f} is evaluated;
     see \sQuote{Details} for additional information.}
   \item{xlim}{either a vector of abscissae where the \code{f} must be evaluated, 
    or a length-two vector with the endpoints of such an interval, 
    in which case \code{npt[1]} equally spaced points will be considered. }
  \item{ylim}{either a vector of ordinates where the \code{f} must be evaluated, 
    or a length-two vector with the endpoints of such an interval, 
    in which case \code{npt[2]} equally spaced points will be considered. }   
   \item{...}{additional parameters to be supplied to \code{f}; these must be
    named as expected by the specification of \code{f}.}    
  \item{npt}{either an integer value or a two-element integer vector with
    the number of equally-spaced points, within the endpoints of \code{xlim} 
    and \code{ylim}, used to set up the grid of points 
    where \code{f} is evaluated; default value: \code{51}.
    When a single value is supplied, this is expanded into a length-2 vector.
    If \code{length(xlim)>2} and \code{length(ylim)>2}, \code{npt} is ignored.
    }     
  \item{grf}{an optional character string with the name of the function  
    which produces the graphical display, selectable among 
    \code{"contour", "filled.contour", "persp", "image" } 
    of package \code{graphics};
    if \code{grf} is unset, \code{"contour"} is used.}
  \item{grpar}{an optional character string with arguments supplied to the
    selected \code{grf} function, with items separated by \code{,} 
    as in a regular call.}    
}
\details{
Function \code{f} will be called with the first argument represented by a 
two-column matrix, where each row represents a point of the grid on the 
Cartesian plane identified by \code{xlim} and \code{ylim};
this set of coordinates is stored in matrix \code{pts} of the returned list. 
If present, arguments supplied as \code{...} are also passed  to \code{f}.
It is assumed that \code{f} accepts this type of call.

The original motivation of \code{plot_fxy} was to plot instances of bivariate 
probability density functions specified by package \code{mnormt}, 
but it can be used for plotting any function fulfilling the above requirements,
as illustrated by some of the examples below. 
}

\value{
an invisible list with the following components:
\tabular{ll}{
  \code{x}\tab a vector of coordinates on the \eqn{x} axis\cr
  \code{y}\tab a vector of coordinates on the \eqn{y} axis\cr  
  \code{pts}\tab a matrix of dimension \code{(npt[1]*npt[2],2)}
            with the coordinates of the evaluation points \eqn{(x,y)} \cr
  \code{f.values} \tab the vector of \code{f} values at the \code{pts} points.             
  }
}

\seealso{
  \code{\link[graphics]{contour}}, \code{\link[graphics]{filled.contour}},
  \code{\link[graphics]{persp}}, \code{\link[graphics]{image}}
}
\examples{
Sigma <- matrix(c(1,1,1,2), 2, 2)
mean <- c(0, -1)
xlim <- c(-3, 5) 
ylim <- c(-5, 3)
#
# multivariate normal density, contour-level plot
gp <- 'col="blue", nlevels=6, main="bivariate normal density"'
u <- plot_fxy(dmnorm, xlim, ylim, mean=mean, varcov=Sigma,  grpar=gp)
cat(str(u))
#---
# multivariate normal density, filled-contour plot
plot_fxy(dmnorm,  xlim, ylim, mean=mean, varcov=Sigma,grf="filled.contour")
#---
# multivariate normal density, perspective plot
gp <- "theta = 10, phi = 25, r = 2.5"
plot_fxy(dmnorm, xlim, ylim, mean=mean, varcov=Sigma, grf="persp", grpar=gp)
#---
# multivariate Student's "t" density;
# the xlim argument passed to function 'grf' overrides the earlier xlim;
# xlim and ylim can be placed after the arguments of 'f', if one prefers so
grp <- 'xlim=c(-1, 3)'
plot_fxy(dmt, mean=mean, S=Sigma, df=8, xlim, ylim, npt=101, 
        grf="filled.contour", grpar=grp)
#---
# multivariate truncated normal density, 'image' plot
low <- c(-3, -5)
hi <- c(1, 0)
plot_fxy(dmtruncnorm, mean=mean, varcov=Sigma, lower=low, upper=hi,
        xlim, ylim, npt=81, grf="image")
#---
# multivariate truncated normal distribution function, 'image' plot;
# hence not a density function 
low <- c(-3, -5)
hi <- c(1, 0)
v <- plot_fxy(pmtruncnorm, mean=mean, varcov=Sigma, lower=low, upper=hi, 
              xlim, ylim, npt=c(61, 81), grf="image")        
#---
# a different sort of 'f' function (lbeta), not a component of this package
funct <- function(z) lbeta(a=z[,1], b=z[,2])
plot_fxy(funct, xlim=c(0.1, 2), ylim=c(0.1, 2), npt=41, 
         grpar='main="function log-beta(a,b)", xlab="a", ylab="b"')
}
\keyword{multivariate}
\keyword{hplot}
