# R code of package 'mnormt', version 1.4-1 (2011-03-21). 
# Author: Adelchi Azzalini (University of Padua, Italy) 

dmnorm <- function(x, mean=rep(0,d), varcov, log=FALSE)
{
  d <- if(is.matrix(varcov)) ncol(varcov) else 1
  if(d==1) return(dnorm(x, mean, sqrt(varcov), log=log))
  x <- if (is.vector(x)) matrix(x, 1, d) else data.matrix(x)
  if (is.vector(mean)) mean <- outer(rep(1, nrow(x)), mean)
  if (is.matrix(mean) && (nrow(mean) != nrow(x) || ncol(mean) != ncol(x)))
     stop("mismatch of dimensions of 'x' and 'mean'")
  if(is.vector(mean)) mean <- outer(rep(1, nrow(x)), mean)
  X  <- t(x - mean)
  conc <- pd.solve(varcov, log.det=TRUE)
  Q  <- apply((conc %*% X)* X, 2, sum) 
  log.det <- attr(conc, "log.det")
  logPDF <- as.vector(Q + d*logb(2*pi)+log.det)/(-2)
  if(log) logPDF else exp(logPDF)
}


rmnorm <- function(n=1, mean=rep(0,d), varcov)
 {
  d <- if(is.matrix(varcov)) ncol(varcov) else 1
  z <- matrix(rnorm(n*d),n,d) %*% chol(varcov)
  y <- t(mean+t(z))
  return(y)
 }


pmnorm <- function(x, mean=rep(0,length(x)), varcov, ...)
  if(length(x) == 2) 
    biv.nt.prob(0, lower=rep(-Inf, 2), upper=x, mean, varcov)  else  
    sadmvn(lower=rep(-Inf, length(x)), upper=x, mean, varcov, ...) 
 

sadmvn <- function(lower, upper, mean, varcov,  
                   maxpts=2000*d, abseps=1e-6, releps=0)
{
  if(any(lower > upper)) stop("lower>upper integration limits")
  if(any(lower == upper)) return(0)
  d <- as.integer(if(is.matrix(varcov)) ncol(varcov) else 1)
  varcov <- matrix(varcov, d, d)
  sd  <- sqrt(diag(varcov))
  rho <- cov2cor(varcov)
  lower <- as.double((lower-mean)/sd)
  upper <- as.double((upper-mean)/sd)
  if(d == 1) return(pnorm(upper)-pnorm(lower))
  infin <- rep(2,d)
  infin <- replace(infin, (upper == Inf) & (lower > -Inf), 1)
  infin <- replace(infin, (upper < Inf) & (lower == -Inf), 0)
  infin <- replace(infin, (upper == Inf) & (lower == -Inf), -1)
  infin <- as.integer(infin)
  if(any(infin == -1)) {
    if(all(infin == -1)) return(1)
    k <- which(infin != -1)
    d <- length(k)
    lower <- lower[k]
    upper <- upper[k]
    if(d == 1) return(pnorm(upper) - pnorm(lower))
    rho <- rho[k, k]
    infin <- infin[k]
    if(d == 2) return(biv.nt.prob(0, lower, upper, rep(0,2), rho))
    }
  lower <- replace(lower, lower == -Inf, 0)
  upper <- replace(upper, upper == Inf, 0)
  correl <- as.double(rho[upper.tri(rho, diag=FALSE)])
  maxpts <- as.integer(maxpts)
  abseps <- as.double(abseps)
  releps <- as.double(releps)
  error  <- as.double(0)
  value  <- as.double(0)
  inform <- as.integer(0)
  result <- .Fortran("sadmvn", d, lower, upper, infin, correl, maxpts,
               abseps, releps, error, value, inform, PACKAGE="mnormt")
  prob <- result[[10]]
  attr(prob,"error")  <- result[[9]]
  attr(prob,"status") <- switch(1+result[[11]], 
                "normal completion", "accuracy non achieved", "oversize")
  return(prob)
}

#----

dmt <- function (x, mean=rep(0,d), S, df = Inf, log = FALSE) 
{
  if (df == Inf)  return(dmnorm(x, mean, S, log = log))
  d <- if(is.matrix(S)) ncol(S) else 1
  if (d==1) {
    y <- dt((x-mean)/sqrt(S), df=df, log=log)
    if(log) y <- y-0.5*logb(S) else y<- y/sqrt(S)
    return(y)
    }
  x <- if (is.vector(x)) matrix(x, 1, d) else data.matrix(x)
  if (is.vector(mean)) mean <- outer(rep(1, nrow(x)), mean)
  if (is.matrix(mean) && (nrow(mean) != nrow(x) || ncol(mean) != ncol(x)))
     stop("mismatch of dimensions of 'x' and 'mean'")
  if(is.vector(mean)) mean <- outer(rep(1, nrow(x)), mean)
  X  <- t(x - mean)
  S.inv <- pd.solve(S, log.det=TRUE)
  Q <- apply((S.inv %*% X) * X, 2, sum)
  logDet <- attr(S.inv, "log.det")
  logPDF <- (lgamma((df + d)/2) - 0.5 * (d * logb(pi * df) + logDet)
             - lgamma(df/2) - 0.5 * (df + d) * logb(1 + Q/df))
  if(log) logPDF else exp(logPDF)
}

rmt <- function(n=1, mean=rep(0,d), S, df=Inf)
{ 
  d <- if(is.matrix(S)) ncol(S) else 1 
  x <- if(df==Inf) 1 else rchisq(n,df)/df
  z <- rmnorm(n, rep(0,d), S)
  y <- t(mean + t(z/sqrt(x)))
  return(y)
}


pmt <- function(x, mean=rep(0,length(x)), S, df=Inf, ...){
  if(length(x) == 2) 
    biv.nt.prob(df, lower=rep(-Inf, 2), upper=x, mean, S) else  
    sadmvt(df, lower=rep(-Inf, length(x)), upper=x, mean, S, ...)  
  }
   


sadmvt <- function(df, lower, upper, mean, S, 
                   maxpts=2000*d, abseps=1e-6, releps=0)
{
  if(df == Inf) return(sadmvn(lower, upper, mean, S, maxpts, abseps, releps))
  if(any(lower > upper)) stop("lower>upper integration limits")
  if(any(lower == upper)) return(0)
  if(round(df) != df) warning("non integer df is rounded to integer") 
  df <- as.integer(round(df))
  d  <- as.integer(if(is.matrix(S)) ncol(S) else 1)
  S  <- matrix(S, d, d)
  sd  <- sqrt(diag(S))
  rho <- cov2cor(S)
  lower <- as.double((lower-mean)/sd)
  upper <- as.double((upper-mean)/sd)
  if(d == 1) return(pt(upper, df) - pt(lower, df))
  infin <- rep(2,d)
  infin <- replace(infin, (upper == Inf) & (lower > -Inf), 1)
  infin <- replace(infin, (upper < Inf) & (lower == -Inf), 0)
  infin <- replace(infin, (upper == Inf) & (lower == -Inf), -1)
  infin <- as.integer(infin)
  if(any(infin == -1)) {
    if(all(infin == -1)) return(1)
    k <- which(infin != -1)
    d <- length(k)
    lower <- lower[k]
    upper <- upper[k]
    if(d == 1) return(pt(upper, df=df) - pt(lower, df=df))
    rho <- rho[k, k]
    infin <- infin[k]
    if(d == 2) return(biv.nt.prob(df, lower, upper, rep(0,2), rho))
    }
  lower <- replace(lower, lower == -Inf, 0)
  upper <- replace(upper, upper == Inf, 0)
  correl <- rho[upper.tri(rho, diag=FALSE)]
  maxpts <- as.integer(maxpts)
  abseps <- as.double(abseps)
  releps <- as.double(releps)
  error  <- as.double(0)
  value  <- as.double(0)
  inform <- as.integer(0)
  result <- .Fortran("sadmvt", d, df, lower, upper, infin, correl, maxpts,
                   abseps, releps, error, value, inform, PACKAGE="mnormt")
  prob <- result[[11]]
  attr(prob,"error")  <- result[[10]]
  attr(prob,"status") <- switch(1+result[[12]], 
                "normal completion", "accuracy non achieved", "oversize")
  return(prob)
}


biv.nt.prob <- function(df, lower, upper, mean, S){
  if(any(dim(S) != c(2,2))) stop("dimensions mismatch")
  if(length(mean) != 2) stop("dimensions mismatch") 
  if(round(df) != df) warning("non integer df is rounded to integer") 
  nu <- if(df<Inf) as.integer(round(df)) else 0
  if(df==Inf) nu <- 0
  sd <- sqrt(diag(S))
  rho <- cov2cor(S)[1,2]
  lower <- as.double((lower-mean)/sd)
  upper <- as.double((upper-mean)/sd)
  if(any(lower > upper)) stop("lower>upper integration limits")
  if(any(lower == upper)) return(0)
  infin <- c(2,2)
  infin <- replace(infin, (upper == Inf) & (lower > -Inf), 1)
  infin <- replace(infin, (upper < Inf) & (lower == -Inf), 0)
  infin <- replace(infin, (upper == Inf) & (lower == -Inf), -1)
  infin <- as.integer(infin)
  if(any(infin == -1)) {
    if(all(infin == -1)) return(1)
    k <- which(infin != -1)
    return(pt(upper[k], df=df) - pt(lower[k], df=df))
    }
  lower <- replace(lower, lower == -Inf, 0)
  upper <- replace(upper, upper == Inf, 0)
  rho   <- as.double(rho)
  prob  <- as.double(0)
  a <- .Fortran("smvbvt", prob, nu, lower, upper, infin, rho, PACKAGE="mnormt")
  return(a[[1]])
  } 
 
pd.solve <- function(x, silent=FALSE, log.det=FALSE)
{
  if(is.null(x)) return(NULL)
  if(any(is.na(x)))
    { if(silent) return (NULL) else stop("NA's in x") } 
  if(length(x) == 1) x <- as.matrix(x)
  if(!is.matrix(x)) 
    {if(silent) return(NULL) else stop("x is not a matrix")}
  if(max(abs(x - t(x))) > .Machine$double.eps) 
    {if(silent) return (NULL) else stop("x appears to be not symmetric") } 
  x <- (x + t(x))/2
  u <- try(chol(x, pivot = FALSE), silent = silent)
  if(class(u) == "try-error") {
     if(silent) return(NULL) else
       stop("x appears to be not positive definite") }
  inv <- chol2inv(u)
  if(log.det) attr(inv, "log.det") <- 2 * sum(log(diag(u)))
  return(inv)
}

.onLoad <- .First.lib <- function(library, pkg)
{ 
   Rv <- R.Version()
   library.dynam("mnormt", pkg, library)
   invisible()
}



