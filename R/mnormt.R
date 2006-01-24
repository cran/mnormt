
# if(.Platform$OS.type=="unix") dyn.load("mnormt/src/sadmvnt.so")
# if(.Platform$OS.type=="windows") dyn.load("mnormt/src/sadmvnt.dll")


dmnorm <- function(x, mean=rep(0,d), varcov, log=FALSE)
{
  d  <- if(is.matrix(varcov)) ncol(varcov) else 1
  if(d>1 & is.vector(x)) x <- matrix(x, 1, d)
  n  <- if(d==1)  length(x) else nrow(x) 
  X  <- t(matrix(x, nrow=n, ncol=d)) - mean
  Q  <- apply((solve(varcov)%*% X)* X, 2, sum) 
  logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
  logPDF <- as.vector(Q + d*logb(2*pi)+logDet)/(-2)
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
   sadmvn(lower=rep(-Inf, length(x)), upper=x, mean, varcov, ...)
 

sadmvn <- function(lower, upper, mean, varcov,  
                   maxpts=2000*d, abseps=1e-6, releps=0)
{
  if(any(lower > upper)) stop("lower>upper integration limits")
  if(any(lower == upper)) return(0)
  d <- as.integer(if(is.matrix(varcov)) ncol(varcov) else 1)
  varcov <- matrix(varcov, d, d)
  sd  <- sqrt(diag(varcov))
  rho <- diag(1/sd,d,d) %*% varcov %*% diag(1/sd,d,d)
  lower <- as.double((lower-mean)/sd)
  upper <- as.double((upper-mean)/sd)
  infin <- rep(2,d)
  infin <- replace(infin, (upper == Inf) & (lower > -Inf), 1)
  infin <- replace(infin, (upper < Inf) & (lower == -Inf), 0)
  infin <- replace(infin, (upper == Inf) & (lower == -Inf), -1)
  infin <- as.integer(infin)
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
  d  <- if(is.matrix(S)) ncol(S) else 1
  if(d>1 & is.vector(x)) x <- matrix(x, 1, d)
  n  <- if(d==1) length(x) else nrow(x) 
  X <- t(matrix(x, nrow = n, ncol = d)) - mean
  Q <- apply((solve(S) %*% X) * X, 2, sum)
  logDet <- sum(logb(abs(diag(qr(S)$qr))))
  logPDF <- (lgamma((df + d)/2) - 0.5 * (d * logb(pi * df) + logDet)
             - lgamma(df/2) - 0.5 * (df + d) * logb(1 + Q/df))
  if(log) logPDF else exp(logPDF)
}

rmt <- function(n=1, mean=rep(0,d), S, df=Inf)
{ 
  d <- if(is.matrix(S)) ncol(S) else 1 
  if(df==Inf) x <-1 else x <- rchisq(n,df)/df
  z <- rmnorm(n, rep(0,d), S)
  y <- t(mean + t(sqrt(x)*z))
  return(y)
}


pmt <- function(x, mean=rep(0,length(x)), S, df=Inf, ...)
   sadmvt(df, lower=rep(-Inf, length(x)), upper=x, mean, S, ...)


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
  rho <- diag(1/sd,d,d) %*% S %*% diag(1/sd,d,d)
  lower <- as.double((lower-mean)/sd)
  upper <- as.double((upper-mean)/sd)
  infin <- rep(2,d)
  infin <- replace(infin, (upper == Inf) & (lower > -Inf), 1)
  infin <- replace(infin, (upper < Inf) & (lower == -Inf), 0)
  infin <- replace(infin, (upper == Inf) & (lower == -Inf), -1)
  infin <- as.integer(infin)
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

.First.lib <- function(library, pkg)
{ 
   Rv <- R.Version()
   if(Rv$major == 0 |(Rv$major == 2 && Rv$minor < 2.0))
        stop("This package requires R 2.2.0 or later")
   fcode <- file.path("../libs", paste("mnormt", .Platform$dynlib.ext, sep=""))
   library.dynam(fcode)
   # if(interactive()) cat("Library 'mnormt'\n")
   invisible()
}



