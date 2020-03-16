########################################################################################################################
# Name         : AuxiliaryFunctions 
#
#
# Author       : Michael Scheuerer (michael.scheuerer@noaa.gov)
# 
#
# Modification : croping funtion
#                Modified version of wgt.md in order not to use the weights
#
#
# By           : Emixi Valdez (emixi-sthefany.valdez-medina.1@ulaval.ca)
# Date         : 12 MAR 2020
#
# Description  : A collection of small functions used by the other R-scripts. 
########################################################################################################################

library(Hmisc)

croping <- function(Date1, Date2, na.rm=FALSE)  {
  # Date1: date of forecast
  # Date2: date of observation

  datesDate1 <- as.numeric(paste(year(Date1),
                                 sprintf("%02d",month(Date1)),
                                 sprintf("%02d",day(Date1)), 
                                 sprintf("%02d",hour(Date1)),"00", sep="")) 
  
  datesDate2 <- as.numeric(paste(year(Date2),
                                 sprintf("%02d",month(Date2)),
                                 sprintf("%02d",day(Date2)), 
                                 sprintf("%02d",hour(Date2)),"00", sep=""))
  
  return(which( datesDate2 %in% datesDate1))
  
}


########################################################################################################################
#                                           CSGD FUNCTIONS
########################################################################################################################



pctg <- function(x, mu, sigma, shift)  {
  return(pgamma(x-shift, scale=(sigma^2)/mu, shape=(mu/sigma)^2))
}


qctg <- function(x, mu, sigma, shift)  {
  pmax(0, shift + qgamma( x, scale=(sigma^2)/mu, shape=(mu/sigma)^2) )
}


gini.md <- function(x, na.rm=FALSE)  {
  if(na.rm & any(is.na(x)))  x <- x[!is.na(x)]
  n <-length(x)
  return(4*sum((1:n)*sort(x, na.last=TRUE))/(n^2)-2*mean(x)*(n+1)/n)
}

wgt.md <- function(x, na.rm=FALSE)  {
  if(na.rm)  x <- x[!is.na(x)]
  x.ord <- order(as.vector(x))
  x.sort <- x[x.ord]
  N <- length(x.sort)
  w <- rep(1/N,N)
  W <- cumsum(w[x.ord])
  2*sum(W[-N]*(1-W[-N])*diff(x.sort))
}


crps.climo <- function(par, obs)  {
  obs <- sort(obs[!is.na(obs)])
  n <- length(obs)
  k0 <- sum(obs==0)
  
  shape <- (par[1]/par[2])^2
  scale <- par[1]/shape
  shift <- par[3]
  
  crps <- numeric(n)
  
  c.std <- -shift/scale
  y.std <- (obs[(k0+1):n]-shift)/scale
  
  F.k.c <- pgamma(c.std, shape=shape)
  F.kp1.c <- pgamma(c.std, shape=shape+1)
  F.2k.2c <- pgamma(2*c.std, shape=2*shape)
  B.05.kp05 <- beta(0.5,shape+0.5)
  
  F.k.y <- pgamma(y.std, shape=shape)
  F.kp1.y <- pgamma(y.std, shape=shape+1)
  
  crps[1:k0]     <- c.std*(2*F.k.c-1) - shape*(2*F.kp1.c-1+F.k.c^2-2*F.kp1.c*F.k.c) - c.std*F.k.c^2 - (shape/pi)*B.05.kp05*(1-F.2k.2c)
  crps[(k0+1):n] <- y.std*(2*F.k.y-1) - shape*(2*F.kp1.y-1+F.k.c^2-2*F.kp1.c*F.k.c) - c.std*F.k.c^2 - (shape/pi)*B.05.kp05*(1-F.2k.2c)
  
  return( scale*mean(crps) )
}

crps.reg <- function(par, obs, enspop, ensmean, ensmeandiff, par.climo)  {
  
  miss <- is.na(obs)
  
  mu.cl    <- par.climo[1]
  sigma.cl <- par.climo[2]
  shift.cl <- par.climo[3]
  
  log.arg <- par[2] + par[3]*enspop[!miss] + par[4]*ensmean[!miss]
  mu      <- mu.cl*log1p(expm1(par[1])*log.arg)/par[1]
  sigma   <- par[5]*sigma.cl*sqrt(mu/mu.cl) + par[6]*sigma.cl*ensmeandiff[!miss]
  shift   <- shift.cl
  
  scale <- sigma^2/mu
  shape <- (mu/sigma)^2
  
  y <- obs[!miss]
  c.std <- -shift/scale
  y.std <- (y-shift)/scale
  
  F.k.c <- pgamma(c.std, shape=shape)
  F.kp1.c <- pgamma(c.std, shape=shape+1)
  F.2k.2c <- pgamma(2*c.std, shape=2*shape)
  B.05.kp05 <- beta(0.5,shape+0.5)
  
  F.k.y <- pgamma(y.std, shape=shape)
  F.kp1.y <- pgamma(y.std, shape=shape+1)
  
  crps <- y.std*(2*F.k.y-1) - shape*(2*F.kp1.y-1+F.k.c^2-2*F.kp1.c*F.k.c) - c.std*F.k.c^2 - (shape/pi)*B.05.kp05*(1-F.2k.2c)
  
  return( mean(scale*crps) )
}


crps.reg.complete <- function(par, obs, enspop, ensmean, ensmeandiff, par.climo)  {
  
  #Individual values of CRPS
  
  miss <- is.na(obs)
  
  mu.cl    <- par.climo[1]
  sigma.cl <- par.climo[2]
  shift.cl <- par.climo[3]
  
  log.arg <- par[2] + par[3]*enspop[!miss] + par[4]*ensmean[!miss]
  mu      <- mu.cl*log1p(expm1(par[1])*log.arg)/par[1]
  sigma   <- par[5]*sigma.cl*sqrt(mu/mu.cl) + par[6]*sigma.cl*ensmeandiff[!miss]
  shift   <- shift.cl
  
  scale <- sigma^2/mu
  shape <- (mu/sigma)^2
  
  y <- obs[!miss]
  c.std <- -shift/scale
  y.std <- (y-shift)/scale
  
  F.k.c <- pgamma(c.std, shape=shape)
  F.kp1.c <- pgamma(c.std, shape=shape+1)
  F.2k.2c <- pgamma(2*c.std, shape=2*shape)
  B.05.kp05 <- beta(0.5,shape+0.5)
  
  F.k.y <- pgamma(y.std, shape=shape)
  F.kp1.y <- pgamma(y.std, shape=shape+1)
  
  crps <- y.std*(2*F.k.y-1) - shape*(2*F.kp1.y-1+F.k.c^2-2*F.kp1.c*F.k.c) - c.std*F.k.c^2 - (shape/pi)*B.05.kp05*(1-F.2k.2c)
  
  return(scale*crps)
}

crps.normal <- function(par, obs, ensmeanano, ensvar, obs.climo)  {
  miss <- is.na(obs)
  
  mu <- obs.climo[!miss] + par[1]*ensmeanano[!miss]
  sigma <- sqrt( par[2] + par[3]*ensvar[!miss] )
  
  obs.stdz <- (obs[!miss]-mu)/sigma
  crps <- sigma * ( obs.stdz*(2*pnorm(obs.stdz)-1) + 2*dnorm(obs.stdz) - 1/sqrt(pi) ) 
  
  return( mean(crps) )
}
