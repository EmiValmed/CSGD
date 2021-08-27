########################################################################################################################
# Name        : CreateSchaakeRanks. Aapoted from Michael Scheuerer's code (paper: 10.1002/2016WR020133)
#
# Author      : Emixi Valdez (emixi-sthefany.valdez-medina.1@ulaval.ca)
# Date        : 13 Mar 2020
#
# Description : This code selects historical trajectories by the method proposed by Scheuerer et Hamill 2017 (Min. Divergence 
#               Schaake Shuffle [MDSS]) and determine the rank order.
# 
########################################################################################################################
#----------------------------------------------------- DECLARATIONS ----------------------------------------------------
rm(list=ls())  #clear all variables
cat("\014")    #clear the console

library(abind)
library(lubridate)

options(scipen = 999)  #Disable scientific notation

########################################################################################################################
#------------------------------------------------ THE ONLY PART TO MODIFY : -------------------------------------------#


#-----------------------------------------------------------------------------------------
setwd("Main/folder/path/")                                                  #Path of the main folder
#-----------------------------------------------------------------------------------------

ts <- '3h'                                                                  #Time step: could be 24h or 3h.  
HH <- 6                                                                     #Equal 6 if ts= 3h and equal 0 if ts= 24h

########################################################################################################################
#-------------------------------------------------- DON'T TOUCH FROM HERE ---------------------------------------------#


# INPUT FOLDERS:
#---------------
file_parameters <- paste(paste("PARAMETERS",ts,sep="/"),"/", sep="")                 #Parameters calculated from CSGD code.
inputs_files <- paste("./ROW_DATA", ts, sep="/")
output_file <- paste(paste("RANKS",ts,sep="/"),"/", sep="")  

## DOWNLOADING DATA:
#-------------------
load(paste(inputs_files, "/DATA_FOR_CSGD_",ts, '.Rdata',sep=''))


# Forecast settings
nbLT <- nbLT             #LeadTime
nbmMet <- nbmMet         #Meteo members

# Catchments
nameC <- nameC
nBV <- nBV


#Data
ENSEMBLE <- Pt_Fcast     # Forecast data.      Dim: (nDays* nMembers* nbLT* nBV) 
OBS <- Pt_Obs            # Observations data.  Dim: (nDays* nbLT* nBV)


#Date settings
nyrs <- nyrs
dates <- dates
years <- years
month.string <- month.string
nmonth <- length(month.string)
yyyy <- dates %/% 100000000
mm <- (dates%/%1000000)%%100
dd <- (dates%/%10000)%%100
hh <- (dates%/%100)%%100

basin.group <- list(A=1, B=2, C=3, D=4, E=5, Fe=6, G=7, H=8, I=9, J=10)	# list which associates the subbasins with one of the basins
#In case there are sub-basins --> basin.group <- list(A=1:5, B=6:10)



for (month in 1:nmonth)  { 
  
  ind.mm <- (mm == month)
  
  
  ## Precipitation parameters:
  #---------------------------
  apcp.mdss.ranks <- integer(nyrs*nBV*nbLT*31*nbmMet)
  dim(apcp.mdss.ranks) <- c(nyrs,nBV,nbLT,31,nbmMet)
  
  load(paste(file_parameters, "parameters_full_",month.string[month],'_',ts,'.Rdata',sep=''))	# load parameters for univariate MAP forecasts
  
  shape.fcst <- (mu.fcst/sigma.fcst)^2								#  here these are the predictive mean, standard dev.
  scale.fcst <- mu.fcst/shape.fcst										#  and shift parameter; convert to shape and scale
  
  k.seq <- rev(50+10*cumsum(0:8))		# number of trajectories retained after each iteration.
  
  
  for (iyear in 1:nyrs)  {
    
    cat(paste("Processing year ", 2007+iyear,"\n"))
    ind.yyyy <- yyyy %in% tail((1997:2016)[(1997:2016)!=(2007+iyear)],nbmMet)
    mid.ind <- which( ind.mm & ind.yyyy & dd==15 & hh==HH )		          # define time window for observation dates with which to start the thinning process
    window.ind <- abs(as.vector(outer(seq(-180,180,4), mid.ind, "+")))	# here: 45 days around the 15th of the month that is currently processed
    
    
    for (iday in dd[ind.mm & yyyy==(2007+iyear) & hh==HH])  {
      
      for (igr in 1:nBV) {                  # loop through the basins
        ngr <- length(basin.group[[igr]])		# number of subbasins within the current basin
        d <- ngr*nbLT                       # dimension of the univariate problem: subbasins * lead times
        
        shape <- shape.fcst[iyear,basin.group[[igr]],,iday]
        scale <- scale.fcst[iyear,basin.group[[igr]],,iday]
        shift <- shift.fcst[iyear,basin.group[[igr]],,iday]
        
        obs.traj.apcp <- array(dim=c(ngr,nbLT,length(window.ind)))
        apcp.fcst <- array(dim=c(ngr,nbLT,nbmMet))
        
        for (j in 1:ngr)  {					
          for (ibLT in 1:nbLT)  {
            
            # subset the full data set of historical observations: chosen time window shifted by lead time
            obs.traj.apcp[j,ibLT,] <- OBS[window.ind,ibLT,basin.group[[igr]][j]]
          }
          
          for (k in 1:nbmMet)  {
            # calculate a systematic sample (i.e. quantiles) of the univariate predictive distributions
            apcp.fcst[j,,k] <- pmax(shift+qgamma((k-0.5)/nbmMet,scale=scale,shape=shape),0)
          }#k
        }
        
        
        # Count the number of dates where historic precipitation are outside the 99% prediction interval
        outside.ind <- pgamma(sweep(obs.traj.apcp,c(1,2),shift,"-"), scale=scale, shape=shape)> 0.995
        nout.apcp <- apply(outside.ind,3,sum)
        
        # Retain at least 500 trajectories, based on how well they are inside the 99% prediction intervals
        thresh <- which(cumsum(tabulate(nout.apcp+1,d))>=500)[1]
        mdss.ind <- which(nout.apcp+1<=thresh)
        nss <- length(mdss.ind)
        e.fcst <- matrix(apcp.fcst,d,nbmMet)
        
        # Discard further trajectories to lower the divergence to the forecast ensemble
        #  (that's the nasty part of the code)
        for (k in k.seq)  {
          e.obs <- matrix(obs.traj.apcp[,,mdss.ind],d,nss)
          o.obs <- t(apply(e.obs,1,order))
          
          eps.m <- matrix(NA,d,nss)
          eps.p <- matrix(NA,d,nss)
          
          for (j in 1:nss)  {
            diff.obs.fcst <- sweep(e.fcst,1,e.obs[,j],'-')
            eps.m[,j] <- apply(pmax(diff.obs.fcst,0),1,mean)
            eps.p[,j] <- apply(pmax(-diff.obs.fcst,0),1,mean)
          }#j
          eps.d <- eps.m - eps.p
          
          alpha <- ((1:(nss-1))-0.5)/(nss-1)
          T1 <- apply(eps.p,2,sum)
          dvg.loo <- rep(NA,nss)
          
          for (j in 1:nss)  {
            T2 <- sapply(1:d, function(i) mean(alpha*eps.d[i,o.obs[i,o.obs[i,]!=j]]))
            dvg.loo[j] <- 2*(mean(T1[-j])+sum(T2))
          }
          mdss.ind <- mdss.ind[tail(order(dvg.loo),k)]
          nss <- length(mdss.ind)
        }#k
        
        # Historical dates have been selected (mdss.ind contains the corresponding indices)
        #  now proceed as for the standard Schaake shuffle and determine ranks of the corresponding observarions
        
        for (j in 1:ngr)  {
          iBV <- basin.group[[igr]][j]
          
          for (ibLT in 1:nbLT)  {
            apcp.obs <- OBS[window.ind[mdss.ind],ibLT,iBV]
            zeros <- (apcp.obs==0)
            # To reduce the amount of random ordering when MAP is zero, calculate weighted average MAP over the
            # entire basin and try ranking again
            if (any(!zeros))  {
              apcp.mdss.ranks[iyear,iBV,ibLT,iday,!zeros] <- rank(apcp.obs,ties.method="random")[!zeros]
            }
            for (i in which(zeros))  {
              apcp.obs[i] <- sum(OBS[window.ind[mdss.ind[i]],ibLT,iBV])
            }
            if (any(zeros))  {
              apcp.mdss.ranks[iyear,iBV,ibLT,iday,zeros] <- rank(apcp.obs[zeros],ties.method="random")
            }
            
          }#ibLT
        }#j
      }#igr
    }#iday
  }# iyear
  
  
  # Export files:
  #-------------
  save(apcp.mdss.ranks, file=paste(output_file,'apcp_mdss_ranks_',month.string[month],'_',ts,'.Rdata',sep=''))
  
}# month
