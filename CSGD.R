########################################################################################################################
# Name        : Post-processing CSGD. Adopted from Michael Scheuerer's code (paper: 10.1002/qj.2183)
#               https://github.com/mscheuerer/PrecipitationFields
#
# Author      : Emixi Valdez (emixi-sthefany.valdez-medina.1@ulaval.ca)
# Date        : 06 Avr 2018
#
# Description : This code is a simplified variant of the method proposed by Scheuerer and Hamill (2015b). We don't use 
#               neighbourhood information from grid-based forecasts.
# 
########################################################################################################################
#----------------------------------------------------- DECLARATIONS ----------------------------------------------------
rm(list=ls())  #clear all variables
cat("\014")    #clear the console

library(abind)
library(R.matlab)
library(lubridate)

options(scipen = 999)  #Disable scientific notation

########################################################################################################################
#------------------------------------------------ THE ONLY PART TO MODIFY : -------------------------------------------#


#-----------------------------------------------------------------------------------------
setwd("C:/Users/Emi Valmed/Documents/COURS LAVAL/THESE/PROJECT/PROGRAMMING/GitHUB/CSGD/") #Path of the main folder
#-----------------------------------------------------------------------------------------

ts <- '24h'                                                                                #Time step: could be: 24h, 3h


########################################################################################################################
#-------------------------------------------------- DON'T TOUCH FROM HERE ---------------------------------------------#

source("AuxiliaryFunctions.r")

# OUTPUT FOLDERS:
#----------------
inputs_files <- paste("./ROW_DATA", ts, sep="/")
file_statistics <- paste(paste("STATISTICS",ts,sep="/"),"/", sep="")
file_parameters <- paste(paste("PARAMETERS",ts,sep="/"),"/", sep="")
output_file <- paste(paste("RESULTS",ts,sep="/"),"/", sep="")


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


########################################################################################################################
#                                            Begenning post-traitment
########################################################################################################################


for (month in 1:nmonth){
  
  #--------------------------------------------------------------------------------------------------------------------#
  # § 4.a of Scheuerer et al. 2015: "ENSEMBLE STATISTICS" 
  #--------------------------------------------------------------------------------------------------------------------#
  
  for (ibLT in 1:nbLT){
    
    if (ts == '24h'){
      cleadb <- formatC((ibLT),width=3,flag="0")       # beginning of the accumulation period
      cleade <- formatC(ibLT+1,width=3,flag="0")       # end of the accumulation period
    }
    
    if (ts=='3h'){ 
      cleadb <- formatC((3*ibLT)+3,width=3,flag="0")            
      cleade <- formatC(3*(ibLT+1)+3,width=3,flag="0")         
    }
    
    cat(paste(month.string[month], " : calcul statistics ", cleadb,'h', " to ", cleade, 'h',"\n", sep=""))
    
    
    # Creating arrays for storage:
    ensmean.verif <- ensmeandiff.verif <- enspop.verif <- obs.verif <- array(dim=c(nyrs, nBV, 31)) 
    par.climo <- array(dim=c(nyrs, nBV, 3))
    ensmean.train <- ensmeandiff.train <- enspop.train <- obs.train <- array(dim=c(nyrs, nBV, 91*(nyrs-1)))  # 91 = max. number of days when surrounding the 
    Fcast <- array(dim=c(31, nbmMet, nbLT, nBV, nyrs))                                                        # 15th of the study month with +-45 days.
    
    
    # Defining training window: for every month, we use the 15th of this month and surrounding
    # 90 days (+-45) during all forecast years for training.
    # ----------------------------------------------------------------------------------------
    mid.ind <- which( ((dates%/%1000000)%%100) %in% month & ((dates%/%10000)%%100) == 15)        
    if(month==1) mid.ind <- c(mid.ind,length(dates)+15)                                         
    if(month==12) mid.ind <- c(-16,mid.ind)                                                      
    date.ind <- as.vector(outer(seq(-45,45,1), mid.ind, '+'))
    date.ind <- date.ind[date.ind >= 1  &  date.ind <= length(dates)]
    
    
    for (iyear in 1:nyrs){
      
      year <- years[iyear]
      
      cat(paste("Calculating ensemble statistics for year", year,"\n"))
      
      
      # Training and verification periods:
      #---------------------------------------------------
      
      # Training : 91 days sliding window. Cross valiadation over all years; leaving out one 
      # year at a time, train with the remaining ones, and verify the left-out year
      train.ind <- date.ind[(dates[date.ind]%/%100000000) != year]
      n.train <- length(train.ind)
      fcst.train <- Pt_Fcast[train.ind, , ibLT,]
      fcst.train[fcst.train<0] <-0
      obs.train[iyear, , 1:n.train] <- Pt_Obs[train.ind, ibLT, ] 
      
      
      # Verification : only the days of the year and month at hand:
      verif.ind <- date.ind[((dates[date.ind]%/%1000000)%%100) %in% month & (dates[date.ind]%/%100000000) == year]
      n.verif <- length(verif.ind)
      fcst.verif<- Pt_Fcast[verif.ind,  , ibLT,]
      fcst.verif[fcst.verif<0] <-0
      obs.verif[iyear, , 1:n.verif] <- Pt_Obs[verif.ind, ibLT, ] 
      Fcast[ 1:n.verif, , , ,iyear] <- Pt_Fcast[verif.ind,  , , ] # For verification purposes

            
      for (iBV in 1:nBV){
        
        cl.avg.fcst <- mean(fcst.train[, , iBV], na.rm=TRUE)    # climatology of the forecasts (all members)
        fcst.bc.train <- fcst.train [, ,iBV ] / cl.avg.fcst
        fcst.bc.verif <- fcst.verif [, ,iBV ] / cl.avg.fcst
        
        ensmean.train[iyear,iBV, 1:n.train] <- apply(fcst.bc.train, 1, mean, na.rm=TRUE)
        ensmeandiff.train[iyear, iBV, 1:n.train] <- apply(fcst.bc.train, 1, wgt.md, na.rm=TRUE)
        enspop.train[iyear, iBV, 1:n.train] <- apply(1*(fcst.bc.train > 0), 1, mean, na.rm=TRUE)
        
        ensmean.verif[iyear, iBV, 1:n.verif] <- apply(fcst.bc.verif, 1, mean, na.rm=TRUE)
        ensmeandiff.verif[iyear, iBV, 1:n.verif] <- apply(fcst.bc.verif, 1, wgt.md, na.rm=TRUE)
        enspop.verif[iyear, iBV, 1:n.verif] <- apply(1*(fcst.bc.verif > 0), 1, mean, na.rm=TRUE)
        
        
      } #iBV
    
      
      #----------------------------------------------------------------------------------------------------------------------------------#
      # § 4.b of Scheuerer et al. 2015: Fit observation climatology (CSGD distribution parameters for the unconditional distribution)
      #----------------------------------------------------------------------------------------------------------------------------------#
      
      for (iBV in 1:nBV){
        
        obs.mean <- mean(obs.train[iyear,iBV,][obs.train[iyear,iBV,] > 0], na.rm=TRUE)
        obs.pop <- mean(obs.train[iyear,iBV,] > 0, na.rm=TRUE)
        sigma <- obs.mean
        
        if (obs.pop < 0.001)  {
          par.climo[iyear,iBV,] <- c(0.0005, 0.0182, -0.00049)
          next
        }
        
        for (mu in (40:1)*(sigma/40))  {
          shape <- (mu/sigma)^2
          scale <- mu/shape
          shift <- -qgamma(obs.pop, shape=shape, scale=scale, lower.tail=FALSE)
          if (shift > -mu/2) break
        }
        par0 <- c(mu, sigma, shift)
        
        if (obs.pop < 0.01)  {
          par.climo[iyear,iBV,] <- par0
          next
        }
        par.climo[iyear,iBV,] <- optim(par0, crps.climo, obs=obs.train[iyear,iBV,], method="L-BFGS-B", lower=par0*c(0.5,0.5,2), upper=par0*c(2,2,0.1))$par
        
      } #iBV
      #---------------------------------------------------
      
      } #iyear
    
    
    save(ensmean.train, ensmeandiff.train, enspop.train, obs.train, ensmean.verif, ensmeandiff.verif, enspop.verif, obs.verif, par.climo, Fcast, 
         file=paste(file_statistics, "statistics_",month.string[month],'_',cleadb,'_',cleade,'_',ts, '.Rdata',sep=''))
  } #ibLT
  
  
  #--------------------------------------------------------------------------------------------------------------------#
  # § 4.c of Scheuerer et al. 2015: REGRESSION EQUATIONS 
  #--------------------------------------------------------------------------------------------------------------------#
  
  
  par.reg <- array(dim=c(nyrs, nBV, nbLT, 6))
  mu.fcst <- sigma.fcst <- shift.fcst <- obs <- crps.save <- array(dim=c(nyrs, nBV, nbLT, 31))
  
  for (ibLT in 1:nbLT){
    
    if (ts == '24h'){
      cleadb <- formatC((ibLT),width=3,flag="0")       # beginning of the accumulation period
      cleade <- formatC(ibLT+1,width=3,flag="0")       # end of the accumulation period
    }
    
    if (ts=='3h'){ 
      cleadb <- formatC((3*ibLT)+3,width=3,flag="0")            
      cleade <- formatC(3*(ibLT+1)+3,width=3,flag="0")         
    }
    
    
    cat(paste(month.string[month], " : calcul regression ", cleadb, " to ", cleade,"\n", sep=""))
    load(paste(file_statistics, "statistics_",month.string[month],'_',cleadb,'_',cleade,'_',ts,'.Rdata',sep=''))
    
    for (iyear in 1:nyrs){
      
      year <- years[iyear]
      
      for (iBV in 1:nBV){
        
        par0 <- c(0.1,0.1,1.0,1.0,0.5,0.5)
        
        opt.res <- optim(par = par0,
                         fn = crps.reg,
                         obs = obs.train[iyear,iBV,],
                         enspop = enspop.train[iyear,iBV,],
                         ensmean = ensmean.train[iyear,iBV,], 
                         ensmeandiff = ensmeandiff.train[iyear,iBV,],
                         par.climo = par.climo[iyear,iBV,],
                         method = "L-BFGS-B",
                         lower = c(0.001, 0.05, 0.0, 0.0, 0.1, 0.0),
                         upper = c(1.0,    1.0, 1.5, 1.5, 1.0, 1.5))
        
        
        par.reg[iyear,iBV,ibLT,] <- opt.res$par
        
        mu.cl    <- par.climo[iyear,iBV,1]
        sigma.cl <- par.climo[iyear,iBV,2]
        shift.cl <- par.climo[iyear,iBV,3]
        
        log.arg <- opt.res$par[2] + opt.res$par[3]*enspop.verif[iyear,iBV,] + opt.res$par[4]*ensmean.verif[iyear,iBV,]
        mu.fcst[iyear,iBV,ibLT,] <- mu.cl*log1p(expm1(opt.res$par[1])*log.arg) / opt.res$par[1]
        sigma.fcst[iyear,iBV,ibLT,] <- opt.res$par[5]*sigma.cl*sqrt(mu.fcst[iyear,iBV,ibLT,]/mu.cl) + opt.res$par[6]*sigma.cl*ensmeandiff.verif[iyear,iBV,]
        shift.fcst[iyear,iBV,ibLT,] <- rep(shift.cl,length(ensmean.verif[iyear,iBV,]))
        
        obs[iyear,iBV,ibLT,] <- obs.verif[iyear,iBV,]
        
        
        # CRPS and MCRPS (by month)
        #---------------------------------------------
        v_crps <- crps.reg.complete(par = opt.res$par, 
                                    obs = obs.verif[iyear,iBV,], 
                                    enspop = enspop.verif[iyear,iBV,], 
                                    ensmean = ensmean.verif[iyear,iBV,], 
                                    ensmeandiff = ensmeandiff.verif[iyear,iBV,], 
                                    par.climo = par.climo[iyear,iBV,])
        crps.save[iyear,iBV,ibLT,1:length(v_crps)] <- v_crps
        
        
      } #iBV
    } #iyear
    
  } #ibLT
  
  save(par.reg, mu.fcst, sigma.fcst, shift.fcst, obs, crps.save, Fcast,
       file=paste(file_parameters, "parameters_full_",month.string[month],'_',ts,'.Rdata',sep=''))
} #month


########################################################################################################################
#                                            Post-traitment finished :) !!
########################################################################################################################






