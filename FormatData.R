########################################################################################################################
# Name        : FormatData
#
# Author      : Emixi Valdez (emixi-sthefany.valdez-medina.1@ulaval.ca)
# Date        : 12 MAR 2020
#
# Description : This code prepares the observations and precipitation forecasts in HOOPLA format(Hydromet_obs & Ens_met_fcst 
#               folders) to work directly on the CSGD method 
#               
########################################################################################################################
#----------------------------------------------------- DECLARATIONS ----------------------------------------------------
rm(list=ls())  #clear all variables
cat("\014")    #clear the console

library(abind)
library(R.matlab)
library(lubridate)

options(scipen = 999) #Disable scientific notation 

########################################################################################################################
#------------------------------------------------ THE ONLY PART TO MODIFY : -------------------------------------------#


#-----------------------------------------------------------------------------------------
setwd("C:/Users/Emi Valmed/Documents/COURS LAVAL/THESE/PROJECT/PROGRAMMING/GitHUB/CSGD/") #Path of the main folder
#-----------------------------------------------------------------------------------------

ts <- '24h'                                                             # Time step: could be 24h or 3h. Otherwise, you must 
                                                                        # create the folder named as your TimeSteps into the 
                                                                        # main folders of the program (the ones in capital letters).
                                                                           
#Forecast settings                                                                       
nbmMet <-50
nbLT <-14

#Date availables
StartDateFcast <- "2008/01/01 00:00"
FinDateFcast <-   "2016/12/15 00:00"

StartDateObs <- "1997/01/01 00:00"
FinDateObs <- "2016/12/31 00:00"

########################################################################################################################
#-------------------------------------------------- DON'T TOUCH FROM HERE ---------------------------------------------#

source("AuxiliaryFunctions.r")

#DIRECTORIES:
#------------
file_data <- paste("./ROW_DATA", ts, sep="/")
output_file <- file_data


#DOWNLOADING THE DATA (HOOPLA variables): 
#---------------------------------------

# Catchments:
pathname<-file.path(file_data,paste("catchment_names.mat",sep=""))
Catchments<-readMat(pathname)
NameC <- Catchments$nameC
nameC <- array(dim=c(length(NameC)))

for (iBV in 1:length(NameC)){
  nameC[iBV] <- NameC[[iBV]]
}
nBV <- length(nameC)


# Dates settings:
datecharFcast <- seq(ymd_hm(StartDateFcast), ymd_hm(FinDateFcast), by = "day")
if (ts=="3h") datecharObs <- seq(ymd_hm(StartDateObs), ymd_hm(FinDateObs), by = "3 hours") 
if (ts=="24h") datecharObs <- seq(ymd_hm(StartDateObs), ymd_hm(FinDateObs), by = "day")   

Index <- croping(datecharFcast, datecharObs) # To match observation's data with those of the forecast

dates <- as.numeric(paste(year(datecharFcast),   
                          sprintf("%02d",month(datecharFcast)),
                          sprintf("%02d",day(datecharFcast)), 
                          sprintf("%02d",hour(datecharFcast)),"00", sep=""))

month.string <- c("Jan","Feb","Mar","Apr","Mai","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
years <- (as.numeric(strsplit(StartDateFcast,'/')[[1]][1])):(as.numeric(strsplit(FinDateFcast,'/')[[1]][1])) # which years are concerned 
nyrs <- length(years)


# Creating arrays for storage:
Pt_Fcast <-array(dim= c(length(dates), nbmMet, nbLT, nBV))
Pt_Obs <-array(dim=c(length(dates), nbLT, nBV))


for (iBV in 1:nBV){
  
  #Extracting forecast:
  pathFcast <-file.path(paste(file_data,"/","Ens_met_fcast/Met_fcast_",nameC[iBV],".mat",sep=''))
  Fcast_Data <- readMat(pathFcast)
  dataList <- Fcast_Data$Met.fcast
  Pt_tmpFcast <- dataList[2,,]
  
  for (ibmMet in 1:nbmMet){
    Pt_Fcast[,ibmMet,,iBV] <- Pt_tmpFcast[[ibmMet]][1:length(dates),1:nbLT]
    
  }
  
  # Extracting observation:
  pathObs <-file.path(paste(file_data,"/","Hydromet_obs/Hydromet_obs_",nameC[iBV],".mat",sep=''))
  Obs_Data <- readMat(pathObs)
  
  
  for (ibLT in 1:nbLT){
    
    Pt_Obs[,ibLT,iBV] <- Obs_Data$Pt[Index+ibLT]

  }
} #nVB

save(Pt_Obs, Pt_Fcast, dates, ts, nameC, nBV, nbLT, nbmMet, month.string, years, nyrs, 
     file=paste(output_file, "/DATA_FOR_CSGD_",ts, '.Rdata',sep=''))



