# CSGD :umbrella: :earth_americas: :computer:
This repository contains files to generate calibrated precipitation ensemble forecasts to be used in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA) directly. 
[CSGD](https://github.com/mscheuerer/PrecipitationFields).

The following gives a brief description of the individual files: 

* **FormatData.r**: This code prepares the observations and precipitation forecasts in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA) format to work directly on the CSGD functions. It reads the Matlab variables in Hydromet_obs & Ens_met_fcst folders. 
* **AuxiliaryFunctions.r**: A collection of small functions (calculating ensemble statistics, data selection, etc.) used by the other R-scripts
* **CSGD**: This code is a simplified variant of the method proposed by Scheuerer and Hamill (2015b): Censored, Shifted Gamma Distributions (CSGD).  
* **CreateSchaakeRanks.r***: This code selects historical trajectories by the method proposed by Scheuerer et Hamill 2017: Min. Divergence Schaake Shuffle (MDSS) and determine the rank order.


## How does it work? :memo:

### Preliminar steps

1. Copy and paste the **Hydromet_obs** & **Ens_met_fcst** folders from [HOOPLA](https://github.com/AntoineThiboult/HOOPLA) into the **RAW_DATA/time step** folder (time step could be 3h or 24h). Folders path in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA):
    >  **Observations:** HOOPLA-master/Data/time step/Hydromet_obs
    >  **Forecasts:** HOOPLA-master/Data/time step/Ens_met_fcst 
2. Copy and paste the **catchment_names.mat** file into the **RAW_DATA/time step** folder. File path in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA):
    > **Catchments names:** HOOPLA-master/Data/time step/Misc

### 
Execute the functions in the following order:
1. **AuxiliaryFunctions.r** 
2. **FormatData.r**
3. **CSGD**
4. **CreateSchaakeRanks.r
    
Before executing the functions, modify the forecasts/observations settings as desired. All parts that can be modified are indicated ("THE ONLY PART TO MODIFY"):
    * setwd("C:/Main/Folder/Path/")  
    * ts = time step
    * nbmMet = No. of ensemble members
    * nbLT = lead times
    * StartDateFcast = 
    * FinDateFcast
    StartDateObs
    FinDateObs
    HH = Equal 6 if ts= 3h and equal 0 if ts= 24h define time window for observation dates


