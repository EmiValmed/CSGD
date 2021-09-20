# CSGD :umbrella: :earth_americas: :computer:
This repository contains files to generate calibrated precipitation ensemble forecasts to be used in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA) using the [CSGD](https://journals.ametsoc.org/view/journals/mwre/143/11/mwr-d-15-0061.1.xml) method.

The [CSGD](https://journals.ametsoc.org/view/journals/mwre/143/11/mwr-d-15-0061.1.xml) post-processor is based on a complex heteroscedastic, nonlinearregression model conceived to address the peculiarities of precipitation (e.g., its intermittent and highly skewed nature andits typically large forecast errors). This method yields full predictive probability distributions for precipitation accumulations based on ensemble model output statistics (EMOS) and censored, shifted gamma distributions. 

[HOOPLA](https://github.com/AntoineThiboult/HOOPLA) is an automatic software that allows to carry out model calibration and obtain hydrological simulations and forecasts at different time steps.

The following gives a brief description of the individual files: 

* **FormatData.r**: This code prepares the observations and precipitation forecasts in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA) format to work directly on the CSGD functions. It reads the Matlab variables in **Hydromet_obs** & **Ens_met_fcst folders**. 
* **AuxiliaryFunctions.r**: A collection of small functions (calculating ensemble statistics, data selection, etc.) used by the other R-scripts
* **CSGD.r**: This code is a simplified variant of the method proposed by Scheuerer and Hamill [2015b](https://journals.ametsoc.org/view/journals/mwre/143/11/mwr-d-15-0061.1.xml): Censored, Shifted Gamma Distributions (CSGD).  
* **CreateSchaakeRanks.r**: This code selects historical trajectories by the method proposed by Scheuerer et al. [2017](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR020133) and determine the rank order: Min. Divergence Schaake Shuffle (MDSS).


## How does it work? :memo:

### Preliminar steps

1. Copy and paste the **Hydromet_obs** & **Ens_met_fcst** folders from [HOOPLA](https://github.com/AntoineThiboult/HOOPLA) into the **RAW_DATA/time step** folder (time step could be 3h or 24h). Folders path in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA):

    >  **Observations:** HOOPLA-master/Data/time step/Hydromet_obs
    
    >  **Forecasts:** HOOPLA-master/Data/time step/Ens_met_fcst 
       
2. Copy and paste the **catchment_names.mat** file into the **RAW_DATA/time step** folder. File path in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA):

    > **Catchments names:** HOOPLA-master/Data/time step/Misc

### Execute the functions in the following order:

1. **AuxiliaryFunctions.r** 
2. **FormatData.r**
3. **CSGD.r**
4. **CreateSchaakeRanks.r**
    
Before executing the functions, modify the forecasts/observations settings as desired. All parts that can be modified in the scripts are indicated ("THE ONLY PART TO MODIFY"):

* setwd("C:/Main/Folder/Path/")  
* ts = time step
* nbmMet = No. of meteorological ensemble members
* nbLT = lead times
* StartDateFcast = first day of the forecast period
* FinDateFcast = last day of the forecast period
* StartDateObs = first day of the observation period
* FinDateObs = last day of the observation period
* HH = Define time window for observation dates. Equal 6 if ts= 3h and equal 0 if ts= 24h 

