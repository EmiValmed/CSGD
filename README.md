# CSGD
This repository contains files with R-code used to generate calibrated precipitation ensemble forecasts to be used in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA) directly. 
[CSGD](https://github.com/mscheuerer/PrecipitationFields)


The following gives a brief description of the individual files: 

* **FormatData.r**: This code prepares the observations and precipitation forecasts in [HOOPLA](https://github.com/AntoineThiboult/HOOPLA) format to work directly on the CSGD functions. It reads the Matlab variables in Hydromet_obs & Ens_met_fcst folders. 


* **AuxiliaryFunctions.r**: A collection of small functions (calculating ensemble statistics, plotting reliability diagrams, etc.) used by the other R-scripts
* **CSGD**:
* **CreateSchaakeRanks.r***: 
