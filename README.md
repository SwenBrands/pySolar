# pySolar
A Python package for the evaluation of surface downwelling shortwave
 radiation time series from model and reanalysis datasets
 against station observations. Each script is described by the comments
 included within. Tested under Python 3.10
 
dependecies: xarray, xskillscore, pandas, numpy, matplotlib, cartopy, os, math

author: Swen Brands (CSIC-UC), brandssf@ifca.unican.es or swen.brands@gmail.com

########################################################################

Description of the files contained in this folder

########################################################################
Each script has a specific task. The scripts have to be run sequentially 
in the following order:

1. download_era5_datasets_1h_pti_rad.py: script to download ERA5 and
ERA5-Land data from CDS.

2. csv2nc.py: Reads the csv.file containing daily AEMET station data and brings
it into the right temporal and spatial order. Generates an output
netCDF file including a 2d data array with the dimensions time x location
containing in-situ global radiation time series at 60 AEMET stations 
in mainland Spain and the Canary Islands providing daily average radiation
intensities in W/m2. The considered time period currently is 12/2017 to 12/2022.
This netCDF file comes with plenty of metadata and should be self-describing.

3. get neighbour.py: Postprocesses the reanalysis / model data so that
it can be compared with the observations stored in point 2. First, the
model / reanalyis data is disaggregated to daily average values in W/m2
and outliers, if present, are set to nan. Then, the time-series at
the grid-box nearest to the respective station location is identified
 and the corresponding data retrieved. A netCDF output file is 
generated which also comes with metadata and is thus easy to understand.

4. validate.py: Validates the model / reanalysis data from point 3. with the
observations stored in point 2. using xskillscore and stores the results
in a netCDF file. One file per dataset is generated. The script currently
validates ERA5 and ERA5-Land.

5. plot_validation_results.py: Plots the verification results obtained
in point 4.

6. disagg_and_get_neighbour.py: Currently experimental, disaggrates
hourly ERA5-Land data aggregated from 00 UTC to the forecast step
(i.e. the accumulation period is reset every 00 UTC step) to normal
hourly accumulations ending at the indicated hour, i.e. the value at 01
UTC contains the accumulation from the previous hour (00 to 01 UTC).
Currently does not interfer with scripts 1 to 5 mentioned above.

7. functions_radiation.py: Contains specific functions used by the 
scripts contained in this folder.

8. launchme.sh: bash script to launch the aforementioned Python scripts
from the shell.


Credits
-------
This ongoing research work is being funded by the Ministry for Ecological Transition and Demographic Challenge (MITECO) and the European Commission NextGenerationEU (Regulation EU 2020/2094), through CSIC's Interdisciplinary Thematic Platform Clima (PTI-Clima).

![alt text](https://pti-clima.csic.es/wp-content/uploads/2023/11/Web-Gob-Min-CSIC-COLOR-LOGO-PNG-RGB-300pppCLIMA.png)
