#!/usr/bin/env python

'''loads example csv files sent by Santiago Begueria contianing point-wise radiation value from AEMET and converts them to netCDF format.
Since the data is a Swiss cheese (no common time axis is delivered, time slice from different stations are stacked on top of each other,
 time order is upside down), the following working stetps are accomplished:

1. Read the CSV file
2. Define a nan matrix M with i = common time instants from the earliest to latest date on which data is available (M has daily resolution)
and j = numer of stations
3. Run through each time instant and station to fill M
4. Save M in netCDF format, including appropriate metadata / attributes
Author: Swen Brands, brandssf@ifca.unican.es
'''

#load packages
import numpy as np
import xarray as xr
#import matplotlib
#matplotlib.use('WebAgg') #define backend to avoid wrong alhpa values when using the default backend
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import os
import pandas as pd
exec(open('functions_radiation.py').read())

#set input parameters
home = os.getenv('HOME')
rundir = home+'/datos/tareas/proyectos/pticlima/pyPTIclima/pySolar' #you should be here when running these scripts via ipython
dir_obs = home+'/datos/tareas/proyectos/pticlima/radiation/srcfiles' #path to input csv files provided by AEMET
dir_rean = home+'/datos/OBSData/era5_land/day/ssrd/2022' #path to reanalysis data used for comparison;currently not in use
dir_netcdf = home+'/datos/tareas/proyectos/pticlima/radiation/netcdf' #path to output netcdf file generated by this script
dir_figs = home+'/datos/tareas/proyectos/pticlima/radiation/figs' #path to output figures file generated by this script
precision = 'float32' #precision of the variable in the output netCDF files
dpival = 300 #resultion of the output figure in dpi
figformat = 'pdf' #format of the output figures: pdf, png, etc.
colormap = 'Spectral_r'

##EXECUTE ##############################################################
#read the necesary columns from the csv file
filename_csv = dir_obs + '/RadiacionDESEMON_example.csv'
df = pd.read_csv(filename_csv, sep=';', header=0, index_col=[0], encoding='latin-1')
years = df['AÑO'].values
months = df['MES'].values
days = df['DIA'].values
dates_file = pd.DatetimeIndex([str(df['AÑO'].values[ii])+'-'+str(df['MES'].values[ii]).zfill(2)+'-'+str(df['DIA'].values[ii]).zfill(2) +'T00:00:00' for ii in np.arange(len(df))])
#dates_unique = pd.DatetimeIndex(np.unique(dates_file))
dates_range = pd.date_range(dates_file.min(),dates_file.max(),freq='D')

indications = df.index.values #INDICATIVO is the index and contains the distinct station names
aemet_code, ind_unique = np.unique(indications,return_index = True)
aemet_code = aemet_code.tolist() # is a list of unique station identifiers / codes as provided by AEMET
lon = df['LONGITUD'][ind_unique].values/100000
lat = df['LATITUD'][ind_unique].values/10000
station_name = df['NOMBRE'][ind_unique].values.tolist()
altitude = df['ALTITUD'][ind_unique].values

#retain the last digit in longitude values provided by AEMET, this last digit inicates degrees East (1) or degrees West (2)
last_digit = [int(str(lon[ii])[-1]) for ii in np.arange(len(lon))]
west_ind = np.array(last_digit) == 2 #find western hemisphere longitudes
lon[west_ind] = lon[west_ind]*-1 #multiply westerly longitudes by -1 to get -180 to 180 format
#remove the last digit
lon = np.floor(lon*10000)/10000

#init output numpy arrays to be filled
rdirdia_out = np.zeros((len(dates_range),len(aemet_code)))
rdirdia_out[:] = np.nan
rdifdia_out = np.copy(rdirdia_out)
rglodia_out = np.copy(rdirdia_out)

#first loop browses through station identifiers (called 'INDICATIVO' in the csv file)
for st in np.arange(len(aemet_code)):
    #select all rows for a given station
    df_subset = df.loc[aemet_code[st]]
    dates_subset = pd.DatetimeIndex([str(df_subset['AÑO'].values[ii])+'-'+str(df_subset['MES'].values[ii]).zfill(2)+'-'+str(df_subset['DIA'].values[ii]).zfill(2) +'T00:00:00' for ii in np.arange(len(df_subset))])
    
    ##correct time ordering, in the source file time ordering is upside down !
    time_ascend_ind = np.argsort(dates_subset)
    dates_subset = dates_subset[time_ascend_ind]
    df_subset = df_subset.iloc[time_ascend_ind]
    
    ##check for duplicate time instants in input csv file
    # dupl_ind = np.where(dates_subset.duplicated())[0] #used to document the duplicated time instants at station <aemet_code[st]>
    # if len(dupl_ind) > 0:
        # print('WARNING: At station '+aemet_code[st]+', duplicated values have been dected in input file '+filename_csv+' on the dates listed below that will be set to nan:')
        # print(dates_subset.values[dupl_ind])
    # else:
        # print('INFO: At station '+aemet_code[st]+', no duplicated values have been detected in input file '+filename_csv)
    
    ##retain only time instants with no duplicates
    no_dupl_ind = np.where(~dates_subset.duplicated())[0]
    df_subset = df_subset.iloc[no_dupl_ind]
    dates_subset = dates_subset[no_dupl_ind]
    
    ##second loop browses through dates available for a specific station and then fills in the data at st and dd
    for dd in np.arange(len(dates_subset)):
        dates_subset_ind = np.where(dates_range == dates_subset[dd])[0]
        bool_dates = np.isin(dates_range,dates_subset)
        rdirdia_out[bool_dates,st] = df_subset['RDIRDIA'].values
        rdifdia_out[bool_dates,st] = df_subset['RDIFDIA'].values
        rglodia_out[bool_dates,st] = df_subset['RGLODIA'].values

#generate xr data arrays containing the original AEMET data
source_info = 'This netCDF file contains daily solar radiation time series from AEMET provided by Santiago Begueria in the framework of the PTIclima project; the original csv file contained the 3 radation variables RDIRDIA, RDIFDIA and RGLODIA.'
rglodia_out = get_xr_arr(rglodia_out, [dates_range, np.arange(len(aemet_code))], 'rglodia', 'global downward shortwave radiation', '?', altitude, station_name, aemet_code, lat, lon, source_info)
rdirdia_out = get_xr_arr(rdirdia_out, [dates_range, np.arange(len(aemet_code))], 'rdirdia', 'direct downward shortwave radiation', '?', altitude, station_name, aemet_code, lat, lon, source_info)
rdifdia_out = get_xr_arr(rdifdia_out, [dates_range, np.arange(len(aemet_code))], 'rdifdia', 'diffuse downward shortwave radiation', '?', altitude, station_name, aemet_code, lat, lon, source_info)

#generate xr data arrays containing RGLODIA data transformed to rsds
rsds = rglodia_out.copy().rename('rsds')
rsds = rsds*10000/86400 #unit conversion
rsds = rsds.astype(precision)
rsds.attrs['name'] = 'rsds'
rsds.attrs['standard_name'] = 'rsds'
rsds.attrs['long_name'] = 'surface_downwelling_shortwave_radiation'
rsds.attrs['units'] = 'W*m-2'
rsds.attrs['conversion_factor'] = 'rsds = RGLODIA*10000/86400'
#rsds.attrs['conversion_factor'] = 'rsds = RGLODIA*10000/3600'
rsds.attrs['source'] = 'This netCDF file contains rsds time series from AEMET provided by Santiago Begueria; rsds was obtained by multiplying the original RGLODIR variable with the factor indicated in the <conversion_factor> attribute.'
rsds.attrs['nan_percentage'] = np.isnan(rsds.values).sum(axis=0)/rsds.shape[0]*100

#save the output netCDF files
if os.path.isdir(dir_netcdf) != True:
    os.makedirs(dir_netcdf)

start_time = str(dates_range.min()).replace('-','').replace(' ','').replace(':','')[0:-6]
end_time = str(dates_range.max()).replace('-','').replace(' ','').replace(':','')[0:-6]
savename_rsds = dir_netcdf+'/rsds_day_aemet_'+start_time+'_'+end_time+'.nc'
print('INFO: saving re-organized AEMET station data at '+savename_rsds) 
rsds.to_netcdf(savename_rsds)

##plot the matrix of daily files (time x location)
fig = plt.figure()
ax = rsds.plot.pcolormesh(cmap = colormap, x = 'time', y = 'location', vmin = rsds.min(), vmax = rsds.max(), add_colorbar=False)
ax.axes.set_yticks(rsds.location.values)
ax.axes.set_yticklabels(station_name,fontsize=2)
plt.xticks(fontsize=5)
plt.xlabel(None)
cbar = plt.colorbar(ax,shrink=0.5,label=rsds.name + ' ('+rsds.units+')')
if figformat == 'pdf': #needed to account for irregular behaviour with the alpha parameter when plotting a pdf file
    fig.set_rasterized(True)
savename = dir_figs+'/rsds_day_aemet_coverage_'+start_time+'_'+end_time+'.'+figformat
plt.savefig(savename,dpi=dpival)
plt.close('all')

rsds.close()
print('INFO: csv2nc.py has been run successfully ! AEMET radiation data in csv format were transformed to netCDF format. The output files are located at '+dir_netcdf)
