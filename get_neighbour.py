#!/usr/bin/env python

'''Retains reanalysis data at the grid-boxes located nearest to the AEMET station data, that was read-in and ordered by csv2nc.py in a previous working step.
Hourly reanlaysis data is aggregated to daily-mean values as indicated by the <temporal_aggregation> attribute defined below.
Author: Swen Brands, brandssf@ifca.unican.es
'''

#load packages
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import os
import pandas as pd
import xskillscore as xs
from math import radians, cos, sin, asin, sqrt #needed to calculate haversine distance
home = os.getenv('HOME')
exec(open('functions_radiation.py').read())
exec(open(home+'/datos/tareas/proyectos/pticlima/seasonal/python/functions_seasonal.py').read())

#set input parameters
model_dataset = 'era5_land' #model dataset to be processed, ERA5-Land or ERA5
rundir = home+'/datos/tareas/proyectos/pticlima/pyPTIclima/pySolar' #script directory, you should be there or point to this directory when running these scripts via python
dir_obs = home+'/datos/tareas/proyectos/pticlima/radiation/netcdf' #path to input netcdf files produced with csv2nc.py containg AEMET station data.
filename_obs = 'rsds_day_aemet_20171201_20221231.nc'
dir_rean = home+'/datos/OBSData/'+model_dataset #path to reanalysis data used for comparison; base directory structure similar to data, will be expanded as a function of the target area (Iberia or Canarias) and year
#dir_rean = '/media/swen/ext_disk2/datos/OBSdata/'+model_dataset 
dir_netcdf = home+'/datos/tareas/proyectos/pticlima/radiation/results/data' #path to output netcdf file generated by this script
dir_figs = home+'/datos/tareas/proyectos/pticlima/radiation/results/validation' #path to output figures file generated by this script
taryears = [2017,2022] #start and end years used for validation
variable = 'ssrd' #variable names are harmonized to ERA5 standard by csv2nc.py
variable_aemet = 'rsds' #variable name after transformation

precision = 'float32' #precision of the variable in the output netCDF files
dpival = 300 #resultion of the output figure in dpi
figformat = 'pdf' #format of the output figures: pdf, png, etc.
colormap = 'Spectral_r'

##EXECUTE ##############################################################
years = np.arange(taryears[0],taryears[1]+1)
listdir_ca = get_nc_path(years,'Canarias',dir_rean,variable)
listdir_sp = get_nc_path(years,'Iberia',dir_rean,variable)
print('The following files will be loaded for the Canarias domain:')
print(listdir_ca)
print('The following files will be loaded for the Iberia domain:')
print(listdir_sp)
nc_ca = xr.open_mfdataset(listdir_ca)
nc_sp = xr.open_mfdataset(listdir_sp)

#create pandas Datetime indices
dates_sp = pd.DatetimeIndex(nc_sp.time.values)
dates_ca = pd.DatetimeIndex(nc_ca.time.values)

#get 00 UTC value for each day, containing the accumlated value of the prior 24 hours in ERA5-Land, then transform toe W/m^2, see https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790 and https://confluence.ecmwf.int/pages/viewpage.action?pageId=155337784
if model_dataset == 'era5_land':
    #hourly ERA5-Land data from CDS are accumulated daily from 00 UTC to the hour indicated in the file and then set to 0 at 00UTC of the next day, when the accumulation starts once again.
    acc_ind_sp = np.where(dates_sp.hour == 0)[0]
    acc_ind_ca = np.where(dates_ca.hour == 0)[0]
    nc_sp = nc_sp.isel(time=acc_ind_sp)/86400
    nc_ca = nc_ca.isel(time=acc_ind_ca)/86400
    #create pandas Datetime indices and shift by one day backward to tack into account the daily accumulation period in case of ERA5-Land (the 00 UTC value at date t contains the 24 hour flux accumlation of the day before, i.e. t-1 day)
    dates_sp = pd.DatetimeIndex(nc_sp.time.values).shift(-1,'D')
    dates_ca = pd.DatetimeIndex(nc_ca.time.values).shift(-1,'D')
    #redefine time dimension and variable name
    nc_sp['time'] = dates_sp
    nc_ca['time'] = dates_ca
elif model_dataset == 'era5':
    #cacluate daily mean values, in contrast to ERA5-Land hourly ERA5 data are normally accumulated for each hour, the value of 00 UTC containing the accumulation of the previous hour, i.e. from 23 UTC of the day before to 00 UTC of the current day. 
    nc_ca['time'] = dates_ca.shift(-1,'H')
    nc_sp['time'] = dates_sp.shift(-1,'H')
    nc_ca = nc_ca.resample(time="D").sum(dim="time")/86400
    nc_sp = nc_sp.resample(time="D").sum(dim="time")/86400
    dates_ca = pd.DatetimeIndex(nc_ca.time.values)
    dates_sp = pd.DatetimeIndex(nc_sp.time.values)
    nc_ca[variable][0,:,:] = np.nan #set first daily value to nan because it was accumulated with only one hour (23 UTC of the day previous to the first day to 00 UTC of the first day) 
    nc_sp[variable][0,:,:] = np.nan
    nc_ca[variable][-1,:,:] = np.nan #set last daily value to nan because the last hour of accumulation is missing (23 UTC to 24 UTC / 00 UTC of the next day) for the accumulated data ending at 23 UTC of the last day.
    nc_sp[variable][-1,:,:] = np.nan
else:
    raise Exception('ERROR: unknown enty for <model_dataset>!')
    
nc_sp = nc_sp[variable].rename(variable_aemet)
nc_ca = nc_ca[variable].rename(variable_aemet)

#load AEMET station data and corresponding metadata
obsfile = dir_obs+'/'+filename_obs
nc_obs = xr.open_dataset(obsfile)
altitude = nc_obs.rsds.location.altitude
station_name = nc_obs.rsds.location.station_name
aemet_code = nc_obs.rsds.location.aemet_code
lat_obs = nc_obs.rsds.location.latitude
lon_obs = nc_obs.rsds.location.longitude
dates_obs = pd.DatetimeIndex(nc_obs.time.values)

#check whether reanalysis dates for the 2 regions (Iberia and Canarias) are identical, if yes, use only one date object thereafter (dates_rean)
if np.all(dates_sp.isin(dates_ca)) != True:
    raise Exception('ERROR: Reanalysis dates for the two regions SP and CA are not identical !')
dates_rean = dates_sp
del(dates_sp,dates_ca)

#get common time period
ind_dates_rean = np.where(dates_rean.isin(dates_obs))[0]
ind_dates_obs = np.where(dates_obs.isin(dates_rean))[0]
nc_sp = nc_sp[ind_dates_rean]
nc_ca = nc_ca[ind_dates_rean]
dates_rean = dates_rean[ind_dates_rean]
#nc_obs[variable_aemet] = nc_obs.rsds[ind_dates_obs]
nc_obs = nc_obs.isel(time=ind_dates_obs)
dates_obs = dates_obs[ind_dates_obs]
nc_sp.values
nc_ca.values
nc_obs.values

# #get nearest neighbour indices
lat_sp = nc_sp.latitude.values
lon_sp = nc_sp.longitude.values
lat_ca = nc_ca.latitude.values
lon_ca = nc_ca.longitude.values

#get nearest neighbour values (neighs) from reanalysis
nanmask_sp = np.transpose(np.isnan(nc_sp.values).sum(axis=0)/nc_sp.values.shape[0])
nanmask_ca = np.transpose(np.isnan(nc_ca.values).sum(axis=0)/nc_ca.values.shape[0])
neighs = np.zeros(nc_obs.rsds.shape)
lat_neighs = np.zeros(nc_obs.rsds.shape[1])
lon_neighs = np.copy(lat_neighs) 
for st in np.arange(nc_obs.rsds.shape[1]):
    print('INFO: processing '+nc_obs.location.station_name[st]+'...')
    #load nearest neighbour reanlaysis data as a function of the domain, either Iberian Peninsula (IB) or Canary Islands (CA)
    if nc_obs.location.domain[st] == 'IB':
        #nc_point = nc_sp.sel(latitude=lat_obs[st], longitude=lon_obs[st], method = 'nearest')
        dist_sp = np.zeros((len(lon_sp),len(lat_sp)))
        for xx in np.arange(len(lon_sp)):
            for yy in np.arange(len(lat_sp)):
                dist_sp[xx,yy] = haversine(lon_obs[st], lat_obs[st], lon_sp[xx], lat_sp[yy])
        dist_sp[nanmask_sp==1] = np.nan
        minind = np.where(dist_sp == np.nanmin(dist_sp)) #2d index pointing to nearest grid box in reanalysis that is not entirely nan (i.e. ocean values in ERA5-Land are excluded)
        nc_point = nc_sp.sel(latitude=lat_sp[minind[1]], longitude=lon_sp[minind[0]])
        lat_neighs[st] = lat_sp[minind[1]]
        lon_neighs[st] = lon_sp[minind[0]]
    elif nc_obs.location.domain[st] == 'CA':
        #nc_point = nc_ca.sel(latitude=lat_obs[st], longitude=lon_obs[st], method = 'nearest')
        dist_ca = np.zeros((len(lon_ca),len(lat_ca)))
        for xx in np.arange(len(lon_ca)):
            for yy in np.arange(len(lat_ca)):
                dist_ca[xx,yy] = haversine(lon_obs[st], lat_obs[st], lon_ca[xx], lat_ca[yy])
        dist_ca[nanmask_ca==1] = np.nan
        minind = np.where(dist_ca == np.nanmin(dist_ca)) #2d index pointing to nearest grid box in reanalysis that is not entirely nan (i.e. ocean values in ERA5-Land are excluded)
        nc_point = nc_ca.sel(latitude=lat_ca[minind[1]], longitude=lon_ca[minind[0]])
        lat_neighs[st] = lat_ca[minind[1]]
        lon_neighs[st] = lon_ca[minind[0]]
    else:
        raise Exception('ERROR: Unknown domain! Please check <location.domain> attibute in the input netCDF file containing AEMET station data !') 
    neighs[:,st] = nc_point.values.squeeze()
    nc_point.close()
    #del nc_point

nc_sp.close()
nc_ca.close()

#create xarray data array
altitude_neigh = 'Will be filled in future versions'
station_name_neigh = ['nearest neighbour grid-box corresponding to '+ii for ii in station_name]
aemet_code_neigh = ['nearest neighbour grid-box corresponding to '+ii for ii in aemet_code]
source_info = 'This netCDF file contains daily solar radiation time series from '+model_dataset+' provided by Copernicus Data Store.'
neighs = get_xr_arr(neighs, [dates_rean, np.arange(neighs.shape[1])], variable_aemet, 'global downward shortwave radiation', 'W/m2', altitude_neigh, station_name_neigh, aemet_code_neigh, lat_neighs, lon_neighs, source_info)
neighs.attrs['temporal_aggregation'] = 'daily mean data caclulated upon hourly '+variable_aemet+' data from '+str(dates_obs.hour.values.min())+' to '+str(dates_obs.hour.values.min())+ 'o clock'
start_time = str(dates_obs.min()).replace('-','').replace(' ','').replace(':','')[0:-6]
end_time = str(dates_obs.max()).replace('-','').replace(' ','').replace(':','')[0:-6]
savename_neigh = dir_netcdf+'/rsds_day_'+model_dataset+'_nn_aemet_'+start_time+'_'+end_time+'.nc' #"nn" in the output file name refers to "nearest neighbour"
neighs.to_netcdf(savename_neigh)

##start verification
obs = nc_obs[variable_aemet]
##calculalate hindcast correlation coefficient for the inter-annual seasonal-mean time series (observations vs. ensemble mean) and corresponding p-values based on the effective sample size
pearson_r = xs.pearson_r(obs,neighs,dim='time',skipna=True).rename('pearson_r')
pearson_pval = xs.pearson_r_p_value(obs,neighs,dim='time',skipna=True).rename('pearson_pval')
pearson_pval_effn = xs.pearson_r_eff_p_value(obs,neighs,dim='time',skipna=True).rename('pearson_pval_effn')
spearman_r = xs.spearman_r(obs,neighs,dim='time',skipna=True).rename('spearman_r')
spearman_pval = xs.spearman_r_p_value(obs,neighs,dim='time',skipna=True).rename('spearman_pval')
spearman_pval_effn = xs.spearman_r_eff_p_value(obs,neighs,dim='time',skipna=True).rename('spearman_pval_effn')

nc_obs.close()
neighs.close()
print('INFO: get_neighbour.py has been run successfully !')
