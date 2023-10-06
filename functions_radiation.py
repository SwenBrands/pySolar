#!/usr/bin/env python

'''functions used by the scripts located in .../pticlima/radiation/python'''

def get_xr_arr(np_arr_f, coords_f, varname_f, longname_f, units_f, altitude_f, station_name_f, aemet_code_f, latitude_f, longitude_f):
    
    '''Obtain xarray DataArray from input parameters and set attributes. Input: <np_arr_f> is a 2d numpy array,
     <coords_f> is a list containing the 2 coordinate variables corresponding to the dimensions 'time' and 'location' (which are hard-coded in the function)
     , varname_f is the variable name to be set in the xarray DataArray'''
    
    xr_arr_f = xr.DataArray(np_arr_f, coords=coords_f, dims=['time', 'location'], name=varname_f)
    xr_arr_f.attrs['standard_name'] = varname_f
    xr_arr_f.attrs['long_name'] = longname_f
    xr_arr_f.attrs['units'] = units_f
    xr_arr_f.location.attrs['standard_name'] = 'location index'
    xr_arr_f.location.attrs['long_name'] = 'index of the station location'
    xr_arr_f.location.attrs['altitude'] = altitude_f
    xr_arr_f.location.attrs['station_name'] = station_name_f
    xr_arr_f.location.attrs['aemet_code'] = aemet_code_f
    xr_arr_f.location.attrs['latitude'] = latitude_f
    xr_arr_f.location.attrs['longitude'] = longitude_f
    xr_arr_f.attrs['source'] = 'This netCDF file contains daily solar radiation time series from AEMET provided by Santiago Begueria in the framework of the PTIclima project; the original csv file contained the 3 radation variables RDIRDIA, RDIFDIA and RGLODIA.'
    xr_arr_f.attrs['nan_percentage'] = np.isnan(np_arr_f).sum(axis=0)/np_arr_f.shape[0]*100
    return(xr_arr_f)
