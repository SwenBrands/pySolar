#!/usr/bin/env python

'''functions used by the scripts located in .../pticlima/radiation/python'''

def get_xr_arr(np_arr_f, coords_f, varname_f, longname_f, units_f, altitude_f, station_name_f, aemet_code_f, latitude_f, longitude_f, source_f):
    
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
    
    #generate domain attribute to assign whether the station lies in Spain or on the Canary Islands
    domain_f = np.repeat('IB',len(latitude_f))
    canaries_ind = np.where((latitude_f < 30) & (longitude_f < -10))[0].tolist()
    domain_f[canaries_ind] = 'CA'
    xr_arr_f.location.attrs['domain'] = domain_f
    xr_arr_f.location.attrs['info'] = 'IB refers to stations located in mainland Spain, Baleares and northern Africa. CA refers to stations located on the Canary Islands. Station altitude is in metres aboves sea level.'

    xr_arr_f.attrs['source'] = source_f
    xr_arr_f.attrs['nan_percentage'] = np.isnan(np_arr_f).sum(axis=0)/np_arr_f.shape[0]*100
    return(xr_arr_f)


def get_nc_path(years_f,region_f,dir_f,variable_f):
    '''Loads ERA5-Land data from local diks. Input: year_f = a numpy array containing a list of years to be loaded; region_f = string with a label referring to the geographical area to be loaded;
    dir_f = string, basic path to the files which will be extended by the function; variable_f = string, name of the variable as provided in the input nc files. Output: get_nc_path = list of 
    nc files to be loaded.'''
    #get list of full paths of nc dataset to be tested; sp = Spain, ca = Canaries
    listdir_f = []
    for yy in np.arange(len(years_f)):
        root_year_f = dir_f+'/'+region_f+'/hour/'+variable_f+'/'+str(years_f[yy])

        listdir_year = os.listdir(root_year_f)
        listdir_year_full = [root_year_f+'/'+listdir_year[ii] for ii in np.arange(len(listdir_year))] #contains the full path to the files
       
        listdir_f = np.append(listdir_f,listdir_year_full,axis=0)

    #remove all detected files that are not nc format from the lists to be loaded, just in case someone put other files into the directory
    dropind_f = []
    for ff in list(range(len(listdir_f))):
        if listdir_f[ff][-3:] != '.nc':
            dropind_f.append(ff)

    #delete non nc files from the lists and sort the individual filenames in ascending temporal order
    listdir_f = sorted(np.delete(listdir_f,dropind_f).tolist())
    return(listdir_f)

def add_location_metadata(xr_ds_f,xr_arr_obs_f,xr_arr_mod_f):
    '''Adds location metadata to xarray dataset xr_ds_f containing various validation scores. The metadata is taken from the xr DataArrays xr_arr_obs_f and xr_arr_mod_f'''
    #assgin location attibutes
    xr_ds_f.location.attrs['standard_name'] = 'location index'
    xr_ds_f.location.attrs['long_name'] = 'index of the station location'
    xr_ds_f.location.attrs['altitude_obs'] = xr_arr_obs_f.location.altitude
    xr_ds_f.location.attrs['station_name'] = xr_arr_obs_f.location.station_name
    xr_ds_f.location.attrs['aemet_code'] = xr_arr_obs_f.location.aemet_code
    xr_ds_f.location.attrs['latitude_obs'] = xr_arr_obs_f.location.latitude
    xr_ds_f.location.attrs['longitude_obs'] = xr_arr_obs_f.location.longitude
    xr_ds_f.location.attrs['latitude_nn'] = xr_arr_mod_f.location.latitude
    xr_ds_f.location.attrs['longitude_nn'] = xr_arr_mod_f.location.longitude
    xr_ds_f.location.attrs['altidue_nn'] = xr_arr_mod_f.location.altitude
    xr_ds_f.location.attrs['info'] = 'obs and nn refer to observations and nearest neighbour model or reanalysis data, respectively.'
    return(xr_ds_f)
    xr_ds_f.close()
    del xr_ds_f
    
def add_season_metadata(xr_ds_f,months_f,months_labels_f):
    '''Adds season metadata to xarray dataset xr_ds_f containing various validation scores. The metadata is taken from the lists months_f and months_labels_f'''
    #assign season attributes
    xr_ds_f.season.attrs['standard_name'] = 'season index'  
    xr_ds_f.season.attrs['long_name'] = 'index of season'
    xr_ds_f.season.attrs['season_label'] = months_labels
    #xr_ds_f.season.attrs['months'] = months #multi-dimensional array attributes are not supported and <months> has 2 dimensions. This is why this line is commented.
    return(xr_ds_f)

def plot_pcolormesh(xr_ds_f,score_f,savename_f,colormap_f,dpival_f):
    '''Plots matrix of the verfication results contained in xarray dataset <xr_ds_f>, indicated by the string <score_f>. Seasons are plotted on the x axis, stations on the y axis.'''
    #set min and max values as a function of the score to be plotted
    if score_f in ('pearson_r','spearman_r'):
        minval_f = 0
        maxval_f = 1
    elif score_f in ('bias','relbias'):
        maxval_f = np.abs(xr_ds_f[score_f]).max().values
        minval_f = maxval_f*-1.
    elif score_f in ('rmse','mae','mape'):
        maxval_f = xr_ds_f[score_f].max().values
        minval_f = xr_ds_f[score_f].min().values
    else:
        raise Excpetion('ERROR: check entry of <score_f> input parameter in the function plot_pcolormesh() !')
    
    fig = plt.figure()
    ax = xr_ds_f[score_f].plot.pcolormesh(cmap = colormap_f, x = 'season', y = 'location', vmin = minval_f, vmax = maxval_f, add_colorbar=False)
    ax.axes.set_yticks(xr_ds_f.location.values)
    ax.axes.set_yticklabels(xr_ds_f.location.station_name,fontsize=4)
    ax.axes.set_xticks(xr_ds_f.season.values)
    ax.axes.set_xticklabels(xr_ds_f.season.season_label,fontsize=2, rotation = 45.)
    ax.axes.set_aspect('auto')
    plt.xticks(fontsize=5)
    plt.xlabel(None)
    plt.ylabel(None)
    cbar = plt.colorbar(ax,shrink=0.5,label=xr_ds_f[score_f].name + ' ('+xr_ds_f[score_f].units+')', orientation = 'horizontal')
    cbar.ax.tick_params(labelsize=6)
    fig.tight_layout()
    if figformat == 'pdf': #needed to account for irregular behaviour with the alpha parameter when plotting a pdf file
       #fig.set_rasterized(True)
       print('Info: There is a problem with the aplha parameter when generating the figure on my local system. Correct this in future versions !')
    plt.savefig(savename_f,dpi=dpival_f)
    plt.close('all')

def disaggregate_rean(xr_ds_f,variable_f):
    #disaggregate ERA5-Land data from accumulations from 0 UTC to hour indicated in the file to hour-to-hour accumulations (i.e. from 23 to 0, 0 to 1, 1 to 2 UTC etc.)
    np_arr_f = np.zeros(xr_ds_f[variable_f].shape) #get numpy array with dimensions of xr data array containing the accumulated rean. data; will be filled below
    np_arr_f[:] = np.nan
    dates_f = pd.DatetimeIndex(xr_ds_f.time.values)
    hours_unique_f = np.unique(dates_f.hour.values) #get unique hours
    xr_arr_f = xr_ds_f[variable]
    for hh in np.arange(len(hours_unique_f)): #loop through every hour
        print('INFO: disaggregating hourly data for '+str(hours_unique_f[hh])+' UTC...')
        hourind_f = np.where(dates_f.hour.values == hours_unique_f[hh])[0] #find index for the hour in the non-lagged DatetimeIndex
        hourind_lag1_f = hourind_f-1 #define lagged indices
        if hours_unique_f[hh] == 0:
            hourind_f = np.delete(hourind_f,0)
            hourind_lag1_f = np.delete(hourind_lag1_f,0)
            fillval_f = xr_arr_f.isel(time=hourind_f).values - xr_arr_f.isel(time=hourind_lag1_f).values
        elif hours_unique_f[hh] == 1:
            fillval_f = xr_arr_f.isel(time=hourind_f).values
        else:
            fillval_f = xr_arr_f.isel(time=hourind_f).values - xr_arr_f.isel(time=hourind_lag1_f).values
        np_arr_f[hourind_f,:,:] = fillval_f
        np_arr_f[np_arr_f < 0 ] = 0.
    
    #nc_ca_orig = nc_ca.copy(deep=True) # make a copy of the original aggegeated xr data array
    xr_ds_f[variable_f][:] = np_arr_f #replace values in aggreget xr data array with disaggregated values
    return(xr_ds_f)
    
    # xr_ds_da_f = xr_ds_f.copy() #da for disaggregated
    # xr_ds_da_f[variable_f][:] = np.nan #set complete data array within to nan
    # dates_f = pd.DatetimeIndex(xr_ds_f.time.values)
    # hours_unique_f = np.unique(dates_f.hour.values) #get unique hours
    # hours_unique_shifted_f = [hours_unique_f[-1]]+list(hours_unique_f)[0:-1] #get unique hours shifted by t-1 to get hourly accumulations from t-1 to 1
    # for hh in np.arange(len(hours_unique_f)): #loop through every hour
        # print('INFO: disaggregating hourly data ending at '+str(hours_unique_f[hh])+' UTC...')
        # hourind_f = np.where(dates_f.hour.values == hours_unique_f[hh])[0] #find index for the hour in non-lagged DatetimeIndex
        # hourind_f = np.delete(hourind_f,0) #remove first entry which is not paired in hourind_lag1
        # hourind_lag1_f = np.where(dates_f.hour.values == hours_unique_shifted_f[hh])[0] #find index for the hour in non-lagged DatetimeIndex
        # hourind_lag1_f = np.delete(hourind_lag1_f,-1) #remove last entry which is not paired in hourind
        # xr_ds_da_f[variable_f][hourind_f] = xr_ds_f[variable_f].isel(time=hourind_f).values - xr_ds_f[variable_f].isel(time=hourind_lag1_f).values
    # return(xr_ds_da_f)
