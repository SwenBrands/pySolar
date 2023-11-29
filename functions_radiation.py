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
    
    #generate domain attribute to assign whether the station lies in Iberia or on the Canary Islands
    domain_f = np.repeat('IB',len(latitude_f))
    canaries_ind = np.where((latitude_f < 30) & (longitude_f < -10))[0].tolist()
    domain_f[canaries_ind] = 'CA'
    xr_arr_f.location.attrs['domain'] = domain_f
    xr_arr_f.location.attrs['info'] = 'IB refers to stations located in mainland Iberia, Baleares and northern Africa. CA refers to stations located on the Canary Islands. Station altitude is in metres aboves sea level.'

    xr_arr_f.attrs['source'] = source_f
    xr_arr_f.attrs['nan_percentage'] = np.isnan(np_arr_f).sum(axis=0)/np_arr_f.shape[0]*100
    return(xr_arr_f)


def get_nc_path(years_f,region_f,dir_f,timescale_f,variable_f):
    '''Loads ERA5-Land data from local diks. Input: years_f = a numpy array containing a list of years to be loaded; region_f = string with a label referring to the geographical area to be loaded;
    dir_f = string, basic path to the files which will be extended by the function; timescale_f = string, temporal aggregation of the files; variable_f = string, name of the variable as provided in the input nc files. Output: get_nc_path = list of 
    nc files to be loaded. Currently only works for hourly aggregation, so an <aggregation> parameter should be included in future versions'''
    #get list of full paths of nc dataset to be tested; sp = Iberia, ca = Canarias
    listdir_f = []
    if len(years_f.shape) == 0: #in case years_f is an interger containing a single year
        root_year_f = dir_f+'/'+region_f+'/'+timescale_f+'/'+variable_f+'/'+str(years_f)
        listdir_year = os.listdir(root_year_f)
        listdir_year_full = [root_year_f+'/'+listdir_year[ii] for ii in np.arange(len(listdir_year))] #contains the full path to the files
        listdir_f = np.append(listdir_f,listdir_year_full,axis=0)
    else: #in case years_f is an numpy.ndarray (test also with a list in the future) containing two or more years
        for yy_f in np.arange(len(years_f)):
            root_year_f = dir_f+'/'+region_f+'/'+timescale_f+'/'+variable_f+'/'+str(years_f[yy_f])
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

def plot_pcolormesh(xr_ds_f,score_f,minval_f,maxval_f,savename_f,colormap_f,dpival_f):
    '''Plots matrix of the verfication results contained in xarray dataset <xr_ds_f>, indicated by the string <score_f>. Seasons are plotted on the x axis, stations on the y axis.'''

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
       print('Info: There is a problem with the alpha parameter when generating the figure on my local system. Correct this in future versions !')
    plt.savefig(savename_f,dpi=dpival_f)
    plt.close('all')

def disaggregate_rean(xr_ds_f,variable_f,accumulation_f):
    '''disaggregates ERA5-Land data from accumulations from 0 UTC to hour indicated in the file to hour-to-hour accumulations (i.e. from 23 to 0, 0 to 1, 1 to 2 UTC etc.)'''
    np_arr_f = np.zeros(xr_ds_f[variable_f].shape) #get numpy array with dimensions of xr data array containing the accumulated rean. data; will be filled below
    np_arr_f[:] = np.nan
    dates_f = pd.DatetimeIndex(xr_ds_f.time.values)
    hours_unique_f = np.unique(dates_f.hour.values) #get unique hours
    xr_arr_f = xr_ds_f[variable_f]
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
        ##set negative values to zero
        neg_mask_f = np_arr_f < 0 
        np_arr_f[neg_mask_f] = 0.
        negvals_f = np.any(neg_mask_f)
        del(neg_mask_f)
    
    #nc_ca_orig = nc_ca.copy(deep=True) # make a copy of the original aggegeated xr data array
    xr_ds_f[variable_f][:] = np_arr_f #replace values in aggregated xr dataset with disaggregated values
    
    #define the accumulation method
    if accumulation_f in ('forward','centred','centered'):
        xr_ds_f[variable_f] = xr_ds_f[variable_f].shift(time=-1)
        time_unit_f = 'hours since '+str(xr_ds_f.time[0].values) #set the time unit to be encoded
        xr_ds_f.time.encoding['units'] = time_unit_f
        if accumulation_f in ('centred','centered'):
            print('As requested by the user, the disaggregated hourly data is '+accumulation_f+' around the time instant indicated in the time dimension (i.e. HH:30) and time bounds are added.')
            xr_ds_f = center_dates(xr_ds_f,'hour') #hour is hardcoded here because disaggregate_rean() only works with hourly data, i.e. daily data would not be passed in the future. However, try to find a more robust solution in future versions of the function.
        else:
            print('As requested by the user, the data is accumulated '+accumulation_f+', i.e. the time instant in the output data array indicates the start of the accumulation period.')             
    elif accumulation_f == 'backward':
        time_unit_f = 'hours since '+str(xr_ds_f.time[0].values) #set the time unit to be encoded
        xr_ds_f.time.encoding['units'] = time_unit_f
        print('As requested by the user, the data is accumulated '+accumulation_f+', i.e. the time instant in the output data array indicates the end of the accumulation period.')
    else:
        raise Exception('ERROR: check entry for <accumulation_f> !')
    return(xr_ds_f,negvals_f)

def clean_directory_content(directory_f):
    '''Cleans the content of the specified directory <directory_f>, which is the only input parameter, the full path to the directory is needed. Output: none except print messages'''
    listdir_f = os.listdir(directory_f)
    if len(listdir_f) > 0:
        print('The following previously generated files in '+directory_f+' will be removed:')
        print(listdir_f)
        print(' ')
        [os.remove(listdir_f[ffi]) for ffi in np.arange(len(listdir_f))]
    else:
        print('The user requested '+directory_f+' to be cleaned but this folder is empty and is thus already clean.')
        print(' ')

def print_array_stats(xr_arr_f,text_f):
    '''prints the nan percentage, min, max and mean values of the input xarray array <xr_arr_f>, preceded by a short text description contained in <text_f> to know what these arrays are about.'''
    print(text_f)
    print('NaN percentage: '+str(np.round(np.sum(np.isnan(xr_arr_f)).values/xr_arr_f.size*100,4)))
    print('Minimum: '+str(np.round(xr_arr_f.min().values,2)))
    print('Mean :'+str(np.round(xr_arr_f.mean().values,2)))
    print('Maximum '+str(np.round(xr_arr_f.max().values,2)))

def xr_ds_to_netcdf(xr_ds_hour_f,xr_ds_day_f,encoding_f,file_style_f,dir_hour_f,dir_day_f,file_hour_f,file_day_f):
    '''saves hourly and daily data stored in <xr_ds_hour_f> and <xr_ds_day_f> to compressed netCDF files in year-to-year or month-to-month files (stored in yearly directories).'''
    if file_style_f == 'yearly': #save to year-to-year files
        #save hourly data
        start_date_hour_f = str(xr_ds_hour_f.time.values[0]).replace('-','').replace(':','')[0:-12]
        end_date_hour_f = str(xr_ds_hour_f.time.values[-1]).replace('-','').replace(':','')[0:-12]
        savename_hour_f = dir_hour+'/'+file_hour+'_'+start_date_hour_f+'_'+end_date_hour_f+'.nc'
        print('saving '+savename_hour_f)
        xr_ds_hour_f.to_netcdf(savename_hour_f,encoding=encoding_f)
        
        #save daily data
        start_date_day_f = str(xr_ds_day_f.time.values[0]).replace('-','').replace(':','')[0:-12]
        end_date_day_f = str(xr_ds_day_f.time.values[-1]).replace('-','').replace(':','')[0:-12]
        savename_day_f = dir_day+'/'+file_day+'_'+start_date_day_f+'_'+end_date_day_f+'.nc'
        print('saving '+savename_day_f)
        xr_ds_day_f.to_netcdf(savename_day_f,encoding=encoding_f)

    elif file_style_f == 'monthly': #save month-to-month files in yearly directories
        #define paths of the directories where the output netCDF data will be saved
        dates_hour_f = pd.DatetimeIndex(xr_ds_hour_f.time)
        dates_day_f = pd.DatetimeIndex(xr_ds_day_f.time)
        year_f = np.unique(dates_day_f.year)
        
        #check whether the xarray dataset <xr_ds_day> contains only a single year as expected
        if len(np.unique(np.array(year_f))) > 1:
            raise Exception('ERROR: the xarray dataset <xr_ds_day> contains values for more than one year, but a single year is expected here !')

        #save one netCDF file per month in directory named <year_f>
        months_f = np.unique(dates_day_f.month)
        for mo in np.arange(len(months_f)):
            print('Subsetting output netCDF data for month '+str(months_f[mo])+' and year '+str(year_f)+'...')
            #save hourly data
            monthind_f = np.where(dates_hour_f.month == months_f[mo])[0]
            nc_month_hour_f = xr_ds_hour_f.isel(time=monthind_f)
            start_date_hour_f = str(nc_month_hour_f.time.values[0]).replace('-','').replace(':','')[0:-12]
            end_date_hour_f = str(nc_month_hour_f.time.values[-1]).replace('-','').replace(':','')[0:-12]
            savename_hour_f = dir_hour+'/'+file_hour+'_'+start_date_hour_f+'_'+end_date_hour_f+'.nc'
            nc_month_hour_f.to_netcdf(savename_hour_f,encoding=encoding_f)
            nc_month_hour_f.close()
            del(nc_month_hour_f)
                
            #save daily data
            monthind_f = np.where(dates_day_f.month == months_f[mo])[0]
            nc_month_day_f = xr_ds_day_f.isel(time=monthind_f)
            start_date_day_f = str(nc_month_day_f.time.values[0]).replace('-','').replace(':','')[0:-12]
            end_date_day_f = str(nc_month_day_f.time.values[-1]).replace('-','').replace(':','')[0:-12]
            savename_day_f = dir_day+'/'+file_day+'_'+start_date_day_f+'_'+end_date_day_f+'.nc'
            nc_month_day_f.to_netcdf(savename_day_f,encoding=encoding_f)
            nc_month_day_f.close()
            del(nc_month_day_f)
    else:
        raise Exception('ERROR: unknown entry for <file_style_f> input parameter !')
    
    xr_ds_hour_f.close()
    xr_ds_day_f.close()
    del(xr_ds_hour_f,xr_ds_day_f)


def get_temporal_aggregation_metadata(variable_f,variable_unit_f,accumulation_f):
    '''get metadata for a given variable, unit and temporal aggregation defined by <variable_f>, <variable_unit_f> and <accumulation_f>, all being character string'''
    if variable_f == 'ssrd':
        if accumulation_f == 'forward':
            meta_hour_f = 'Accumulated flux data per second (in '+variable_unit_f+') assumed to be constant from hour h to h+1 (e.g. from 00 to 01 UTC of a given day), with h being indicated in the <time> dimension.'
            meta_day_f = 'Accumulated flux data per second (in '+variable_unit_f+') assumed to be constant from hour 00 to 24 UTC of day d, with d being indicated in the <time> dimension.'
        elif accumulation_f == 'backward':
            meta_hour_f = 'Accumulated flux data per second (in '+variable_unit_f+') assumed to be constant from hour h-1 to h (e.g. from 00 to 01 UTC of a given day), with h being indicated in the <time> dimension.'
            meta_day_f = 'Accumulated flux data per second (in '+variable_unit_f+') assumed to be constant from hour 23 UTC of day d-1 to hour 23 UTC of day d, with d being indicated in the <time> dimension.'
        elif accumulation_f in ('centred','centered'):
            meta_hour_f = 'Accumulated flux data per second (in '+variable_unit_f+') assumed to be constant from hour HH to HH+1 (e.g. from 00 to 01 UTC of a given day), with HH:30:00 UTC (e.g. 00:30:00) indicated in the <time> dimension.'
            meta_day_f = 'Accumulated flux data per second (in '+variable_unit_f+') assumed to be constant from hour 00 to 24 UTC of day d, with 12:00:00 UTC indicated in the <time> dimension.'
        else:
            raise Exception('Error: Unknown entry for input parameter named <accumulation_f> !')
    elif variable_f == 'pvpot':
        if accumulation_f == 'forward':
            meta_hour_f = 'hourly index data calculated upon rsds data from h to h+1 with h being the time instant indicated in the time dimension pointing to the start of the accumulation period. Underlying u10, v10 and tas are instantaneous values sampled at h.'
            meta_day_f = 'daily mean of the hourly index data, valid for the time period 00 UTC to 24 UTC of day d, with d being indicated in the time dimension.'
        elif accumulation_f == 'backward':
            meta_hour_f = 'hourly index data calculated upon rsds data from h-1 to h with h being the time instant indicated in the time dimension pointing to the end of the accumulation period. Underlying u10, v10 and tas are instantaneous values sampled at h.'
            meta_day_f = 'daily mean of the hourly index data, valid for the time period 23 UTC of day d-1 to 23 UTC of day d, with d being indicated in the time dimension.'
        elif accumulation_f in ('centred','centered'):
            meta_hour_f = 'hourly index data assumed to be constant from hour HH to HH+1 (e.g. from 00 to 01 UTC of a given day), with HH:30:00 UTC (e.g. 00:30:00) indicated in the <time> dimension.'
            meta_day_f = 'daily mean of the hourly index data, accumulated from hour 00 to 24 UTC of day d, with 12:00:00 UTC indicated in the <time> dimension.'
        else:
            raise Exception('Error: Unknown entry for input parameter named <accumulation_f> !')
    elif variable_f == 'tp':
        if accumulation_f == 'forward':
            meta_hour_f = 'Accumulated precipitation amount in '+variable_unit_f+' from hour h to h+1 (e.g. from 00 to 01 UTC of a given day), with h being indicated in the <time> dimension.'
            meta_day_f = 'Accumulated precipitation amount in '+variable_unit_f+' from hour 00 to 24 UTC of day d, with d being indicated in the <time> dimension.'
        elif accumulation_f == 'backward':
            meta_hour_f = 'Accumulated precipitation amount in '+variable_unit_f+' from hour h-1 to h (e.g. from 00 to 01 UTC of a given day), with h being indicated in the <time> dimension.'
            meta_day_f = 'Accumulated precipitation amount in '+variable_unit_f+' from hour 23 UTC of day d-1 to hour 23 UTC of day d, with d being indicated in the <time> dimension.'
        elif accumulation_f in ('centred','centered'):
            meta_hour_f = 'Accumulated precipitation amount in '+variable_unit_f+' from hour HH to HH+1 (e.g. from 00 to 01 UTC of a given day), with HH:30:00 UTC (e.g. 00:30:00) indicated in the <time> dimension.'
            meta_day_f = 'Accumulated precipitation amount in '+variable_unit_f+' from hour 00 to 24 UTC of day d, with 12:00:00 UTC indicated in the <time> dimension.'
        else:
            raise Exception('Error: Unknown entry for input parameter named <accumulation_f> !')
    else:
        raise Exception('ERROR: unknown entry for input parameter <variable_f> !')
    return(meta_hour_f,meta_day_f)

def center_dates(xr_ds_f,timescale_f):
    ''' centers the dates along the time dimension in the input xarray dataset <xr_ds_f> to HH:30 for <timescale_f> = 'hourly' or to 12:00 for <timescale_f> = 'daily', respectively.'''
    dates_hour_f = pd.DatetimeIndex(xr_ds_f.time.values)
    if timescale_f == 'hour':
        dates_transformed_f = pd.DatetimeIndex([dates_hour_f[ii].replace(minute=30) for ii in np.arange(len(dates_hour_f))]) #set hour from HH:00 to HH:30 format to indicate the centre of the accumulation period
        time_unit_f = 'hours since '+str(dates_transformed_f[0]) #set the time units to be encoded below
    elif timescale_f == 'day':
        dates_transformed_f = pd.DatetimeIndex([dates_hour_f[ii].replace(hour=12) for ii in np.arange(len(dates_hour_f))]) #set hour 00 to hour 12 to indicate the centre of the accumulation period
        time_unit_f = 'days since '+str(dates_transformed_f[0]) #set the time units to be encoded below
    else:
        raise Excpetion('ERROR: unknown entry for <timescale_f> !')
    #time_unit_f = xr_ds_f.time.encoding['units'] #catch the time units passed to function via <xr_ds_f>
    xr_ds_f = xr_ds_f.assign_coords({"time":dates_transformed_f}) #reassign time dimension
    xr_ds_f = xr_ds_f.cf.add_bounds("time") # add time bounds https://cf-xarray.readthedocs.io/en/latest/generated/xarray.Dataset.cf.add_bounds.html

    #xr_ds_f.time.attrs['units'] = time_unit_f 
    xr_ds_f.time.attrs['standard_name'] = 'time'
    xr_ds_f.time.encoding['units'] = time_unit_f #set the time unit of the newly generated/transformed time dimension, https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/ch04s04.html
    #xr_ds_f.time_bounds.attrs['units'] = time_unit_f 
    xr_ds_f.time_bounds.attrs['standard_name'] = 'time'
    xr_ds_f.time_bounds.encoding['units'] = time_unit_f 
    return(xr_ds_f)
