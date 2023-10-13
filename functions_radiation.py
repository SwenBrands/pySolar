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
    
def add_season_metadata(xr_ds_f,months_f,months_labels_f):
    '''Adds season metadata to xarray dataset xr_ds_f containing various validation scores. The metadata is taken from the lists months_f and months_labels_f'''
    #assign season attributes
    xr_ds_f.season.attrs['standard_name'] = 'season index'  
    xr_ds_f.season.attrs['long_name'] = 'index of season'
    xr_ds_f.season.attrs['monhts'] = months
    xr_ds_f.season.attrs['season_label'] = months_labels
    return(xr_ds_f)

def plot_pcolormesh(xr_ds_f,score_f,savename_f,colormap_f,dpival_f):
    '''Plots matrix of the verfication results contained in xarray dataset <xr_ds_f>, indicated by the string <score_f>. Seasons are plotted on the x axis, stations on the y axis.'''
    #set min and max values as a function of the score to be plotted
    if score_f in ('pearson_r','spearman_r'):
        minval_f = 0
        maxval_f = 1
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
