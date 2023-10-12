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
    '''loads ERA5-Land data from local diks. Input: year_f = a numpy array containing a list of years to be loaded; region_f = string with a label referring to the geographical area to be loaded;
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
