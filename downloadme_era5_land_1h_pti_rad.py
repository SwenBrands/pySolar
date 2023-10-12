import cdsapi
import numpy as np
import os
c = cdsapi.Client()

startyears = list(range(2017,2023,1))
#startyears = [2017]
months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
endyears = np.copy(startyears) #currently not in use hereafter
varalias = "ssrd" #t2m,rsds and psl
region = 'Spain'
tardir = '/media/swen/ext_disk2/datos/OBSdata/era5_land/'+region+'/hour/'+varalias
fileformat = 'netcdf' #netcdf or netcdf.zip
#tardir = '/home/swen/datos/OBSData/era5_land/1h'

## EXECUTE #############################################################
if fileformat == 'netcdf':
    file_ending = 'nc'
elif fileformat == 'netcdf.zip':
    file_ending = 'zip'
else:
    raise Exception ('ERROR: unknown entry for <fileformat> input parameter !')

if varalias == 'psl':
    variable = "mean_sea_level_pressure"
elif varalias == 't2m':
    variable = '2m_temperature'
elif varalias == 'ssrd':
    variable = 'surface_solar_radiation_downwards'
else:
    raise Exception('ERROR: unknown entry for <varalias>')

if region == 'Spain':
    domain = [44, -10, 34, 6,]
elif region == 'Canaries':
    domain = [29.6, -18.5, 27.3, -13,]
else:
    raise Exception('ERROR: unknown entry for <region>')

for yy in np.arange(len(startyears)):    
    #create target directory if missing
    if os.path.isdir(tardir+'/'+str(startyears[yy])) != True:
        os.makedirs(tardir+'/'+str(startyears[yy]))#create target directory if missing

    for mm in np.arange(len(months)):
        print('INFO: Requesting ERA5-Land data for '+region+', year '+str(startyears[yy])+' and month '+months[mm])
        c.retrieve(
            'reanalysis-era5-land',
            {
                'variable': variable,
                'year': str(startyears[yy]),
                'month': months[mm],
                'day': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                    '13', '14', '15',
                    '16', '17', '18',
                    '19', '20', '21',
                    '22', '23', '24',
                    '25', '26', '27',
                    '28', '29', '30',
                    '31',
                ],
                'time': [
                    '00:00', '01:00', '02:00',
                    '03:00', '04:00', '05:00',
                    '06:00', '07:00', '08:00',
                    '09:00', '10:00', '11:00',
                    '12:00', '13:00', '14:00',
                    '15:00', '16:00', '17:00',
                    '18:00', '19:00', '20:00',
                    '21:00', '22:00', '23:00',
                ],
                'area': domain,
                'format': fileformat, #netcdf.zip is recommended
            },
            tardir+'/'+str(startyears[yy])+'/era5_land_1h_'+varalias+'_'+region+'_'+str(startyears[yy])+str(months[mm])+'.'+file_ending)

print('INFO: downloadme_era5_land_1h_pti_rad.py has been run successfully!')
