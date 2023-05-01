import datetime
import numpy as np
import xarray as xr
from netCDF4 import Dataset as netcdf_dataset

#These need to be installed.
import sys
sys.path.append('M:/job_access/python/dataset_assembly/era5/')
from era5_requests import era5_hourly_single_level_request
sys.path.append('M:/job_access/python/dataset_assembly/access_io/')
from access_output import get_access_output_filename,append_var_to_daily_tb_netcdf
sys.path.append('M:/job_access/python/dataset_assembly/util/')
from access_interpolators import time_interpolate_synoptic_maps_ACCESS


def add_ERA5_single_level_variable_to_ACCESS_output(*,year:int,
                                                      month:int,
                                                      day:int,
                                                      variable:tuple,
                                                      satellite:str,
                                                      dataroot:str,
                                                      verbose:bool=False):

    # Get the maps of observation times from the existing output file that already contains
    # times and Tbs

    filename = get_access_output_filename(year=year,month=month,day=day,satellite=satellite,dataroot=dataroot)

    try: 
        root_grp = netcdf_dataset(filename, 'r', format='NETCDF4')
    except:
        raise FileNotFoundError(f'{filename} not found')
    try:
        times = root_grp.variables['second_since_midnight'][:,:,:].filled(fill_value=-999)
    except:
        raise ValueError(f'Error finding "second_since_midnight" in {filename}')
    finally:
        root_grp.close()

    # Download ERA5 data from ECMWF for all 24 hours of day, and the first hour of the next day.

    date_to_load = datetime.datetime(year,month,day)
    next_day = date_to_load + datetime.timedelta(hours=24)

    try:
        file1 = era5_hourly_single_level_request(year=date_to_load.year, 
                                    month=date_to_load.month,
                                    day=date_to_load.day,
                                    variable=variable[0],
                                    target_path='L:/access/_temp/',
                                    full_day=True)
        file2 = era5_hourly_single_level_request(year=next_day.year, 
                                    month=next_day.month,
                                    day=next_day.day,
                                    variable=variable[0],
                                    target_path='L:/access/_temp/',
                                    full_day=False)
    except:
        raise RuntimeError('Problem downloading ERA5 data using cdsapi')

    #open the files, and combine the two files into a 25-map array
    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)
    skt_first_day = ds1[variable[1]].values
    skt_next_day  = ds2[variable[1]].values
    skt = np.concatenate((skt_first_day,skt_next_day),axis=0)

    #ERA-5 files are upside down relative to RSS convention.
    for i in range(0,25):
        skt[i,:,:] = np.flipud(skt[i,:,:])

    #interpolate the array of skt maps to the times in the "times" maps
    if verbose:
        print('Interpolating...')

    #list of times, each hour.
    skt_times = np.arange(0.0,86401.0,3600.0)

    #create output array
    skt_by_hour = np.zeros_like(times)*np.nan

    for hour_index in range(0,24):
        time_map = times[:,:,hour_index]
        skt_at_time_map = time_interpolate_synoptic_maps_ACCESS(skt,skt_times,time_map)
        skt_by_hour[:,:,hour_index] = skt_at_time_map

    #write the results to the existing ouotput file
    append_var_to_daily_tb_netcdf(  year=year,
                                    month=month,
                                    day=day,
                                    satellite=satellite,
                                    var=skt_by_hour,
                                    var_name='skt',
                                    standard_name='skin_temperature',
                                    long_name='skin temperature interpolated from ERA5',
                                    valid_min=150.0,
                                    valid_max=400.0,
                                    units='degrees kelvin',
                                    v_fill = -999.0,
                                    dataroot='L:/access/',
                                    overwrite=True)





if __name__ =="__main__":
    
    year = 2012
    month = 7
    day = 11
    variable = ('Skin temperature','skt') # need this because var name for the ERA5 request is
                                        # not that same as the variable name in the nc file 
                                        # that is provided/downloaded
    satellite = 'AMSR2'
    verbose = True
    dataroot = 'L:/access/'

    add_ERA5_single_level_variable_to_ACCESS_output(year=year,
                                                    month=month,
                                                    day=day,
                                                    variable=variable,
                                                    satellite=satellite,
                                                    dataroot=dataroot,
                                                    verbose=True)
