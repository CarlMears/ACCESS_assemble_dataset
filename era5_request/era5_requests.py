#!/usr/bin/env python3
"""
This a collection of subroutines that download one day for ERA5 data (including the first hour of the
next day).
In order to use it, the Copernicus Data Server (cdsapi) needs to be installed in your version of python.

see https://cds.climate.copernicus.eu/api-how-to

If you use conda:
conda install -c conda-forge cdsapi

Then you need to get a key from ECMWF and put in a .cdsapirc file in your {user} folder

"""

import calendar
import os.path
import os, errno
import datetime



def era5_hourly_single_level_request(*,year, month, day, variable, target_path,full_day=True):
    import cdsapi
    c = cdsapi.Client()
    from calendar import monthrange



    target = target_path + "ERA5_Skin_Temperature_%04d_%02d.nc" % (year, month)

    if full_day:
        target = f'{target_path}ERA5_Skin_Temperature_{year:04d}_{month:02d}_{day:02d}.full.nc'
        times = [f"{h:02d}:00" for h in range(0,24)]
    else:
        target = f'{target_path}ERA5_Skin_Temperature_{year:04d}_{month:02d}_{day:02d}.1st_hour.nc'
        times = '00:00'

    temp_file =  target_path + "temp.nc"

    file_already_exists = os.path.isfile(target)

    if file_already_exists:
        print('File: ' + target + ' already exists, skipping')
    else:
        print('Getting: '+target)
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
				'grid':'0.25/0.25',
                'variable': 'Skin temperature',
                'year': f"{year:04d}",
                'month': f"{month:02d}",
                'day': f"{day:02d}",
                'time': times
            }, temp_file)
        os.rename(temp_file,target)

    return target

def era5_hourly_pressure_level_request(*,year, 
                                         month, 
                                         day, 
                                         variable, 
                                         levels = ['775','875','975'],
                                         target_path,
                                         full_day=True):
    import cdsapi
    c = cdsapi.Client()
    from calendar import monthrange


    if full_day:
        target = f'{target_path}ERA5_{variable}_{year:04d}_{month:02d}_{day:02d}.full.nc'
        times = [f"{h:02d}:00" for h in range(0,24)]
    else:
        target = f'{target_path}ERA5_{variable}_{year:04d}_{month:02d}_{day:02d}.1st_hour.nc'
        times = '00:00'

    temp_file =  target_path + "temp.nc"

    file_already_exists = os.path.isfile(target)

    if file_already_exists:
        print('File: ' + target + ' already exists, skipping')
    else:
        print('Getting: '+target)
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
				'grid':'0.25/0.25',
                'variable': variable,
                'pressure_level': levels,
                'year': f"{year:04d}",
                'month': f"{month:02d}",
                'day': f"{day:02d}",
                'time': times
            }, temp_file)
        os.rename(temp_file,target)

    return target


if __name__ == "__main__":

    year = 2012
    month = 7
    day = 11
    variable = 'temperature'
    target_path = "L:/access/_temp/"

    levels = ['1','2','3','5',
              '7','10','20','30',
              '50','70','100','125',
              '150','175','200','225',
              '250',' 300','350','400',
              '450','500','550','600',
              '650','700','750','775',
              '800','825','850','875',
              '900','925','950','975','1000']

    file = era5_hourly_pressure_level_request(year = year, 
                                         month = month, 
                                         day = day, 
                                         variable=variable, 
                                         levels = levels,
                                         target_path=target_path,
                                         full_day=True)

    print(file)
    print







