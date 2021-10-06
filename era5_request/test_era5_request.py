# -*- coding: utf-8 -*-
"""
Created on Mon Dec 04 09:18:09 2018

This script retrieves 10m neutral stability  and regular U amd V winds from the ERA5 reanalysis and stores them on N:/data/.

In order to use it, the New Copernicus Data Server needs to be installed in your version of python.

see https://cds.climate.copernicus.eu/api-how-to

@author: mears
"""

# !/usr/bin/env python
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






year = 2012
month = 7
day = 11
variable = 'Skin temperature'
date_to_load = datetime.datetime(year,month,day)
next_day = date_to_load + datetime.timedelta(hours=24)

file1 = era5_hourly_single_level_request(year=date_to_load.year, 
                              month=date_to_load.month,
                              day=date_to_load.day,
                              variable=variable,
                              target_path='L:/access/_temp/',
                              full_day=True)
file2 = era5_hourly_single_level_request(year=next_day.year, 
                              month=next_day.month,
                              day=next_day.day,
                              variable=variable,
                              target_path='L:/access/_temp/',
                              full_day=False)

print(file1,file2)
print