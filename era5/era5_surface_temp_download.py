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



def era5_hourly_skin_temp_request(year, month, target_path):
    import cdsapi
    c = cdsapi.Client()
    from calendar import monthrange

    month_rng = monthrange(year, month)
    num_days_in_month = month_rng[1]
    day_list = []
    for day in range(1, num_days_in_month + 1):
        day_list.append(str(day).zfill(2))

    target = target_path + "ERA5_Skin_Temperature_%04d_%02d.nc" % (year, month)
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
                'year': str(year).zfill(4),
                'month': str(month).zfill(2),
                'day': day_list,
                'time': [
                    '00:00',
                    '03:00',
                    '06:00',
                    '09:00',
                    '12:00',
                    '15:00',
                    '18:00',
                    '21:00'
                ]
            }, temp_file)
        os.rename(temp_file,target)






yearStart = 2016
yearEnd = 2020
monthStart = 1
monthEnd = 12
for year in list(range(yearStart, yearEnd + 1)):
    nc_path = f"L:/access/era5/skin_temperature/{year:04d}/"
    os.makedirs(nc_path,exist_ok = True)

    for month in list(range(monthStart, monthEnd + 1)):
        era5_hourly_skin_temp_request(year, month, nc_path)
