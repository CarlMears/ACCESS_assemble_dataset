import numpy as np
import xarray as xr


import matplotlib.pyplot as plt

from global_map import global_map  #must be installed from rss_plotting package

from access_io.access_output import get_access_output_filename,write_daily_tb_netcdf
from resampled_tbs.read_resampled_orbit import get_orbit_range,read_resampled_tbs
from util.numpy_date_utils import calendar_dates_from_datetime64,convert_to_sec_in_day,convert_to_np_datetime64
from util.orbit_times_amsr2 import find_orbits_in_day,read_amsr2_orbit_times

def make_daily_ACCESS_tb_file(*,year,month,day,satellite,dataroot,channel_list,verbose=False,plot_example_map=True):
    
    if satellite == 'amsr2':
        orbit_times = read_amsr2_orbit_times()
    else:
        raise ValueError(f'Orbit Times for {satellite} no omplemented yet')
   
    num_lats = 721
    num_lons = 1440
    num_hours = 24
    num_channels = 14 #all possible AMSR2 channels 
    list_of_channels=['time','6V','6H','7V','7H','11V','11H','19V','19H',
                     '24V','24H','37V','37H','89V','89H']

    
    #initialize arrays for daily data

    at_least_one_orbit = False
    tb_array_by_hour = np.zeros((721,1440,24,num_channels),dtype=np.float32)*np.nan
    time_array_by_hour = np.zeros((721,1440,24))*np.nan
    
    file_list = []
    orbits_to_do = find_orbits_in_day(times_np64=orbit_times,year=year,month=month,day=day)
    print(f'Processing {year}/{month:02d}/{day:02d}, orbit: ',end = '')
    for orbit in orbits_to_do:
        print(f'{orbit} ',end='')
        #get the times for this orbit, and convert to time within the day...
        try:
            ob_time,filename = read_resampled_tbs(satellite=satellite,channel='time',orbit=orbit,verbose=False)
        except FileNotFoundError:
            print(f'No time file found for orbit: {orbit}')
            continue
        
        file_list.append(filename)
        obtime_in_day = convert_to_sec_in_day(ob_time,year,month,day)
        for hour in range(0,24):         
            time_slice =  time_array_by_hour[:,:,hour]
            start_time_sec = hour*3600.0
            end_time_sec = start_time_sec+3600.0
            ok = np.all([(obtime_in_day > start_time_sec),(obtime_in_day <= end_time_sec)],axis=0)
            num_ok = np.sum(ok)
            if num_ok > 0:
                time_slice[ok] = obtime_in_day[ok]

        if verbose:
            print('reading resampled tb files')
        for channel in channel_list:
            try:
                tbs,filename = read_resampled_tbs(satellite=satellite,channel=channel,orbit=orbit)
            except FileNotFoundError:
                print(f'No file found for orbit: {orbit}, channel: {channel}')
                continue
            file_list.append(filename)
            at_least_one_orbit = True
            #loop through the hours, saving the tbs and times in hour-long slices
            for hour in range(0,24):
                tb_slice = tb_array_by_hour[:,:,hour,channel-1]
                start_time_sec = hour*3600.0
                end_time_sec = start_time_sec+3600.0
                ok = np.all([(obtime_in_day > start_time_sec),(obtime_in_day <= end_time_sec)],axis=0)
                num_ok = np.sum(ok)
                if num_ok > 0:
                    tb_slice[ok] = tbs[ok]
                    if verbose:
                        print(f'orbit:{orbit}, channel:{list_of_channels[channel]}, time_range({start_time_sec}:{end_time_sec}) Num Obs: {num_ok}')
    print() 
    if at_least_one_orbit:
        if plot_example_map:
            global_map(tb_array_by_hour[:,:,0,5],vmin=0.0,vmax=330)
    
        write_daily_tb_netcdf(year=year,month=month,day=day,satellite=satellite,
                    tb_array_by_hour=tb_array_by_hour,
                    time_array_by_hour=time_array_by_hour,
                    file_list=file_list,dataroot=dataroot)

    return file_list

if __name__ == "__main__":

    year = 2012
    month = 7
    channel_list = range(5,13)
    satellite='amsr2'
    dataroot = 'L:/access/'

    for day in range(13,17):
        make_daily_ACCESS_tb_file(  year=year,
                                    month=month,
                                    day=day,
                                    satellite=satellite,
                                    dataroot=dataroot,
                                    channel_list=channel_list,
                                    verbose=False,
                                    plot_example_map=True)

    plt.show()
    print()

