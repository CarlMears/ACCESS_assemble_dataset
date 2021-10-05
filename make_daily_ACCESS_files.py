if __name__ == "__main__":

    import numpy as np
    import xarray as xr

    from global_map import global_map
    import matplotlib.pyplot as plt

    import sys
    sys.path.append('C:/Users/Carl/python_packages_cam/')
    from land_masks.land_fraction import read_land_fraction_1440_720

    sys.path.append('./access_io/')
    sys.path.append('./resampled_tbs/')
    sys.path.append('./util/')

    from access_output import get_access_output_filename,append_lf_daily_tb_netcdf,write_daily_tb_netcdf
    from read_resampled_orbit import get_orbit_range,read_resampled_tbs
    from numpy_date_utils import calendar_dates_from_datetime64,convert_to_sec_in_day,convert_to_np_datetime64


    lf_1440_720 = read_land_fraction_1440_720()
    lf_1440_721 = np.zeros((721,1440))
    temp = np.zeros((719,1440))
    temp = 0.25*(lf_1440_720[0:719,:]+
                 lf_1440_720[1:720,:]+
                 np.roll(lf_1440_720[0:719,:],1,axis=1)+
                 np.roll(lf_1440_720[1:720,:],1,axis=1))
    lf_1440_721[1:720,:] = temp
    lf_1440_721[0,:] = lf_1440_721[1,:]
    lf_1440_721[720,:] = lf_1440_721[719,:]


    

    verbose = False
    satellite = 'amsr2'
    orbit = 792+14
    channel = '11H'

    num_lats = 721
    num_lons = 1440
    num_hours = 24
    num_channels = 14 #all possible AMSR2 channels 
    list_of_channels=['time','6V','6H','7V','7H','11V','11H','19V','19H',
                     '24V','24H','37V','37H','89V','89H']

    tb_array_by_hour = np.zeros((721,1440,24,num_channels),dtype=np.float32)*np.nan
    time_array_by_hour = np.zeros((721,1440,24))*np.nan

    first_orbit = True
    year = 2012
    month = 7
    day = 11
    file_list = []
    for orbit in range(792,807):

        #get the times for this orbit, and convert to time within the day...
        ob_time,filename = read_resampled_tbs(satellite=satellite,channel='time',orbit=orbit)
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

        for channel in range(5,13):
            tbs,filename = read_resampled_tbs(satellite=satellite,channel=channel,orbit=orbit)
            file_list.append(filename)
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
        
    

    print
    global_map(tb_array_by_hour[:,:,0,5],vmin=0.0,vmax=330)
    global_map(tb_array_by_hour[:,:,1,6],vmin=0.0,vmax=330)

    write_daily_tb_netcdf(year=year,month=month,day=day,satellite=satellite,
                       tb_array_by_hour=tb_array_by_hour,
                       time_array_by_hour=time_array_by_hour,
                       file_list=file_list)


    append_lf_daily_tb_netcdf(year=year,month=month,day=day,satellite=satellite,land_fraction = lf_1440_721)

    filename = get_access_output_filename(year=year,month=month,day=day,satellite=satellite,dataroot='L:/access/')
    ds = xr.open_dataset(filename)
    print

    plt.show()



    print