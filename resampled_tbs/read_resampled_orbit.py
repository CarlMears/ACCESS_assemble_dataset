
import numpy as np
import xarray as xr
import datetime
import os

def get_orbit_range(orbit:int = 10000):
    j = int((orbit-1)/5000)
    orbit_lower = 1 + j*5000
    orbit_upper = orbit_lower+4999

    return orbit_lower,orbit_upper

def read_resampled_tbs(*,satellite,channel,orbit,dataroot='L:/access/'):

    list_of_channels=['time','6V','6H','7V','7H','11V','11H','19V','19H',
                     '24V','24H','37V','37H','89V','89H']

    list_of_implemented_satellites=['amsr2']
    if satellite not in list_of_implemented_satellites:
        raise ValueError(f'Sattelite {satellite} is not implemented')
    if isinstance(channel, int):
        if ((channel < 1) or (channel > 14)):
            raise ValueError(f'channel {channel} is out of range')
        channel_str = list_of_channels[channel]
    else:
        if channel in list_of_channels:
            channel_str = channel
        else:
            raise ValueError(f'Channel {channel} not valid')
    orbit_lower,orbit_upper = get_orbit_range(orbit)
    if channel_str == 'time':
        filename = f'{dataroot}{satellite}_tb_orbits/r{orbit_lower:05d}_{orbit_upper:05d}/r{orbit:05d}.time.nc'
    else:
        filename = f'{dataroot}{satellite}_tb_orbits/r{orbit_lower:05d}_{orbit_upper:05d}/r{orbit:05d}.gridded_tbs.{channel_str}.nc'
    print(filename)
    ds = xr.open_dataset(filename)
    return ds.Data.values,filename
    print





    

