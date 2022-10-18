"""Run the RTM and add the atmospheric terms to the resampled TB file.

ERA5 data is downloaded if missing.
"""

import argparse
from http.client import OK
import os
from datetime import date
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from netCDF4 import Dataset
import datetime

from access_io.access_output import get_access_output_filename_daily_folder
from access_io.access_output import set_or_create_attr
from access_io.access_attr_define import common_global_attributes_access
from access_io.access_attr_define import atm_pars_era5_attributes_access

from util.access_interpolators import time_interpolate_synoptic_maps_ACCESS
# Reference frequencies (in GHz) to use
REF_FREQ = np.array([6.9, 7.3, 10.7, 18.7, 23.8, 37.0], np.float32)

REF_FREQ_mapping = np.array([1,1,2,3,4,5],np.int32)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
REF_EIA = np.array([53.0, 53.0, 53.0, 53.0, 53.0, 53.0], np.float32)

class OkToSkipDay(Exception):
    pass


class DailyRtm:
    """RTM results for the entire day."""

    def __init__(
                self,
                date_to_load: datetime.datetime,data_root: Path):

        path_to_data = data_root / f'era5_tbs_{date_to_load.year}-{date_to_load.month:02d}-{date_to_load.day:02d}.nc'
        with Dataset(path_to_data, "r") as f:
            tb_down = f['tb_down'][:,:,:,:]
            tb_up = f['tb_up'][:,:,:,:]
            trans = f['tran'][:,:,:,:]

        dims = tb_down.shape
        num_lats = dims[1]
        num_lons = dims[2]
        num_hours = dims[0]
        num_freqs = dims[3]

        shape_rtm_4d = (num_lats,num_lons,num_hours+1,num_freqs)

        self.downwelling_tb = np.full(shape_rtm_4d,np.nan,dtype=np.float32)
        self.upwelling_tb = np.full(shape_rtm_4d,np.nan,dtype=np.float32)
        self.transmissivity = np.full(shape_rtm_4d,np.nan,dtype=np.float32)

        for hour in range(24):
            for freq in range(num_freqs-1):
                self.downwelling_tb[:,:,hour,freq] = tb_down[hour,:,:,freq]
                self.upwelling_tb[:,:,hour,freq] = tb_up[hour,:,:,freq]
                self.transmissivity[:,:,hour,freq] = trans[hour,:,:,freq]

        date_to_load_plus_one = date_to_load + datetime.timedelta(days=1)
        path_to_data_plus_one = data_root / f'era5_tbs_{date_to_load_plus_one.year}-{date_to_load_plus_one.month:02d}-{date_to_load_plus_one.day:02d}.nc'
        
        with Dataset(path_to_data_plus_one, "r") as f:
            tb_down = f['tb_down'][:,:,:,:]
            tb_up = f['tb_up'][:,:,:,:]
            trans = f['tran'][:,:,:,:]

        for hour in [0]:
            for freq in range(num_freqs-1):
                self.downwelling_tb[:,:,hour+24,freq] = tb_down[hour,:,:,freq]
                self.upwelling_tb[:,:,hour+24,freq] = tb_up[hour,:,:,freq]
                self.transmissivity[:,:,hour+24,freq] = trans[hour,:,:,freq]

        self.time_in_day = np.arange(0,25)*3600.0

def write_atmosphere_to_daily_ACCESS(
    current_day: date,
    satellite: str,
    dataroot: Path,
    temproot: Path,
    verbose: bool = False,
) -> None:
    """Append the atmospheric terms to daily ACCESS resampled-TB dataset.

    The file is read in order to determine the valid points where the RTM should
    be called. The resulting values are written by appending the new variables
    to the file.

    ERA5 data, both on the surface and as profiles, is required.

    If the resampled TB file doesn't exist, a `FileNotFound` error will be
    raised.

    The configured `downloader` is used to download the ERA5 data, if needed.
    """
    if verbose:
        print(f"Opening data for {satellite} on {current_day} in {dataroot}")

    base_filename = get_access_output_filename_daily_folder(
                        current_day, 
                        satellite, 
                        dataroot, 
                        "resamp_tbs"
                        )
    try:
    
        with Dataset(base_filename, "r") as root_grp:

            if verbose:
                print(f"Reading ERA5 computed RTM data {satellite} on {current_day} in {temproot}") 

            num_lats = len(root_grp["latitude"][:])
            num_lons = len(root_grp["longitude"][:])
            num_hours = len(root_grp["hours"][:])
            num_chan = len(REF_FREQ)
                
            rtm_data = DailyRtm(current_day,temproot)

            glb_attrs_common = common_global_attributes_access(current_day, version="v00r00")
            glb_attrs_atm_to_add = atm_pars_era5_attributes_access(satellite, version="v00r00")
            glb_attrs_atm = glb_attrs_common | glb_attrs_atm_to_add
            atm_filename = get_access_output_filename_daily_folder(
                current_day, satellite, dataroot, "atm_par_era5_temp"
                )
            atm_filename_final = get_access_output_filename_daily_folder(
                current_day, satellite, dataroot, "atm_par_era5"
                )


            with Dataset(atm_filename, mode="w") as trg:

                for name, dim in root_grp.dimensions.items():
                    if name in ["hours", "latitude", "longitude"]:
                        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

                trg.createDimension('freq', len(REF_FREQ))

                # Copy the global attributes from the base file
                trg.setncatts({a: root_grp.getncattr(a) for a in root_grp.ncattrs()})

                for key in glb_attrs_atm.keys():
                    value = glb_attrs_atm[key]
                    set_or_create_attr(trg, key, value)

                for var_name in ["time", "hours", "longitude", "latitude"]:
                    # Create the time and dimension variables in the output file
                    var_in = root_grp[var_name]
                    trg.createVariable(var_name, var_in.dtype, var_in.dimensions, zlib=True)
                    
                    # Copy the attributes
                    trg.variables[var_name].setncatts(
                        {a: var_in.getncattr(a) for a in var_in.ncattrs()}
                    )
                    trg[var_name][:] = var_in[:]
            
                dimensions_out = ('latitude','longitude','hours','freq')

                for varname, long_name, units in [
                        ("transmissivity", "atmospheric transmissivity", None),
                        ("upwelling_tb", "upwelling brightness temperature", "kelvin"),
                        ("downwelling_tb", "downwelling brightness temperature", "kelvin"),
                        ]:
                        print(f'starting writing {varname}')
                        if varname == "transmissivity":
                            least_significant_digit = 3
                        else:
                            least_significant_digit = 2
                        trg.createVariable(
                            varname, 
                            np.float32, 
                            dimensions_out, 
                            zlib=True,
                            least_significant_digit=least_significant_digit
                        )
                        for freq_index,freq in enumerate(REF_FREQ):
                            var = getattr(rtm_data,varname)[:,:,:,REF_FREQ_mapping[freq_index]]
                            var = np.moveaxis(var, -1, 0)
                            var_times = rtm_data.time_in_day
                            for hour_index in range(len(root_grp["hours"][:])):
                                time_map = root_grp['time'][:, :, hour_index]
                                time_map = (time_map - (current_day - datetime.date(1900, 1, 1)).total_seconds())
                                var_at_time_map = time_interpolate_synoptic_maps_ACCESS(
                                    var, var_times, time_map
                                    )
                                trg[varname][:, :, hour_index,freq_index] = var_at_time_map
                        trg[varname].long_name = long_name
                        trg[varname].coordinates = "latitude longitude"
                        trg[varname].freq_names = "6, 7, 11, 19, 24, 37"
                        if units is None:
                            trg[varname].units = 'dimensionless'
                        else:
                            trg[varname].units = units
                        print(f'finished writing {varname}')
    except FileNotFoundError:
        raise OkToSkipDay

    print()
    try:
        atm_filename_final.unlink()
    except FileNotFoundError:
        pass

    os.rename(atm_filename,atm_filename_final)
  


    # trg = Dataset(atm_filename, mode="w")

    # glb_attrs_common = common_global_attributes_access(date, version="v00r00")

    # glb_attrs_tbup_to_add = tbup_era5_attributes_access(satellite,version="v00r00")
    # glb_attrs_tbup = glb_attrs_common | glb_attrs_tbup_to_add

    # glb_attrs_tbdown_to_add = tbdown_era5_attributes_access(satellite,version="v00r00")
    # glb_attrs_tbdown = glb_attrs_common | glb_attrs_tbdown_to_add

    # glb_attrs_trans_to_add = trans_era5_attributes_access(satellite,version="v00r00")
    # glb_attrs_trans = glb_attrs_common | glb_attrs_trans_to_add


    
    #         # Copy the values
    #         trg.variables[var_name][:] = root_grp.variables[var_name][:]


    #     for hour in range(24):
    #     if verbose:
    #         print(f"Processing hour {hour+1}/24")



if __name__ == "__main__":
    cds_help = (
        "For downloading ERA5 data from CDS, the UID and API key "
        "must be set as arguments or in the 'CDS_UID' and 'CDS_API_KEY` "
        "environment variables"
    )
    parser = argparse.ArgumentParser(
        description=(
            "Compute and append atmospheric RTM terms to ACCESS output file. "
            "ERA5 data is downloaded if required."
        ),
        epilog=cds_help,
    )
    parser.add_argument(
        "access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument(
        "start_date", type=date.fromisoformat, help="Day to process, as YYYY-MM-DD"
    )
    parser.add_argument(
        "end_date", type=date.fromisoformat, help="Day to process, as YYYY-MM-DD"
    )
    parser.add_argument("sensor", choices=["amsr2"], help="Microwave sensor to use")
    
    args = parser.parse_args()

    access_root: Path = args.access_root
    rtm_dir = access_root.parent / "_temp" / "rtm"
    date_to_do = args.start_date   
    while date_to_do <= args.end_date:
        try: 
            write_atmosphere_to_daily_ACCESS(
                date_to_do, args.sensor, access_root, rtm_dir, verbose=True
                )
        except OkToSkipDay:
            print(f'Problem finding file...skipping {date_to_do}')
        date_to_do = date_to_do + datetime.timedelta(days=1)
