import argparse
import datetime
import os
from pathlib import Path
from typing import Tuple, Union

from netCDF4 import Dataset as netcdf_dataset
import numpy as np
from nptyping import NDArray
import xarray as xr

#import sys
#sys.path.append("/mnt/ops1p-ren/m/job_access/python/")

# from era5_request.era5_requests import era5_hourly_single_level_request
from access_io.access_output import get_access_output_filename_daily_folder
from satellite_definitions.amsr2 import SAT_NAME,REF_FREQ,REF_FREQ_mapping,REF_EIA
from rss_lock.locked_dataset import LockedDataset
# from access_io.access_output import write_daily_ancillary_var_netcdf

from access_io.access_output import set_or_create_attr
from access_io.access_attr_define import common_global_attributes_access
from access_io.access_attr_define import atm_pars_era5_attributes_access
from access_io.access_attr_define import ocean_emiss_era5_attributes_access #old version
#from access_io.access_attr_define import ocean_emiss_era5_attributes_access2 # new version, reads json file
 
from geomod10 import wind_emiss #python wrapper for geomod10b and geomod10c

# from access_io.access_attr_define import rr_imerg_attributes_access
# from access_io.access_attr_define import (
#     tcwv_era5_attributes_access,
#     tclw_era5_attributes_access,
#     skt_era5_attributes_access,
# )
# from access_io.access_attr_define import (
#     u10n_era5_attributes_access,
#     v10n_era5_attributes_access,
# )

# from util.access_interpolators import time_interpolate_synoptic_maps_ACCESS

# these are for debugging only

from rss_plotting.global_map import plot_global_map
import matplotlib.pyplot as plt

def calc_emissivity_maps(*,date,wind_source,sst_source,target_size):

    print(f"{date}")

    #Get wind info from data repository
    anc_name = f"u10n_{wind_source}"
    var_filename_final = get_access_output_filename_daily_folder(
            date, satellite.lower(), target_size, access_root, anc_name
            )
    #ds = xr.open_dataset(var_filename_final)
    with LockedDataset(var_filename_final , "r", 60) as root_grp:
        # read in existing u10n
        try:
            u10n = root_grp.variables[anc_name][:, :, :]
        except KeyError:
            raise KeyError(f"Variable: {var_name} in not dataset, can not replace")

    #Get wind info from data repository
    anc_name = f"v10n_{wind_source}"
    var_filename_final = get_access_output_filename_daily_folder(
            date, satellite.lower(), target_size, access_root, anc_name
            )
        
    with LockedDataset(var_filename_final , "r", 60) as root_grp:
        # read in existing v10n
        try:
            v10n = root_grp.variables[anc_name][:,:,:]
        except KeyError:
            raise KeyError(f"Variable: {var_name} in not dataset, can not replace")

    w10n = np.sqrt(np.square(u10n)+np.square(v10n))

    

    anc_name = f"skt_{wind_source}"
    var_filename_final = get_access_output_filename_daily_folder(
            date, satellite.lower(), target_size, access_root, anc_name
            )
    with LockedDataset(var_filename_final , "r", 60) as root_grp:
        try:
            sst = root_grp.variables[anc_name][:,:,:]
        except KeyError:
            raise KeyError(f"Variable: {var_name} in not dataset, can not replace")

    wind_ok = np.isfinite(w10n)
    w10n_ok = w10n[wind_ok]

    sst = sst - 273.15
    sst[sst < -3.0] = 3.0
    sst[sst > 34.0] = 34.0
    sst_ok = sst[wind_ok]

    phir = -999.0
    n = len(w10n_ok)
    phir = np.full_like(sst_ok,-999.0)

    num_freq = len(REF_FREQ)
    num_pol = 2
    sz = w10n.shape
    num_lat = sz[0]
    num_lon = sz[1]
    num_time = sz[2]
    shape_5d = (num_lat, num_lon, num_time, num_freq, num_pol)
    shape_3d = shape_5d[0:3]

    emiss_maps = np.full(shape_5d,np.nan,dtype=np.float32)

    for ifreq,freq in enumerate(REF_FREQ):
        print(f"Processing {freq:.3f}")
        tht=REF_EIA[ifreq]
        emiss_this_freq = wind_emiss(freq,tht,sst_ok,w10n_ok,phir)
        emiss_this_freq_v = emiss_this_freq[:,0]
        emiss_this_freq_h = emiss_this_freq[:,1]

        for ipol in [0,1]:
            maps_temp = np.full(shape_3d,np.nan,dtype=np.float32)
            maps_temp[wind_ok] = emiss_this_freq[:,ipol]
            emiss_maps[:,:,:,ifreq,ipol] = maps_temp

    return emiss_maps

def write_ocean_emiss_to_daily_ACCESS(*,
    ocean_emiss: NDArray,
    current_day: datetime.date,
    satellite: str,
    dataroot: Path,
    outputroot: Path,
    verbose: bool = False,
) -> None:
    """Add the ocean emissivity to daily ACCESS dataset.

    The base file is read in order to determine the global attributes.

    Some source of wind speed data is required.

    Currently, the only choice is era5

    If the resampled TB file doesn't exist, a `FileNotFound` error will be
    raised.

    The configured `downloader` is used to download the ERA5 data, if needed.
    """

    
    if verbose:
        print(f"Opening base file for {satellite} on {current_day} in {dataroot}")

    target_size = 0
    base_filename = get_access_output_filename_daily_folder(
                        current_day, 
                        satellite, 
                        target_size,
                        dataroot, 
                        "resamp_tbs"
                        )
    try:
    
        with netcdf_dataset(base_filename, "r") as root_grp:
            num_lats = len(root_grp["latitude"][:])
            num_lons = len(root_grp["longitude"][:])
            num_hours = len(root_grp["hours"][:])
            num_chan = len(REF_FREQ)

            glb_attrs_common = common_global_attributes_access(current_day, version="v00r00")
            glb_attrs_ocean_emiss_to_add = ocean_emiss_era5_attributes_access(satellite, version="v00r00")
            glb_attrs_ocean_emiss_to_add2 = ocean_emiss_era5_attributes_access2(satellite, version="v00r00")


            for key in glb_attrs_ocean_emiss_to_add:
                 print(f'{key}: {glb_attrs_ocean_emiss_to_add[key]}')
                
            for key in glb_attrs_ocean_emiss_to_add2:
                 print(f'{key}: {glb_attrs_ocean_emiss_to_add2[key]}')
            
            glb_attrs_ocean_emiss = glb_attrs_common | glb_attrs_ocean_emiss_to_add

            emiss_filename = get_access_output_filename_daily_folder(
                current_day, satellite, target_size,outputroot, "ocean_emiss_era5_temp"
                )
            emiss_filename_final = get_access_output_filename_daily_folder(
                current_day, satellite, target_size,outputroot, "ocean_emiss_era5"
                )
            os.makedirs(emiss_filename.parent,exist_ok=True)
            with netcdf_dataset(emiss_filename, mode="w") as trg:

                for name, dim in root_grp.dimensions.items():
                    if name in ["hours", "latitude", "longitude"]:
                        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

                # create the new dimensions needed for the emissivity data
                trg.createDimension('freq', len(REF_FREQ))
                trg.createDimension('pol', 2)

                # Copy the global attributes from the base file
                trg.setncatts({a: root_grp.getncattr(a) for a in root_grp.ncattrs()})

                for key in glb_attrs_ocean_emiss.keys():
                    value = glb_attrs_ocean_emiss[key]
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

                trg.createVariable("pol",np.int32,'pol',zlib=True)
                trg["pol"][:] = np.array([0,1],dtype=np.int32)
                trg.variables["pol"].setncatts({"pol mapping": "0 = V-pol, 1 = H-pol"})
            
                dimensions_out = ('latitude','longitude','hours','freq','pol')

                for varname, long_name, units in [
                        ("emissivity", "ocean surface emissivity", None),
                        ]:
                        print(f'starting writing {varname}')
                        least_significant_digit = 3
                        trg.createVariable(
                            varname, 
                            np.float32, 
                            dimensions_out, 
                            zlib=True,
                            least_significant_digit=least_significant_digit
                        )

                        trg[varname][:, :, :, :, :] = ocean_emiss
                        trg[varname].long_name = long_name
                        trg[varname].coordinates = "latitude longitude"
                        trg[varname].freq_names = "6, 7, 11, 19, 24, 37"
                        if units is None:
                            trg[varname].units = 'dimensionless'
                        else:
                            trg[varname].units = units
                        print(f'finished writing {varname}')
                        print()
    except FileNotFoundError:
        raise OkToSkipDay

    print()
    try:
        emiss_filename_final.unlink()
    except FileNotFoundError:
        pass

    os.rename(emiss_filename,emiss_filename_final)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=(
            "Calculate and Add Wind emissivity based on wind variables alrealy present in the ACCESS dataset."
        ),
    )
    parser.add_argument(
        "access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument(
        "access_root_output", type=Path, 
        help="directory to output data for ACCESS project.  If not provided, use access_root"
    )

    parser.add_argument("wind",choices=["era5"],help="Source of wind info")
    parser.add_argument("sst",choices=["era5"],help="Source of wind info")



    parser.add_argument(
        "start_date",
        type=datetime.date.fromisoformat,
        help="First Day to process, as YYYY-MM-DD",
    )
    parser.add_argument(
        "end_date",
        type=datetime.date.fromisoformat,
        help="Last Day to process, as YYYY-MM-DD",
    )
    parser.add_argument("sensor", choices=["amsr2"], help="Microwave sensor to use")
    parser.add_argument(
        "--verbose", help="enable more verbose screen output", action="store_true"
    )

    args = parser.parse_args()

    access_root: Path = args.access_root
    output_root: Path = args.access_root_output
    wind_source = args.wind
    sst_source = args.sst

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()

    target_size=0
    date = START_DAY
while date <= END_DAY:

    ocean_emiss = calc_emissivity_maps(date=date,
                    wind_source='era5',
                    sst_source='era5',
                    target_size=target_size)

    write_ocean_emiss_to_daily_ACCESS(
        ocean_emiss=ocean_emiss,
        current_day=date,
        satellite=satellite,
        dataroot=access_root,
        outputroot=output_root,
        verbose=True)

    date = date + datetime.timedelta(days=1.0)


    # plot_global_map(w10n[:,:,11],vmin=0.0,vmax=25.,cmap='viridis',plt_colorbar=True,title='Wind Speed (m/s)')
    # plot_global_map(sst[:,:,11],vmin=-3.0,vmax=30.,cmap='viridis',plt_colorbar=True,title='Wind Speed (m/s)')
        
    # for ifreq,freq in enumerate(REF_FREQ):
    #     plot_global_map(emiss_maps[:,:,11,ifreq,0],vmin=0.5,vmax=0.75,cmap='viridis',plt_colorbar=True,title=f'{freq:.2f}, V-Pol')
    #     plot_global_map(emiss_maps[:,:,11,ifreq,1],vmin=0.2,vmax=0.45,cmap='viridis',plt_colorbar=True,title=f'{freq:.2f}, H-Pol')
        
    # plt.show()
    print()



