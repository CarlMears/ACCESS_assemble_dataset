import datetime
import os
from pathlib import Path
from typing import Tuple

import numpy as np
import xarray as xr
from netCDF4 import Dataset as netcdf_dataset


from era5_request.era5_requests import era5_hourly_single_level_request
from access_io.access_output import (
    get_access_output_filename,
    append_var_to_daily_tb_netcdf,
)
from util.access_interpolators import time_interpolate_synoptic_maps_ACCESS


def add_ERA5_single_level_variable_to_ACCESS_output(
    *,
    current_day: datetime.date,
    variable: Tuple[str, str],
    satellite: str,
    dataroot: Path,
    verbose: bool = False,
):
    # Get the maps of observation times from the existing output file that
    # already contains times and Tbs
    filename = get_access_output_filename(current_day, satellite, dataroot)

    with netcdf_dataset(filename, "r") as root_grp:
        try:
            times = root_grp.variables["second_since_midnight"][:, :, :].filled(
                fill_value=-999
            )
        except KeyError:
            raise ValueError(f'Error finding "second_since_midnight" in {filename}')

    # Download ERA5 data from ECMWF for all 24 hours of day, and the first hour
    # of the next day.
    next_day = current_day + datetime.timedelta(hours=24)
    try:
        file1 = era5_hourly_single_level_request(
            date=current_day,
            variable=variable[0],
            target_path=dataroot / "_temp",
            full_day=True,
        )
        file2 = era5_hourly_single_level_request(
            date=current_day,
            variable=variable[0],
            target_path=dataroot / "_temp",
            full_day=False,
        )
    except Exception:
        raise RuntimeError("Problem downloading ERA5 data using cdsapi")

    # open the files, and combine the two files into a 25-map array
    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)
    skt_first_day = ds1[variable[1]].values
    skt_next_day = ds2[variable[1]].values
    skt = np.concatenate((skt_first_day, skt_next_day), axis=0)

    # ERA-5 files are upside down relative to RSS convention.
    # TODO: I think you can just do skt = skt[:, ::-1, :] and avoid the loop
    for i in range(0, 25):
        skt[i, :, :] = np.flipud(skt[i, :, :])

    # interpolate the array of skt maps to the times in the "times" maps
    if verbose:
        print("Interpolating...")

    # list of times, each hour.
    skt_times = np.arange(0.0, 86401.0, 3600.0)

    # create output array
    skt_by_hour = np.full_like(times, np.nan)

    for hour_index in range(0, 24):
        time_map = times[:, :, hour_index]
        skt_at_time_map = time_interpolate_synoptic_maps_ACCESS(
            skt, skt_times, time_map
        )
        skt_by_hour[:, :, hour_index] = skt_at_time_map

    # write the results to the existing output file
    append_var_to_daily_tb_netcdf(
        date=current_day,
        satellite=satellite,
        var=skt_by_hour,
        var_name="skt",
        standard_name="skin_temperature",
        long_name="skin temperature interpolated from ERA5",
        valid_min=150.0,
        valid_max=400.0,
        units="Kelvin",
        v_fill=-999.0,
        dataroot=dataroot,
        overwrite=True,
    )


if __name__ == "__main__":
    date = datetime.date(2012, 7, 11)
    variable = (
        "Skin temperature",
        "skt",
    )  # need this because var name for the ERA5 request is
    # not that same as the variable name in the nc file
    # that is provided/downloaded
    satellite = "AMSR2"
    verbose = True
    if os.name == "nt":
        dataroot = Path("L:/access/")
    elif os.name == "posix":
        dataroot = Path("/mnt/ops1p-ren/l/access")

    add_ERA5_single_level_variable_to_ACCESS_output(
        current_day=date,
        variable=variable,
        satellite=satellite,
        dataroot=dataroot,
        verbose=True,
    )
