import argparse
import datetime
import os
from pathlib import Path
from typing import Tuple

from netCDF4 import Dataset as netcdf_dataset
import numpy as np

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
    temproot: Path,
    verbose: bool = False,
    force_overwrite: bool = False,
) -> None:
    # Get the maps of observation times from the existing output file that
    # already contains times and Tbs
    filename = get_access_output_filename(current_day, satellite, dataroot)

    try:
        with netcdf_dataset(filename, "r") as root_grp:
            # check to see if variable already exists
            if not force_overwrite:
                try:
                    root_grp.variables[variable[1]][:, :, :].filled(fill_value=-999)
                    print(
                        f"var {variable[0]} ({variable[1]}) already exists.  skipping.."
                    )
                    return
                except KeyError:
                    # we expect a key error if variable is needed
                    pass
            try:
                times = root_grp.variables["time"][:, :, :]
                times = (
                    times - (current_day - datetime.date(1900, 1, 1)).total_seconds()
                )
            except KeyError:
                raise ValueError(f'Error finding "second_since_midnight" in {filename}')
    except FileNotFoundError:
        print(f"File: {filename} not found, skipping")
        return
    # Download ERA5 data from ECMWF for all 24 hours of day, and the first hour
    # of the next day.
    next_day = current_day + datetime.timedelta(hours=24)
    # try:
    os.makedirs(temproot, exist_ok=True)
    file1 = era5_hourly_single_level_request(
        date=current_day,
        variable=variable[0],
        target_path=temproot,
        full_day=True,
        full_month=True,
    )

    # if next day is in the same month, this second request
    # refers to the same file, so no second download will be done
    file2 = era5_hourly_single_level_request(
        date=next_day,
        variable=variable[0],
        target_path=temproot,
        full_day=True,
        full_month=True,
    )

    # except Exception:
    #    raise RuntimeError("Problem downloading ERA5 data using cdsapi")

    # open the file(s), and combine the two files into a 25-map array
    if current_day.month == next_day.month:
        hour_index1 = 24 * (current_day.day - 1)
        hour_index2 = hour_index1 + 25
        ds1 = netcdf_dataset(file1)
        skt = ds1[variable[1]][hour_index1:hour_index2, :, :]
    else:
        hour_index1 = 24 * (current_day.day - 1)
        hour_index2 = hour_index1 + 24
        ds1 = netcdf_dataset(file1)
        skt_first_day = ds1[variable[1]][hour_index1:hour_index2, :, :]

        ds2 = netcdf_dataset(file2)
        skt_next_day = ds2[variable[1]][0, :, :]
        skt = np.concatenate((skt_first_day, skt_next_day[np.newaxis, :, :]), axis=0)

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
        units="kelvin",
        v_fill=-999.0,
        dataroot=dataroot,
        overwrite=True,
    )


if __name__ == "__main__":

    cds_help = (
        "For downloading ERA5 data from CDS, the UID and API key "
        "must be set as arguments or in the 'CDS_UID' and 'CDS_API_KEY` "
        "environment variables"
    )
    parser = argparse.ArgumentParser(
        description=(
            "Interpolate and append surface temperature to ACCESS output file. "
            "ERA5 data is downloaded if required."
        ),
        epilog=cds_help,
    )
    parser.add_argument(
        "access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument(
        "temp_root", type=Path, help="Root directory store temporary files"
    )
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
    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()

    date = START_DAY
    while date <= END_DAY:
        print(f"{date}")

        # need this because var name for the ERA5 request is not that same as
        # the variable name in the nc file that is provided/downloaded
        variable = ("Skin temperature", "skt")
        add_ERA5_single_level_variable_to_ACCESS_output(
            current_day=date,
            variable=variable,
            satellite=satellite,
            dataroot=access_root,
            temproot=temp_root,
            verbose=args.verbose,
        )

        date += datetime.timedelta(days=1)
