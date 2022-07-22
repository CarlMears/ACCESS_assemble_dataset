import argparse
import datetime
import os
from pathlib import Path
from typing import Tuple

from netCDF4 import Dataset as netcdf_dataset
import numpy as np

from era5_request.era5_requests import era5_hourly_single_level_request
from access_io.access_output import get_access_output_filename_daily_folder
from access_io.access_output import write_daily_ancillary_var_netcdf
from util.access_interpolators import time_interpolate_synoptic_maps_ACCESS


def add_ERA5_single_level_variable_to_ACCESS_output(
    *,
    current_day: datetime.date,
    variable: Tuple[str, str],
    var_attrs: dict,
    satellite: str,
    dataroot: Path,
    temproot: Path,
    verbose: bool = False,
    force_overwrite: bool = False,
) -> None:
    # Get the maps of observation times from the existing output file that
    # already contains times and Tbs
    base_filename = get_access_output_filename_daily_folder(
        current_day, satellite.lower(), dataroot, "resamp_tbs"
    )

    anc_name = f"{variable[1]}_era5"
    var_filename_final = get_access_output_filename_daily_folder(
        current_day, satellite.lower(), dataroot, anc_name
    )

    if not base_filename.is_file():
        print(f"base file for {current_day} does not exist, skipping")
        return

    if var_filename_final.is_file():
        if not force_overwrite:
            print(f"{variable[0]} file for {current_day} exists, skipping to next day")
            return

    try:
        # read in base file and extract dimensions and make sure time is available
        with netcdf_dataset(base_filename, "r") as root_grp:
            try:
                times = root_grp.variables["time"][:, :, :]
                times = (
                    times - (current_day - datetime.date(1900, 1, 1)).total_seconds()
                )
            except KeyError:
                raise ValueError(f'Error finding "time" in {base_filename}')
    except FileNotFoundError:
        print(f"File: {base_filename} not found, skipping")
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

    # open the file(s), and combine the two files into a 25-map array for the day being processed

    if current_day.month == next_day.month:
        hour_index1 = 24 * (current_day.day - 1)
        hour_index2 = hour_index1 + 25
        ds1 = netcdf_dataset(file1)
        var = ds1[variable[1]][hour_index1:hour_index2, :, :]
        
    else:
        #This is the case when the 25th hour is in the next month
        hour_index1 = 24 * (current_day.day - 1)
        hour_index2 = hour_index1 + 24
        ds1 = netcdf_dataset(file1)
        var_first_day = ds1[variable[1]][hour_index1:hour_index2, :, :]

        ds2 = netcdf_dataset(file2)
        var_next_day = ds2[variable[1]][0, :, :]
        var = np.concatenate((var_first_day, var_next_day[np.newaxis, :, :]), axis=0)

    #file1 modification time as a datetime.datetime object
    mod_time = datetime.datetime.utcfromtimestamp(file1.stat().st_mtime)

    # ERA-5 files are upside down relative to RSS convention.
    # TODO: I think you can just do var = var[:, ::-1, :] and avoid the loop
    for i in range(0, 25):
        var[i, :, :] = np.flipud(var[i, :, :])

    # interpolate the array of var maps to the times in the "times" maps

    print(f"Interpolating {variable[0]}")

    # list of times, each hour.
    var_times = np.arange(0.0, 86401.0, 3600.0)

    # create output array
    var_by_hour = np.full_like(times, np.nan)

    for hour_index in range(0, 24):
        time_map = times[:, :, hour_index]
        var_at_time_map = time_interpolate_synoptic_maps_ACCESS(
            var, var_times, time_map
        )
        var_by_hour[:, :, hour_index] = var_at_time_map

    var_attrs['Date accessed'] = f'{mod_time}'

    # write the results to the existing output file
    write_daily_ancillary_var_netcdf(
        date=current_day,
        satellite=satellite,
        anc_data=var_by_hour,
        anc_name=anc_name,
        anc_attrs=var_attrs,
        dataroot=dataroot
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
        source = (
            "Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz-Sabater, J., et al. "
               "(2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), "
               "1999–2049. https://doi.org/10.1002/qj.3803. "
               "ERA5 hourly data on single levels from 1959 to present. Skin Temperature"
               "0.25 degree x 0.25 degree gridded data downloaded from the Copernicus Climate Data Store. "
               "Dataset DOI: 10.24381/cds.adbb2d47 ")

        var_attrs = dict(standard_name="surface_temperature",
                     long_name="skin temperature interpolated from hourly ERA5 output",
                     valid_min=150.0,
                     valid_max=400.0,
                     units="kelvin",
                     v_fill=-999.0,
                     source=source
                        )
        add_ERA5_single_level_variable_to_ACCESS_output(
            current_day=date,
            variable=variable,
            var_attrs=var_attrs,
            satellite=satellite,
            dataroot=access_root,
            temproot=temp_root,
            verbose=args.verbose,
            force_overwrite=True
        )

        variable = ("Total column water vapour", "tcwv")
        
        source = (
               "Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz-Sabater, J., et al. "
               "(2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), "
               "1999–2049. https://doi.org/10.1002/qj.3803. "
               "ERA5 hourly data on single levels from 1959 to present. Total column water vapour"
               "0.25 degree x 0.25 degree gridded data downloaded from the Copernicus Climate Data Store. "
               "Dataset DOI: 10.24381/cds.adbb2d47 ")

        var_attrs = dict(standard_name="atmosphere_mass_content_of_water_vapor",
                     long_name="Total column water vapour interpolated from hourly ERA5 output",
                     valid_min=0.0,
                     valid_max=120.0,
                     units="kg/m^2",
                     v_fill=-999999.0,
                     source=source
                        )
        add_ERA5_single_level_variable_to_ACCESS_output(
            current_day=date,
            variable=variable,
            var_attrs=var_attrs,
            satellite=satellite,
            dataroot=access_root,
            temproot=temp_root,
            verbose=args.verbose,
            force_overwrite=True
        )

        variable = ("total_column_cloud_liquid_water", "tclw")
        
        source = (
               "Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz-Sabater, J., et al. "
               "(2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), "
               "1999–2049. https://doi.org/10.1002/qj.3803. "
               "ERA5 hourly data on single levels from 1959 to present. Total column cloud liquid water"
               "0.25 degree x 0.25 degree gridded data downloaded from the Copernicus Climate Data Store. "
               "Dataset DOI: 10.24381/cds.adbb2d47 ")

        var_attrs = dict(standard_name="atmosphere_mass_content_of_cloud_liquid_water",
                     long_name="Total column cloud liquid water interpolated from hourly ERA5 output",
                     valid_min=0.0,
                     valid_max=20.0,
                     units="kg/m^2",
                     v_fill=-999999.0,
                     source=source
                        )
        add_ERA5_single_level_variable_to_ACCESS_output(
            current_day=date,
            variable=variable,
            var_attrs=var_attrs,
            satellite=satellite,
            dataroot=access_root,
            temproot=temp_root,
            verbose=args.verbose,
            force_overwrite=True
        )

        variable = ("10m_u_component_of_neutral_wind", "u10n")
        
        source = (
               "Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz-Sabater, J., et al. "
               "(2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), "
               "1999–2049. https://doi.org/10.1002/qj.3803. "
               "ERA5 hourly data on single levels from 1959 to present. 10m_u_component_of_neutral_wind"
               "0.25 degree x 0.25 degree gridded data downloaded from the Copernicus Climate Data Store. "
               "Dataset DOI: 10.24381/cds.adbb2d47 ")

        var_attrs = dict(standard_name="eastward_wind",
                     long_name="10m u component of neutral wind interpolated from hourly ERA5 output",
                     valid_min=-100.0,
                     valid_max=100.0,
                     units="m/s",
                     v_fill=-999999.0,
                     source=source
                        )
        add_ERA5_single_level_variable_to_ACCESS_output(
            current_day=date,
            variable=variable,
            var_attrs=var_attrs,
            satellite=satellite,
            dataroot=access_root,
            temproot=temp_root,
            verbose=args.verbose,
            force_overwrite=True
        )

        variable = ("10m_v_component_of_neutral_wind", "v10n")
        
        source = (
               "Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz-Sabater, J., et al. "
               "(2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), "
               "1999–2049. https://doi.org/10.1002/qj.3803. "
               "ERA5 hourly data on single levels from 1959 to present. 10m_v_component_of_neutral_wind"
               "0.25 degree x 0.25 degree gridded data downloaded from the Copernicus Climate Data Store. "
               "Dataset DOI: 10.24381/cds.adbb2d47 ")

        var_attrs = dict(standard_name="northward_wind",
                     long_name="10m v component of neutral wind interpolated from hourly ERA5 output",
                     valid_min=-100.0,
                     valid_max=100.0,
                     units="m/s",
                     v_fill=-999999.0,
                     source=source
                        )
        add_ERA5_single_level_variable_to_ACCESS_output(
            current_day=date,
            variable=variable,
            var_attrs=var_attrs,
            satellite=satellite,
            dataroot=access_root,
            temproot=temp_root,
            verbose=args.verbose,
            force_overwrite=True
        )


        date += datetime.timedelta(days=1)
