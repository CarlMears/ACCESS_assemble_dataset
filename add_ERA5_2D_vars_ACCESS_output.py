import argparse
import datetime
import os
from contextlib import suppress
from pathlib import Path
from netCDF4 import Dataset as netcdf_dataset
import numpy as np

from Era5_requests.era5_requests import era5_hourly_single_level_request

# from access_io.access_output import get_access_output_filename_daily_folder
# from access_io.access_output import write_daily_ancillary_var_netcdf
from access_io.access_output_polar import write_daily_ancillary_var_netcdf_polar
from typing import Any, Tuple

import git

# from netCDF4 import Dataset as netcdf_dataset

from access_io.access_attr_define import (
    anc_var_attributes_access,
    common_global_attributes_access,
)
from access_io.access_output import (
    get_access_output_filename_daily_folder,
    write_daily_ancillary_var_netcdf,
)

from util.access_interpolators import time_interpolate_synoptic_maps_ACCESS
from util.file_times import need_to_process

from era5.resample_ERA5 import ResampleERA5


def add_ERA5_single_level_variable_to_ACCESS_output(
    *,
    current_day: datetime.date,
    variable: Tuple[str, str],
    glb_attrs: dict[str, Any],
    var_attrs: dict[str, Any],
    satellite: str,
    ksat: str,
    target_size: int,
    look: int,
    region: str,
    dataroot: Path,
    outputroot: Path,
    temproot: Path,
    verbose: bool = False,
    overwrite: bool = False,
    update: bool = False,
    script_name: str,
    commit: str,
    resampler=None,
) -> None:
    # Do some logic about the region and file names

    if region == "global":
        base_filename = get_access_output_filename_daily_folder(
            current_day, satellite, target_size, dataroot, "resamp_tbs",ksat=ksat,look=look,
        )
        anc_name = f"{variable[0]}_era5"
        var_filename_final = get_access_output_filename_daily_folder(
            current_day, satellite, target_size, dataroot, anc_name,ksat=ksat,look=look,
        )
        grid_type = "equirectangular"
        pole = "None"
    elif region in ["north", "south"]:
        pole = region
        base_filename = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            target_size,
            dataroot,
            "resamp_tbs",
            grid_type="ease2",
            pole=pole,
            ksat=ksat,
            look=look
        )
        anc_name = f"{variable[0]}_era5"
        var_filename_final = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            target_size,
            dataroot,
            anc_name,
            grid_type="ease2",
            pole=pole,
            ksat=ksat,
            look=look
        )
        grid_type = "ease2"
        pole = region
    else:
        raise ValueError(f"region {region} not valid")

    # Get the maps of observation times from the existing output file that
    # already contains times and Tbs
    # base_filename = get_access_output_filename_daily_folder(
    #     current_day, satellite.lower(), target_size, dataroot, "resamp_tbs"
    # )

    # anc_name = f"{variable[0]}_era5"
    # var_filename_final = get_access_output_filename_daily_folder(
    #     current_day, satellite.lower(), target_size, dataroot, anc_name
    # )

    if need_to_process(
        date=current_day,
        satellite=satellite,
        ksat=ksat,
        target_size=target_size,
        dataroot=dataroot,
        outputroot=outputroot,
        var=anc_name,
        overwrite=overwrite,
        update=update,
        grid_type=grid_type,
        pole=pole,
        look=look,
        verbose=True,
    ):
        if not base_filename.is_file():
            print(f"base file for {current_day} does not exist, skipping")
            return

        if not (overwrite or update):
            if var_filename_final.is_file():
                print(
                    f"{variable[0]} file for {current_day} exists, skipping to next day"
                )
                return
        else:
            with suppress(FileNotFoundError):
                var_filename_final.unlink()

        print(f"Processing: {anc_name} for {current_day}")
        try:
            # read in base file and extract dimensions and make sure time is available
            with netcdf_dataset(base_filename, "r") as root_grp:
                try:
                    times = root_grp.variables["time"][:, :, :]
                    times = (
                        times
                        - (current_day - datetime.date(1900, 1, 1)).total_seconds()
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
            variable=variable[1],
            target_path=temproot,
            full_day=True,
            full_month=True,
        )

        # if next day is in the same month, this second request
        # refers to the same file, so no second download will be done
        file2 = era5_hourly_single_level_request(
            date=next_day,
            variable=variable[1],
            target_path=temproot,
            full_day=True,
            full_month=True,
        )

        # except Exception:
        #    raise RuntimeError("Problem downloading ERA5 data using cdsapi")

        # open the file(s), and combine the two files into a
        # 25-map array for the day being processed

        if current_day.month == next_day.month:
            hour_index1 = 24 * (current_day.day - 1)
            hour_index2 = hour_index1 + 25
            ds1 = netcdf_dataset(file1)
            var = ds1[variable[0]][hour_index1:hour_index2, :, :]

        else:
            # This is the case when the 25th hour is in the next month
            hour_index1 = 24 * (current_day.day - 1)
            hour_index2 = hour_index1 + 24
            ds1 = netcdf_dataset(file1)
            var_first_day = ds1[variable[0]][hour_index1:hour_index2, :, :]

            ds2 = netcdf_dataset(file2)
            var_next_day = ds2[variable[0]][0, :, :]
            var = np.concatenate(
                (var_first_day, var_next_day[np.newaxis, :, :]), axis=0
            )

        var = np.flip(var, 1)
        # file1 modification time as a datetime.datetime object
        mod_time = datetime.datetime.utcfromtimestamp(file1.stat().st_mtime)

        resample_required = True
        if (target_size == 30) and (grid_type == "equirectangular"):
            resample_required = False

        if resample_required:
            if resampler is None:
                # should be done at a higher level, but just in case....
                resampler = ResampleERA5(target_size=target_size, region=region)

            # time_begin = datetime.datetime.now()
            var = resampler.resample_fortran(var)
            # time = datetime.datetime.now() - time_begin
            # print(time)

        # list of times, each hour.
        var_times = np.arange(0.0, 86401.0, 3600.0)

        # create output array

        var_by_hour = np.full_like(times, np.nan).filled()

        # time_begin = datetime.datetime.now()
        for hour_index in range(0, 24):
            time_map = times[:, :, hour_index]
            var_at_time_map = time_interpolate_synoptic_maps_ACCESS(
                var, var_times, time_map
            )
            var_by_hour[:, :, hour_index] = var_at_time_map
        # time = datetime.datetime.now() - time_begin

        if "_FillValue" in var_attrs.keys():
            var_by_hour[~np.isfinite(var_by_hour)] = var_attrs["_FillValue"]

        glb_attrs["date_accessed"] = f"{mod_time}"
        glb_attrs["id"] = var_filename_final.name
        glb_attrs["corresponding_resampled_Tb_file"] = base_filename.name
        glb_attrs["script_name"] = script_name
        glb_attrs["commit"] = commit

        # write the results to a separate output file
        if grid_type == "equirectangular":
            write_daily_ancillary_var_netcdf(
                date=current_day,
                satellite=satellite,
                target_size=target_size,
                look=look,
                ksat=ksat,
                anc_data=var_by_hour,
                anc_name=anc_name,
                anc_attrs=var_attrs,
                global_attrs=glb_attrs,
                dataroot=outputroot,
            )
        elif grid_type == "ease2":
            write_daily_ancillary_var_netcdf_polar(
                date=current_day,
                satellite=satellite,
                ksat=ksat,
                target_size=target_size,
                look=look,
                grid_type=grid_type,
                pole=pole,
                anc_data=var_by_hour,
                anc_name=anc_name,
                anc_attrs=var_attrs,
                global_attrs=glb_attrs,
                dataroot=outputroot,
            )
        else:
            raise ValueError(f"Grid Type {grid_type} is not valid")
    else:
        print(f"No processing needed for {anc_name} on {current_day}")


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
        "--access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument(
        "--output_root", type=Path, help="Root directory to write output"
    )
    parser.add_argument(
        "--temp_root", type=Path, help="Root directory store temporary files"
    )
    parser.add_argument(
        "--start_date",
        type=datetime.date.fromisoformat,
        help="First Day to process, as YYYY-MM-DD",
    )
    parser.add_argument(
        "--end_date",
        type=datetime.date.fromisoformat,
        help="Last Day to process, as YYYY-MM-DD",
    )
    parser.add_argument("--sensor", choices=["amsr2","ssmi", "smap"], help="Microwave sensor to use")
    parser.add_argument("--ksat", choices=["13", "15"],help="SSMI Satellite Number")
    parser.add_argument(
        "--target_size", choices=["30", "70"], help="Size of target footprint in km"
    )
    parser.add_argument("--look", choices=["0", "1"], default="0", help="Look direction")
    parser.add_argument(
        "--region", help="region to process", choices=["global", "north", "south"]
    )
    parser.add_argument("--version", help="version sting - e.g. v01r00")

    parser.add_argument("-v", "--variables", nargs="*", default=[])

    parser.add_argument(
        "--overwrite", help="overwrite exisitng files", action="store_true"
    )
    parser.add_argument(
        "--update",
        help="overwrite existing files is older than base file",
        action="store_true",
    )
    parser.add_argument(
        "--verbose", help="enable more verbose screen output", action="store_true"
    )

    args = parser.parse_args()

    access_root: Path = args.access_root
    output_root: Path = args.output_root
    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()
    ksat = args.ksat
    target_size = int(args.target_size)
    look = int(args.look)
    region = args.region
    version = args.version
    var_list = args.variables
    overwrite = args.overwrite
    update = args.update
    print(var_list)
    script_name = parser.prog
    repo = git.Repo(search_parent_directories=True)
    commit = repo.head.object.hexsha

    for var in var_list:
        print("Adding ERA5 2D variables")
        print(f"ACCESS Root: {access_root}")
        print(f"Date Range: {START_DAY} - {END_DAY}")
        print(f"Satellite: {satellite}")
        print(f"Target Size: {target_size}")
        print(f"Look: {look}")
        print(f"Version: {version}")
        print(f"Variable: {var}")
        print()

        var_dict = {
            "skt": "Skin temperature",
            "tcwv": "Total column water vapour",
            "tclw": "total_column_cloud_liquid_water",
            "u10n": "10m_u_component_of_neutral_wind",
            "v10n": "10m_v_component_of_neutral_wind",
        }
        try:
            variable = (var, var_dict[var])
        except KeyError:
            print(f"Variable {var} not defined - skipping")
            continue

        resample_required = True
        if (target_size == 30) and (region == "global"):
            resample_required = False

        resampler = None
        if resample_required:
            if resampler is None:
                resampler = ResampleERA5(target_size=target_size, region=region)

        date = START_DAY
        while date <= END_DAY:
            # common global_attributes for the project
            glb_attrs = common_global_attributes_access(
                date, satellite, target_size, version
            )

            # variable-specific attributes
            var_attrs = anc_var_attributes_access(satellite, var + "_era5", version)

            # add the global part of these to the global_attrs
            glb_attrs.update(var_attrs["global"])

            # keep the variable decription parts in var_attrs

            var_attrs = var_attrs["var"]

            add_ERA5_single_level_variable_to_ACCESS_output(
                current_day=date,
                variable=variable,
                var_attrs=var_attrs,
                glb_attrs=glb_attrs,
                satellite=satellite,
                ksat=ksat,
                target_size=target_size,
                look=look,
                region=region,
                dataroot=access_root,
                outputroot=output_root,
                temproot=temp_root,
                verbose=args.verbose,
                overwrite=overwrite,
                update=update,
                script_name=script_name,
                commit=commit,
                resampler=resampler,
            )

            date += datetime.timedelta(days=1)
