"""Add IMERG rain rates to an existing daily ACCESS data file."""

import argparse
from contextlib import suppress
import datetime
from pathlib import Path
from typing import Optional
import git
import numpy as np
from netCDF4 import Dataset as netcdf_dataset

from access_io.access_attr_define import (
    anc_var_attributes_access,
    common_global_attributes_access,
)
from access_io.access_output import (
    write_daily_ancillary_var_netcdf,
    edit_attrs_daily_ancillary_var_netcdf,
    get_access_output_filename_daily_folder,
    write_daily_ancillary_var_netcdf,
)
from access_io.access_output_polar import write_daily_ancillary_var_netcdf_polar
from imerg_request.imerg_requests import imerg_half_hourly_request
from resampling_utils.resample_imerg_polar import ResampleIMERG
from resampling_utils.imerg_resampling_routines import resample_imerg_day
from util.file_times import need_to_process
import git


def redo_imerg_rain_rate_attrs_ACCESS(
    *,
    current_day: datetime.date,
    satellite: str,
    dataroot: Path,
    outputroot: Path,
    temproot: Path,
    footprint_diameter_km: int,
    overwrite: bool = False,
    script_name: str,
    commit: str,
) -> None:
    imerge_filename_final = get_access_output_filename_daily_folder(
        current_day, satellite, footprint_diameter_km, outputroot, "rainfall_rate"
    )

    if imerge_filename_final.is_file():
        version = "v01r00"
        rr_attrs = anc_var_attributes_access(
            satellite, "rain_rate_imerg", version=version
        )
        global_attrs = common_global_attributes_access(
            current_day, satellite, footprint_diameter_km, version=version
        )

        # add in the global part of the var-specific attributes
        global_attrs.update(rr_attrs["global"])

        # get the variable description parts of the var-specific attributes
        var_attrs = rr_attrs["var"]

        edit_attrs_daily_ancillary_var_netcdf(
            date=current_day,
            satellite=satellite,
            target_size=footprint_diameter_km,
            anc_data=None,
            anc_name="rainfall_rate",
            anc_attrs=var_attrs,
            global_attrs=global_attrs,
            dataroot=dataroot,
        )


def write_imerg_rain_rate_for_ACCESS(
    *,
    current_day: datetime.date,
    satellite: str,
    dataroot: Path,
    outputroot: Path,
    temproot: Path,
    footprint_diameter_km: int,
    region: str,
    overwrite: Optional[bool] = False,
    update: Optional[bool] = False,
    script_name: str,
    commit: str,
    resampler: Optional[ResampleIMERG],
) -> None:
    
    if region == "global":
        base_filename = get_access_output_filename_daily_folder(
            current_day, satellite, footprint_diameter_km, dataroot, "resamp_tbs")
        var = "rainfall_rate"
        imerge_filename_final = get_access_output_filename_daily_folder(
            current_day, satellite, footprint_diameter_km, outputroot, var
        )
        grid_type = "equirectangular"
        pole = None
    elif region in ["north", "south"]:
        pole = region
        base_filename = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            footprint_diameter_km,
            dataroot,
            "resamp_tbs",
            grid_type="ease2",
            pole=pole,
        )
        var = "rainfall_rate"
        imerge_filename_final = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            footprint_diameter_km,
            dataroot,
            var,
            grid_type="ease2",
            pole=pole,
        )
        grid_type = "ease2"
        pole = region
    else:
        raise ValueError(f"region {region} not valid")

    if need_to_process(
        date=current_day,
        satellite=satellite,
        target_size=footprint_diameter_km,
        dataroot=dataroot,
        outputroot=outputroot,
        var=var,
        overwrite=overwrite,
        update=update,
        grid_type=grid_type,
        pole=pole,
    ):
        with suppress(FileNotFoundError):
            imerge_filename_final.unlink()

        # try:
        #     # read in base file and extract dimensions and metadata
        with netcdf_dataset(base_filename, "r") as root_grp:
            try:
                times = root_grp.variables["time"][:, :, :]
                times = (
                    times - (current_day - datetime.date(1900, 1, 1)).total_seconds()
                )
            except KeyError:
                raise ValueError(f'Error finding "time" in {base_filename}')
        # except FileNotFoundError:
        #     print(f"File: {base_filename} not found, skipping")
        #     return

        # Downloding all IMERG files for the day
        try:
            imerg_half_hourly_request(
                date=current_day,
                target_path=temproot / "imerg",
            )

        except Exception as e:
            raise RuntimeError("Problem downloading IMERG data") from e

        # An array of hour times in seconds
        hourly_intervals = np.arange(0, 86401, 3600)

        # Return resampled rain rate maps for each hour of the day
        rr_for_access, mod_time = resample_imerg_day(
            np.roll(times, 720, axis=1),
            hourly_intervals,
            date,
            footprint_diameter_km,
            region,
            target_path=temproot / "imerg",
            resampler=resampler,
        )

        rr_for_access = np.roll(rr_for_access, 720, axis=1)

        version = "v01r00"
        rr_attrs = anc_var_attributes_access(
            satellite, "rain_rate_imerg", version=version
        )
        global_attrs = common_global_attributes_access(
            date, satellite, footprint_diameter_km, version=version
        )

        # add in the global part of the var-specific attributes
        global_attrs.update(rr_attrs["global"])

        # get the variable description parts of the var-specific attributes
        var_attrs = rr_attrs["var"]

        global_attrs["date_accessed"] = f"{datetime.datetime.today()}"
        global_attrs["id"] = imerge_filename_final.name
        global_attrs["corresponding_resampled_Tb_file"] = base_filename.name
        global_attrs["cell_method"] = (
            "time: closest 30-min IMERG file; "
            f"area: weighted average over {footprint_diameter_km}km Gaussian footprint"
        )
        global_attrs["corresponding_resampled_tb_file"] = base_filename.name
        global_attrs["commit"] = commit
        global_attrs["script_name"] = script_name

        # write the results to a separate output file
        if grid_type == "equirectangular":
            write_daily_ancillary_var_netcdf(
                date=date,
                satellite=satellite,
                target_size=footprint_diameter_km,
                anc_data=rr_for_access,
                anc_name="rainfall_rate",
                anc_attrs=var_attrs,
                global_attrs=global_attrs,
                dataroot=dataroot,
            )
        elif grid_type == "ease2":

            write_daily_ancillary_var_netcdf_polar(
                date=date,
                satellite=satellite,
                target_size=footprint_diameter_km,
                grid_type=grid_type,
                pole=pole,
                anc_data=rr_for_access,
                anc_name="rainfall_rate",
                anc_attrs=var_attrs,
                global_attrs=global_attrs,
                dataroot=dataroot,
            )
    else:
        print(f"No processing needed for  {var} on {date}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Interpolate and append IMERG rainfall data to ACCESS output file. "
            "IMERG data is downloaded if required."
        )
    )
    parser.add_argument(
        "access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument("output_root", type=Path, help="Root directory to write files")
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
        "footprint_diameter",
        type=int,
        help="Diameter of resampling footprint (in km). Default=30km",
        nargs="?",
        default=30,
    )

    parser.add_argument(
        "region", help="region to process", choices=["global", "north", "south"]
    )
    parser.add_argument(
        "--overwrite", help="force overwrite if file exists", action="store_true"
    )
    parser.add_argument(
        "--update",
        help="force overwrite if file older than base file",
        action="store_true",
    )

    parser.add_argument(
        "--redo_attrs", help="rewrite attrs if file exists", action="store_true"
    )

    args = parser.parse_args()

    script_name = parser.prog
    # commit = str(subprocess.check_output(["git", "rev-parse", "HEAD"]))

    repo = git.Repo(search_parent_directories=True)
    commit = repo.head.object.hexsha

    args = parser.parse_args()

    access_root: Path = args.access_root
    output_root: Path = args.output_root
    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()
    footprint_diameter_km = args.footprint_diameter
    region = args.region

    overwrite = args.overwrite
    update = args.update

    redo_attrs = args.redo_attrs

    #Initialize resampler
    if region in ["north", "south"]:
        resampler = ResampleIMERG(target_size=footprint_diameter_km, region=region)
    else:
        resampler = None

    date = START_DAY
    while date <= END_DAY:
        print(f"{date}")
        if redo_attrs:
            redo_imerg_rain_rate_attrs_ACCESS(
                current_day=date,
                satellite=satellite,
                dataroot=access_root,
                outputroot=output_root,
                temproot=temp_root,
                footprint_diameter_km=footprint_diameter_km,
                script_name=script_name,
                commit=commit
            )
        else:
            write_imerg_rain_rate_for_ACCESS(
                current_day=date,
                satellite=satellite,
                dataroot=access_root,
                outputroot=output_root,
                temproot=temp_root,
                footprint_diameter_km=footprint_diameter_km,
                region=region,
                overwrite=overwrite,
                update=update,
                script_name=script_name,
                commit=commit,
                resampler=resampler,
            )
        date += datetime.timedelta(days=1)
