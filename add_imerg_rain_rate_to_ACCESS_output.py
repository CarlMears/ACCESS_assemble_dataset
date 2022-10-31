"""Add IMERG rain rates to an existing daily ACCESS data file."""

import datetime
from pathlib import Path
import numpy as np
from netCDF4 import Dataset as netcdf_dataset

from access_io.access_output import (
    write_daily_ancillary_var_netcdf,
    get_access_output_filename_daily_folder,
)

from access_io.access_attr_define import common_global_attributes_access
from access_io.access_attr_define import anc_var_attributes_access
from imerg_request.imerg_requests import imerg_half_hourly_request
from resampling_utils.imerg_resampling_routines import resample_imerg_day

import git


def write_imerg_rain_rate_for_ACCESS(
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

    base_filename = get_access_output_filename_daily_folder(
        current_day, satellite, footprint_diameter_km, dataroot, "resamp_tbs"
    )
    imerge_filename_final = get_access_output_filename_daily_folder(
        current_day, satellite, footprint_diameter_km, outputroot, "rain_rate_imerge"
    )

    if not base_filename.is_file():
        print(f"base file for {current_day} does not exist, skipping")
        return

    if imerge_filename_final.is_file():
        if not overwrite:
            print(f"imerge_filename for {current_day} exists, skipping to next day")
            return
        else:
            imerge_filename_final.unlink()

    try:
        # read in base file and extract dimensions and metadata
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
        target_path=temproot / "imerg",
    )

    rr_for_access = np.roll(rr_for_access, 720, axis=1)

    version = "v01r00"
    rr_attrs = anc_var_attributes_access(satellite, "rain_rate_imerg", version=version)
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
        f"area: weighted average over {footprint_diameter_km}km gaussian footprint"
    )
    global_attrs["corresponding_resampled_tb_file"] = base_filename.name
    global_attrs["commit"] = commit
    global_attrs["script_name"] = script_name

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
        "--overwrite", help="force overwrite if file exists", action="store_true"
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

    date = START_DAY
    while date <= END_DAY:
        print(f"{date}")

        write_imerg_rain_rate_for_ACCESS(
            current_day=date,
            satellite=satellite,
            dataroot=access_root,
            outputroot=output_root,
            temproot=temp_root,
            footprint_diameter_km=footprint_diameter_km,
            script_name=script_name,
            commit=commit,
        )
        date += datetime.timedelta(days=1)
