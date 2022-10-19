"""Add IMERG rain rates to an existing daily ACCESS data file."""

import datetime
from pathlib import Path
import os

import numpy as np
from netCDF4 import Dataset as netcdf_dataset

from access_io.access_output import (
    write_daily_ancillary_var_netcdf,
    get_access_output_filename_daily_folder,  
)
from access_io.access_output import set_or_create_attr
from access_io.access_attr_define import common_global_attributes_access
from access_io.access_attr_define import anc_var_attributes_access
from imerg_request.imerg_requests import imerg_half_hourly_request
from resampling_utils.imerg_resampling_routines import resample_imerg_day


def write_imerg_rain_rate_for_ACCESS(
    *,
    current_day: datetime.date,
    satellite: str,
    dataroot: Path,
    temproot: Path,
    footprint_diameter_km: int,
    force_overwrite: bool = False,
) -> None:

    base_filename = get_access_output_filename_daily_folder(
        current_day, satellite, footprint_diameter_km, dataroot, "resamp_tbs"
    )
    imerge_filename_final = get_access_output_filename_daily_folder(
        current_day, satellite, footprint_diameter_km, dataroot, "rain_rate_imerge"
    )

    if not base_filename.is_file():
        print(f"base file for {current_day} does not exist, skipping")
        return

    if imerge_filename_final.is_file():
        if not force_overwrite:
            print(f"imerge_filename for {current_day} exists, skipping to next day")
            return

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
    # write the results to the existing output file
    today = datetime.date.today()

    version = "v00r00"
    rr_attrs = anc_var_attributes_access(satellite,"imerg_rr",version=version)
    global_attrs = common_global_attributes_access(date,satellite,footprint_diameter_km,version=version)
    global_attrs.update(rr_attrs['global'])
    var_attrs = rr_attrs['var']

    global_attrs['source'] = (
            "Huffman, G.J., E.F. Stocker, D.T. Bolvin, E.J. Nelkin, Jackson Tan "
            "(2019), GPM IMERG Final Precipitation L3 Half Hourly "
            "0.1 degree x 0.1 degree V06, Greenbelt, MD, Goddard Earth Sciences "
            "Data and Information Services Center (GES DISC), "
            f"Accessed: {today.strftime('%m/%d/%Y')}, 10.5067/GPM/IMERG/3B-HH/06"
            )
    global_attrs['cell_method']=(
            "time: closest 30-min IMERG file; "
            f"area: weighted average over {footprint_diameter_km}km footprint"
            )

    write_daily_ancillary_var_netcdf(
        date=date,
        satellite=satellite,
        target_size=footprint_diameter_km,
        anc_data = rr_for_access,
        anc_name = "rainfall_rate",
        anc_attrs = var_attrs,
        global_attrs = global_attrs,
        dataroot = dataroot
    )

    # append_var_to_daily_tb_netcdf(
    #     date=date,
    #     satellite=satellite,
    #     var=rr_for_access,
    #     var_name="rainfall_rate",
    #     standard_name="rainfall_rate",
    #     long_name="rainfall rates from 30-minute IMERG",
    #     valid_min=0.0,
    #     valid_max=50.0,
    #     units="mm/hr",
    #     v_fill=-999.0,
    #     dataroot=dataroot,
    #     overwrite=True,
    # )


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
    parser.add_argument("footprint_diameter", type=int, help="Diameter of resampling footprint (in km). Default=30km", nargs='?', default=30)

    args = parser.parse_args()

    access_root: Path = args.access_root
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
            temproot=temp_root,
            footprint_diameter_km=footprint_diameter_km,
        )
        date += datetime.timedelta(days=1)
