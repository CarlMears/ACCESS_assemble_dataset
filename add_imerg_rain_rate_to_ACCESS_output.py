"""Add IMERG rain rates to an existing daily ACCESS data file."""

import datetime
from pathlib import Path

import numpy as np
from netCDF4 import Dataset as netcdf_dataset

from access_io.access_output import (
    append_var_to_daily_tb_netcdf,
    get_access_output_filename,
)
from imerg_request.imerg_requests import imerg_half_hourly_request
from resampling_utils.imerg_resampling_routines import resample_imerg_day


def add_imerg_rain_rate_to_ACCESS_output(
    *,
    current_day: datetime.date,
    satellite: str,
    dataroot: Path,
    temproot: Path,
    force_overwrite: bool = False,
) -> None:

    filename = get_access_output_filename(current_day, satellite, dataroot)
    try:
        with netcdf_dataset(filename, "r") as root_grp:
            if not force_overwrite:
                try:
                    root_grp.variables["rainfall_rate"][:, :, :].filled(fill_value=-999)
                    print(f"var rainfall rate already exists for {str(current_day)}.")
                    print("skipping to next day")
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
                raise ValueError(f'Error finding "time" in {filename}')
    except FileNotFoundError:
        print(f"File: {filename} not found, skipping")
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
    rr_for_access = resample_imerg_day(
        np.roll(times, 720, axis=1),
        hourly_intervals,
        date,
        target_path=temproot / "imerg",
    )
    rr_for_access = np.roll(rr_for_access, 720, axis=1)
    # write the results to the existing output file
    today = datetime.date.today()
    append_var_to_daily_tb_netcdf(
        date=date,
        satellite=satellite,
        var=rr_for_access,
        var_name="rainfall_rate",
        standard_name="rainfall_rate",
        long_name="rainfall rates from 30-minute IMERG",
        valid_min=0.0,
        valid_max=50.0,
        units="mm/hr",
        source=(
            "Huffman, G.J., E.F. Stocker, D.T. Bolvin, E.J. Nelkin, Jackson Tan "
            "(2019), GPM IMERG Final Precipitation L3 Half Hourly "
            "0.1 degree x 0.1 degree V06, Greenbelt, MD, Goddard Earth Sciences "
            "Data and Information Services Center (GES DISC), "
            f"Accessed: {today.strftime('%m/%d/%Y')}, 10.5067/GPM/IMERG/3B-HH/06"
        ),
        cell_method=(
            "time: closest 30-min IMERG file; "
            "area: weighted average over 30km footprint"
        ),
        v_fill=-999.0,
        dataroot=dataroot,
        overwrite=True,
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

    args = parser.parse_args()

    access_root: Path = args.access_root
    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()

    date = START_DAY
    while date <= END_DAY:
        print(f"{date}")

        add_imerg_rain_rate_to_ACCESS_output(
            current_day=date,
            satellite=satellite,
            dataroot=access_root,
            temproot=temp_root,
        )
        date += datetime.timedelta(days=1)
