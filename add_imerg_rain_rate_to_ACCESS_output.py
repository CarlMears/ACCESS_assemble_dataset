"""Add IMERG rain rates to an existing daily ACCESS data file."""

import datetime
from pathlib import Path
import os

import numpy as np
from netCDF4 import Dataset as netcdf_dataset

from access_io.access_output import (
    append_var_to_daily_tb_netcdf,
    get_access_output_filename_daily_folder,
)
from imerg_request.imerg_requests import imerg_half_hourly_request
from resampling_utils.imerg_resampling_routines import resample_imerg_day


def write_imerg_rain_rate_for_ACCESS(
    *,
    current_day: datetime.date,
    satellite: str,
    dataroot: Path,
    temproot: Path,
    force_overwrite: bool = False,
) -> None:

    base_filename = get_access_output_filename_daily_folder(
        current_day, satellite, dataroot, "resamp_tbs"
    )
    imerge_filename = get_access_output_filename_daily_folder(
        current_day, satellite, dataroot, "rain_rate_imerge_temp"
    )
    imerge_filename_final = get_access_output_filename_daily_folder(
        current_day, satellite, dataroot, "rain_rate_imerge"
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
        target_path=temproot / "imerg",
    )
    rr_for_access = np.roll(rr_for_access, 720, axis=1)
    # write the results to the existing output file
    # today = datetime.date.today()

    trg = netcdf_dataset(imerge_filename, mode="w")

    with netcdf_dataset(base_filename, "r") as root_grp:
        # Create the dimensions of the file
        for name, dim in root_grp.dimensions.items():
            if name in ["hours", "latitude", "longitude"]:
                trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

        # Copy the global attributes
        trg.setncatts({a: root_grp.getncattr(a) for a in root_grp.ncattrs()})

        for var_name in ["time", "hours", "longitude", "latitude"]:
            # Create the time variables in the output file
            var_in = root_grp[var_name]
            trg.createVariable(var_name, var_in.dtype, var_in.dimensions, zlib=True)

            # Copy the time attributes
            trg.variables[var_name].setncatts(
                {a: var_in.getncattr(a) for a in var_in.ncattrs()}
            )

            # Copy the time values
            trg.variables[var_name][:] = root_grp.variables[var_name][:]

        # make the rain rate variable with the same dimensions as the time variable in the base file
        rr = trg.createVariable(
            "rainfall_rate", np.float32, trg.variables["time"].dimensions, zlib=True
        )
        rr.var_name = ("rainfall_rate",)
        rr.standard_name = ("rainfall_rate",)
        rr.long_name = ("rainfall rates from 30-minute IMERG",)
        rr.valid_min = (0.0,)
        rr.valid_max = (50.0,)
        rr.units = ("mm/hr",)
        rr.source = (
            (
                "Huffman, G.J., E.F. Stocker, D.T. Bolvin, E.J. Nelkin, Jackson Tan "
                "(2019), GPM IMERG Final Precipitation L3 Half Hourly "
                "0.1 degree x 0.1 degree V06, Greenbelt, MD, Goddard Earth Sciences "
                "Data and Information Services Center (GES DISC), "
                f"Accessed: {datetime.datetime.utcfromtimestamp(mod_time)}"
                "doi: http://doi.org/10.5067/GPM/IMERG/3B-HH/06"
            ),
        )
        rr.cell_method = (
            (
                "time: closest 30-min IMERG file; "
                "area: weighted average over 30km footprint"
            ),
        )
        rr.v_fill = (-999.0,)
        rr.coordinates = "latitude longitude"

        rr[:, :, :] = rr_for_access

        trg.close()
        imerge_filename.rename(imerge_filename_final)


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
    satellite = args.sensor  # .upper()

    date = START_DAY
    while date <= END_DAY:
        print(f"{date}")

        write_imerg_rain_rate_for_ACCESS(
            current_day=date,
            satellite=satellite,
            dataroot=access_root,
            temproot=temp_root,
            force_overwrite=False,
        )
        date += datetime.timedelta(days=1)
