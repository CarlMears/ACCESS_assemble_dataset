"""Fix Time variable in ACCESS output files"""

import datetime
import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset as netcdf_dataset

from access_io.access_output import get_access_output_filename_daily_folder


def fix_time_in_ACCESS_output(
    *, current_day: datetime.date, satellite: str, dataroot: Path
) -> None:
    filename = get_access_output_filename_daily_folder(
        current_day, satellite, 0, dataroot, "UNKNOWN"
    )
    try:
        with netcdf_dataset(filename, "r+") as root_grp:
            time_sec_midnight = root_grp.variables["second_since_midnight"]
            time_np = time_sec_midnight[:, :, :]
            t = np.array([(current_day - datetime.date(1900, 1, 1)).total_seconds()])
            time_np2 = (time_np + t).astype(np.int64).filled(-999999)

            if "time" in root_grp.variables:
                print('Variable "time" already exists')
                time = root_grp.variables["time"]
            else:
                time = root_grp.createVariable(
                    "time",
                    "i8",
                    (
                        "latitude",
                        "longitude",
                        "hours",
                    ),
                    zlib=True,
                    fill_value=-999999,
                )

            time[:, :, :] = time_np2
            time.standard_name = "time"
            time.long_name = "time of satellite observation"
            time.units = "seconds since 1900-01-01 00:00:00.0"
            time.missing = -999999
            time.valid_range = 0
            time.valid_max = 200 * 366 * 24 * 3600
            time.coordinates = "latitude longitude hours"
            print(f"Wrote time for {current_day}")

    except FileNotFoundError:
        print(f"File: {filename} not found, skipping")
        return


if __name__ == "__main__":
    import calendar

    satellite = "AMSR2"
    verbose = True
    if os.name == "nt":
        dataroot = Path("L:/access/amsr2_out")
    elif os.name == "posix":
        dataroot = Path("/mnt/l/access/amsr2_out")

    for year in range(2012, 2022):
        for month in range(1, 13):
            for day in range(1, calendar.monthrange(year, month)[1] + 1):
                date = datetime.date(year, month, day)
                fix_time_in_ACCESS_output(
                    current_day=date, satellite=satellite, dataroot=dataroot
                )
