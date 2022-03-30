#!/usr/bin/env python3
"""Download one day for ERA5 data (including the first hour of the next day).

In order to use it, the Copernicus Data Server (cdsapi) needs to be installed in
your version of python.

See https://cds.climate.copernicus.eu/api-how-to

If you use conda: conda install -c conda-forge cdsapi

Then you need to get a key from ECMWF and put in a .cdsapirc file in your {user}
folder.
"""

import datetime
import os
from pathlib import Path
from typing import Sequence
from calendar import monthrange

import cdsapi


def era5_hourly_single_level_request(
    *,
    date: datetime.date,
    variable: str,
    target_path: Path,
    full_day: bool = True,
    full_month: bool = False,
) -> Path:

    c = cdsapi.Client()

    # target = target_path / f"ERA5_Skin_Temperature_{date:%Y_%m}.nc"
    if full_month:
        target = target_path / f"ERA5_Skin_Temperature_{date:%Y_%m}.full.nc"
        times = [f"{h:02d}:00" for h in range(0, 24)]
        day_list = range(1, monthrange(date.year, date.month)[1] + 1)
        days = [f"{day:02d}" for day in day_list]
    elif full_day:
        target = target_path / f"ERA5_Skin_Temperature_{date:%Y_%m_%d}.full.nc"
        times = [f"{h:02d}:00" for h in range(0, 24)]
        days = [f"{date:%d}"]
    else:
        target = target_path / f"ERA5_Skin_Temperature_{date:%Y_%m_%d}.1st_hour.nc"
        times = ["00:00"]
        days = [f"{date:%d}"]

    temp_file = target_path / "temp.nc"

    if target.exists():
        print(f"File: {target} already exists, skipping")
    else:
        print(f"Getting: {target}")
        c.retrieve(
            "reanalysis-era5-single-levels",
            {
                "product_type": "reanalysis",
                "format": "netcdf",
                "grid": "0.25/0.25",
                "variable": "Skin temperature",
                "year": f"{date:%Y}",
                "month": f"{date:%m}",
                "day": days,
                "time": times,
            },
            temp_file,
        )
        temp_file.rename(target)
    return target


def era5_hourly_single_level_request_entire_month(
    *, date: datetime.date, variable: str, target_path: Path
) -> Path:
    c = cdsapi.Client()

    # target = target_path / f"ERA5_Skin_Temperature_{date:%Y_%m}.nc"

    target = target_path / f"ERA5_Skin_Temperature_{date:%Y_%m}.full_month.nc"
    times = [f"{h:02d}:00" for h in range(0, 24)]

    year = int(f"{date:%Y}")
    month = int(f"{date:%m}")
    num_days_in_month = monthrange(year, month)[1]
    days = [f"{d:02d}" for d in range(1, 1 + num_days_in_month)]

    temp_file = target_path / "temp.nc"

    if target.exists():
        print(f"File: {target} already exists, skipping")
    else:
        print(f"Getting: {target}")
        c.retrieve(
            "reanalysis-era5-single-levels",
            {
                "product_type": "reanalysis",
                "format": "netcdf",
                "grid": "0.25/0.25",
                "variable": "Skin temperature",
                "year": f"{date:%Y}",
                "month": f"{date:%m}",
                "day": days,
                "time": times,
            },
            temp_file,
        )
        temp_file.rename(target)
    return target


def era5_hourly_pressure_level_request(
    *,
    date: datetime.date,
    variable: str,
    levels: Sequence[str] = ["775", "875", "975"],
    target_path: Path,
    full_day: bool = True,
) -> Path:
    c = cdsapi.Client()

    if full_day:
        target = target_path / f"ERA5_{variable}_{date:%Y_%m_%d}.full.nc"
        times = [f"{h:02d}:00" for h in range(0, 24)]
    else:
        target = target_path / f"ERA5_{variable}_{date:%Y_%m_%d}.1st_hour.nc"
        times = ["00:00"]

    temp_file = target_path / "temp.nc"

    if target.exists():
        print(f"File: {target} already exists, skipping")
    else:
        print(f"Getting: {target}")
        c.retrieve(
            "reanalysis-era5-pressure-levels",
            {
                "product_type": "reanalysis",
                "format": "netcdf",
                "grid": "0.25/0.25",
                "variable": variable,
                "pressure_level": levels,
                "year": f"{date:%Y}",
                "month": f"{date:%m}",
                "day": f"{date:%d}",
                "time": times,
            },
            temp_file,
        )
        temp_file.rename(target)
    return target


if __name__ == "__main__":
    date = datetime.date(2012, 7, 11)
    variable = "temperature"
    if os.name == "nt":
        target_path = Path("L:/access/_temp")
    elif os.name == "posix":
        target_path = Path("/mnt/ops1p-ren/l/access/_temp")

    era5_hourly_single_level_request_entire_month(
        date=date, variable="skt", target_path=target_path
    )

    levels = [
        "1",
        "2",
        "3",
        "5",
        "7",
        "10",
        "20",
        "30",
        "50",
        "70",
        "100",
        "125",
        "150",
        "175",
        "200",
        "225",
        "250",
        "300",
        "350",
        "400",
        "450",
        "500",
        "550",
        "600",
        "650",
        "700",
        "750",
        "775",
        "800",
        "825",
        "850",
        "875",
        "900",
        "925",
        "950",
        "975",
        "1000",
    ]

    file = era5_hourly_pressure_level_request(
        date=date,
        variable=variable,
        levels=levels,
        target_path=target_path,
        full_day=True,
    )

    print(file)
