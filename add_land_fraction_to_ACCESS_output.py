import argparse
import datetime
import numpy as np
import os
from pathlib import Path
import subprocess
import xarray as xr

from access_io.access_output import write_daily_lf_netcdf


def add_land_fraction_to_ACCESS_output(
    *,
    date: datetime.date,
    satellite: str,
    target_size: int,
    version: str,
    dataroot: Path,
    script_name: str,
    commit: str,
    overwrite: bool,
    lf_version: str = "modis",
) -> None:

    if os.name == "nt":
        land_path = Path("L:/access/land_water")
    elif os.name == "posix":
        land_path = Path("/mnt/ops1p-ren/l/access/land_water")

    if lf_version.lower() == "combined_hansen":
        land_file = land_path / "land_fraction_1440_721_30km.combined_hansen_nsidc.nc"
        land_fraction_xr = xr.open_dataset(land_file)
        land_fraction_np = land_fraction_xr["land_fraction"].values
        var_name = "land_area_fraction_hansen"
    elif lf_version.lower() == "modis":
        land_file = land_path / f"resampled.modislandwater.{target_size}km.nc"
        land_fraction_xr = xr.open_dataset(land_file)
        land_fraction_np = land_fraction_xr["resampled_modis_land_water_mask"].values
        var_name = "land_area_fraction_modis"

        # we need the other dataset to fill in a few spots....
        land_file = land_path / "land_fraction_1440_721_30km.combined_hansen_nsidc.nc"
        land_fraction_xr = xr.open_dataset(land_file)
        land_fraction_np_hansen = land_fraction_xr["land_fraction"].values

        # for areas of missing data, use the hansen map
        land_fraction_np[~np.isfinite(land_fraction_np)] = land_fraction_np_hansen[
            ~np.isfinite(land_fraction_np)
        ]
        # also use for lats southward of -75.0
        land_fraction_np[0:60, :] = land_fraction_np_hansen[0:60, :]

    else:
        raise KeyError(f"Land version {lf_version} not supported")

    try:

        write_daily_lf_netcdf(
            date=date,
            satellite=satellite,
            target_size=target_size,
            version=version,
            lf_version=lf_version,
            land_fraction=land_fraction_np,
            dataroot=dataroot,
            overwrite=overwrite,
        )

    except FileNotFoundError:
        print("File Not Found for {date} - skipping day")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Interpolate and append surface temperature to ACCESS output file. "
            "ERA5 data is downloaded if required."
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
    parser.add_argument("target_size")
    parser.add_argument("version")
    parser.add_argument(
        "lf_version",
        choices=["modis", "combined_hansen"],
        help="Microwave sensor to use",
    )
    parser.add_argument(
        "--verbose", help="enable more verbose screen output", action="store_true"
    )
    parser.add_argument(
        "--overwrite", help="force overwrite if file exists", action="store_true"
    )

    args = parser.parse_args()

    script_name = parser.prog
    commit = str(subprocess.check_output(["git", "rev-parse", "HEAD"]))

    access_root: Path = args.access_root
    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()
    target_size = int(args.target_size)

    day_to_do = START_DAY
    while day_to_do <= END_DAY:
        print(f"{day_to_do}")
        add_land_fraction_to_ACCESS_output(
            date=day_to_do,
            satellite=satellite,
            target_size=target_size,
            version=args.version,
            lf_version=args.lf_version,
            dataroot=access_root,
            overwrite=args.overwrite,
            script_name=script_name,
            commit=commit,
        )
        day_to_do += datetime.timedelta(days=1)
