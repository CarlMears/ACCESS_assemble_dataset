import argparse
import datetime
import os
import subprocess
from pathlib import Path

import numpy as np
import xarray as xr

from access_io.access_output import write_daily_lf_netcdf
from access_io.access_output_polar import write_daily_lf_netcdf_polar


def add_land_fraction_to_ACCESS_output(
    *,
    date: datetime.date,
    satellite: str,
    ksat: str,
    target_size: int,
    region: str,
    look: int,
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
        land_path = Path("/mnt/l/access/land_water")

    if region == "global":
        grid_type = "equirectangular"
        pole = None
    elif region in ["north", "south"]:
        grid_type = "ease2"
        pole = region
        if region == "north":
            assert "NP" in str(dataroot)

        if region == "south":
            assert "SP" in str(dataroot)

    else:
        raise ValueError(f"region {region} not valid")

    if region in ["north", "south"]:
        land_file = land_path / (
            f"land_fraction_1440_721_{target_size}km."
            f"from_nsidc_3km_mask.{region}.ease25.v5.nc"
        )
        try:
            land_fraction_xr = xr.open_dataset(land_file)
            land_fraction_np = land_fraction_xr["land_fraction"].values
        except FileNotFoundError:
            raise

    elif region == "global":
        if lf_version.lower() == "combined_hansen":
            land_file = (
                land_path / "land_fraction_1440_721_30km.combined_hansen_nsidc.nc"
            )
            land_fraction_xr = xr.open_dataset(land_file)
            land_fraction_np = land_fraction_xr["land_fraction"].values

        elif lf_version.lower() == "modis":
            land_file = land_path / f"resampled.modislandwater.{target_size}km.nc"
            land_fraction_xr = xr.open_dataset(land_file)
            land_fraction_np = land_fraction_xr[
                "resampled_modis_land_water_mask"
            ].values

            # we need the hansen/nsidc dataset to fill in a few spots that are missing in
            # the modis dataset

            land_file = (
                land_path
                / f"land_fraction_1440_721_{target_size}km.combined_hansen_nsidc.nc"
            )
            land_fraction_xr = xr.open_dataset(land_file)
            land_fraction_np_hansen = land_fraction_xr["land_fraction"].values

            # for areas of missing data, use the hansen/nsidc map
            land_fraction_np[~np.isfinite(land_fraction_np)] = land_fraction_np_hansen[
                ~np.isfinite(land_fraction_np)
            ]
            # also use for lats southward of -75.0, where there is a lot of missing modis data
            # and the modis data is a little suspect because of sea ice.  We assume the NSIDC knows
            # what it is doing in these areas.

            land_fraction_np[0:60, :] = land_fraction_np_hansen[0:60, :]

        else:
            raise KeyError(f"Land version {lf_version} not supported")
    else:
        raise ValueError(f"Region {region} not valid")

    if region in ["north", "south"]:
        try:
            write_daily_lf_netcdf_polar(
                date=date,
                satellite=satellite,
                ksat=ksat,
                target_size=target_size,
                look=look,
                pole=pole,
                grid_type=grid_type,
                version=version,
                lf_version=lf_version,
                land_fraction=land_fraction_np,
                dataroot=dataroot,
                overwrite=overwrite,
                script_name=script_name,
                commit=commit,
            )
        except FileNotFoundError:
            print("File Not Found for {date} - skipping day")
    elif region == "global":
        try:
            write_daily_lf_netcdf(
                date=date,
                satellite=satellite,
                target_size=target_size,
                look=look,
                ksat=ksat,
                version=version,
                lf_version=lf_version,
                land_fraction=land_fraction_np,
                dataroot=dataroot,
                overwrite=overwrite,
                script_name=script_name,
                commit=commit,
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
        "--output_root", type=Path, help="Root directory to write output data"
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
    parser.add_argument(
        "--sensor", choices=["amsr2", "ssmi", "smap"], help="Microwave sensor to use"
    )
    parser.add_argument("--ksat", choices=["13","15"], help="Satellite Number for SSMI")
    parser.add_argument("--target_size")
    parser.add_argument(
        "--look", choices=["0", "1"], default="0", help="Look direction"
    )
    parser.add_argument("--version")
    parser.add_argument("--region", help="region to process")
    parser.add_argument(
        "--lf_version",
        choices=["modis", "combined_hansen", "NSIDC"],
        default="modis",
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

    output_root: Path = args.output_root
    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()
    target_size = int(args.target_size)
    look = int(args.look)

    day_to_do = START_DAY
    while day_to_do <= END_DAY:
        print(f"{day_to_do}")
        add_land_fraction_to_ACCESS_output(
            date=day_to_do,
            satellite=satellite,
            ksat=args.ksat,
            target_size=target_size,
            look=look,
            region=args.region,
            version=args.version,
            lf_version=args.lf_version,
            dataroot=output_root,
            overwrite=args.overwrite,
            script_name=script_name,
            commit=commit,
        )
        day_to_do += datetime.timedelta(days=1)
