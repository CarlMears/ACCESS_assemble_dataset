import argparse
import datetime
import numpy as np
from pathlib import Path
import xarray as xr

from access_io.access_output import append_const_var_to_daily_tb_netcdf


def add_land_fraction_to_ACCESS_output(
    *,
    date: datetime.date,
    satellite: str,
    dataroot: Path,
    overwrite: bool,
    version: str = "modis",
) -> None:

    if version.lower() == "combined_hansen":
        land_path = land_path = Path("L:/access/land_water")
        land_file = land_path / "land_fraction_1440_721_30km.combined_hansen_nsidc.nc"
        land_fraction_xr = xr.open_dataset(land_file)
        land_fraction_np = land_fraction_xr["land_fraction"].values
        var_name = "land_area_fraction_hansen"
    elif version.lower() == "modis":
        land_path = land_path = Path("L:/access/land_water")
        land_file = land_path / "resampled.modislandwater.nc"
        land_fraction_xr = xr.open_dataset(land_file)
        land_fraction_np = land_fraction_xr["resampled_modis_land_water_mask"].values
        var_name = "land_area_fraction_modis"

        # we need the other dataset to fill in a few spots....
        land_path = land_path = Path("L:/access/land_water")
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
        raise KeyError(f"Land version {version} not supported")

    append_const_var_to_daily_tb_netcdf(
        date=date,
        satellite=satellite,
        var=land_fraction_np,
        var_name=var_name,
        standard_name="land_area_fraction",
        long_name="land fraction averaged over gaussian footprint",
        valid_min=0.0,
        valid_max=1.0,
        units="dimensionless",
        v_fill=-999.0,
        dataroot=dataroot,
        overwrite=True,
        verbose=True,
        lock_stale_time=30.0,
    )


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
    parser.add_argument(
        "version", choices=["modis", "combined_hansen"], help="Microwave sensor to use"
    )
    parser.add_argument(
        "--verbose", help="enable more verbose screen output", action="store_true"
    )
    parser.add_argument(
        "--overwrite", help="force overwrite if file exists", action="store_true"
    )

    args = parser.parse_args()

    access_root: Path = args.access_root
    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()

    day_to_do = START_DAY
    while day_to_do <= END_DAY:
        print(f"{day_to_do}")
        add_land_fraction_to_ACCESS_output(
            date=day_to_do,
            satellite=satellite,
            version=args.version,
            dataroot=access_root,
            overwrite=args.overwrite,
        )
        day_to_do += datetime.timedelta(days=1)
