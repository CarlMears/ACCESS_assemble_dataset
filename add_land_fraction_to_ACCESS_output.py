import datetime
import os
from pathlib import Path

import xarray as xr

from access_io.access_output import append_const_var_to_daily_tb_netcdf


def add_land_fraction_to_ACCESS_output(
    *, date: datetime.date, satellite: str, dataroot: Path, overwrite: bool
) -> None:

    land_path = land_path = Path("L:/access/land_water")
    combined_land_file = (
        land_path / "land_fraction_1440_721_30km.combined_hansen_nsidc.nc"
    )
    land_fraction_xr = xr.open_dataset(combined_land_file)
    land_fraction_np = land_fraction_xr["land_fraction"].values

    append_const_var_to_daily_tb_netcdf(
        date=date,
        satellite=satellite,
        var=land_fraction_np,
        var_name="land_area_fraction",
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
    import calendar

    year = 2021
    satellite = "AMSR2"
    verbose = True
    if os.name == "nt":
        dataroot = Path("L:/access/amsr2_daily_test")
    elif os.name == "posix":
        dataroot = Path("/mnt/ops1p-ren/l/access/amsr2_daily_test")

    for month in range(7, 8):
        for day in range(1, calendar.monthrange(year, month)[1] + 1):
            date = datetime.date(year, month, day)
            print(f"{date}")
            add_land_fraction_to_ACCESS_output(
                date=date,
                satellite=satellite,
                dataroot=dataroot,
                overwrite=True,
            )
