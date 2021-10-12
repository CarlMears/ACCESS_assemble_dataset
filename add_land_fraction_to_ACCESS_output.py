import datetime
import os
from pathlib import Path

import numpy as np

from access_io.access_output import append_const_var_to_daily_tb_netcdf
from util.land_fraction import read_land_fraction_1440_720


def add_land_fraction_to_ACCESS_output(
    *, date: datetime.date, satellite: str, dataroot: Path, overwrite: bool
) -> None:
    # This a temporary kludge until the more accurate land mask is available
    lf_1440_720 = read_land_fraction_1440_720()
    lf_1440_721 = np.zeros((721, 1440))
    temp = np.zeros((719, 1440))
    temp = 0.25 * (
        lf_1440_720[0:719, :]
        + lf_1440_720[1:720, :]
        + np.roll(lf_1440_720[0:719, :], 1, axis=1)
        + np.roll(lf_1440_720[1:720, :], 1, axis=1)
    )
    lf_1440_721[1:720, :] = temp
    lf_1440_721[0, :] = lf_1440_721[1, :]
    lf_1440_721[720, :] = lf_1440_721[719, :]

    append_const_var_to_daily_tb_netcdf(
        date=date,
        satellite=satellite,
        var=lf_1440_721,
        var_name="land_area_fraction",
        standard_name="land_area_fraction",
        long_name="land fraction averaged over gaussian footprint",
        valid_min=0.0,
        valid_max=1.0,
        units="dimensionless",
        v_fill=-999.0,
        dataroot=dataroot,
        overwrite=True,
    )


if __name__ == "__main__":
    date = datetime.date(2012, 7, 11)
    satellite = "AMSR2"
    verbose = True
    if os.name == "nt":
        dataroot = Path("L:/access/")
    elif os.name == "posix":
        dataroot = Path("/mnt/ops1p-ren/l/access")

    add_land_fraction_to_ACCESS_output(
        date=date,
        satellite=satellite,
        dataroot=dataroot,
        overwrite=True,
    )
