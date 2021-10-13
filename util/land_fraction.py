import os
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike

if os.name == "nt":
    LAND_FRACTION_PATH = Path("//ops1p-to-be-renamed/M/Module_Data/Land_Fraction/")
    SEA_ICE_LAND_FRACTION_DIR = Path("l:/sea_ice/land_fraction/")
elif os.name == "posix":
    LAND_FRACTION_PATH = Path("/mnt/ops1p-ren/m/Module_Data/Land_Fraction/")
    SEA_ICE_LAND_FRACTION_DIR = Path("/mnt/ops1p-ren/l/sea_ice/land_fraction/")


def read_land_fraction_1440_720(path: Path = LAND_FRACTION_PATH) -> ArrayLike:
    bin_file = LAND_FRACTION_PATH / "land_fraction_1440x720.dat"

    with bin_file.open(mode="rb") as file:
        land_frac_raw = np.fromfile(file, dtype=np.float32)

    land_frac = np.reshape(land_frac_raw, (720, 1440))
    return land_frac


def read_land_fraction_polar_stereographic(pole: str = "north") -> ArrayLike:
    if pole == "north":
        return read_land_fraction_polar_stereographic_NP()
    elif pole == "south":
        return read_land_fraction_polar_stereographic_SP()
    else:
        raise ValueError('arg "pole" should be "north" or "south"')


def read_land_fraction_polar_stereographic_NP() -> ArrayLike:
    land_fraction_file = (
        SEA_ICE_LAND_FRACTION_DIR
        / "nsidc_polar_stereographic_land_fraction_north_pole.dat"
    )
    lf = np.fromfile(land_fraction_file, dtype="float64")
    lf = np.reshape(lf, (448, 304))
    return lf


def read_land_fraction_polar_stereographic_SP() -> ArrayLike:
    land_fraction_file = (
        SEA_ICE_LAND_FRACTION_DIR
        / "nsidc_polar_stereographic_land_fraction_south_pole.dat"
    )
    lf = np.fromfile(land_fraction_file, dtype="float64")
    lf = np.reshape(lf, (332, 316))
    return lf
