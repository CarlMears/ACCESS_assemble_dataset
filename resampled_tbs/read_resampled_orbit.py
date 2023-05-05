import os
from pathlib import Path
from typing import Any, Tuple, Union

import xarray as xr

if os.name == "nt":
    ACCESS_ROOT = Path("L:/access")
elif os.name == "posix":
    ACCESS_ROOT = Path("/mnt/ops1p-ren/l/access")

IMPLEMENTED_SATELLITES = ["amsr2"]

AVAILABLE_CHANNELS = [
    "time",
    "6V",
    "6H",
    "7V",
    "7H",
    "11V",
    "11H",
    "19V",
    "19H",
    "24V",
    "24H",
    "37V",
    "37H",
    "89V",
    "89H",
]


def get_resampled_file_name(
    satellite: str,
    channel: Union[str, int],
    target_size: int,
    orbit: int,
    grid_type: str = "equirectangular",
    pole: str = "north",
    dataroot: Path = ACCESS_ROOT,
):
    if grid_type == "equirectangular":
        folder = f"GL_{target_size}"
    elif grid_type == "ease2":
        if pole == "north":
            folder = f"NP_{target_size}"
        elif pole == "south":
            folder = f"NP_{target_size}"
        else:
            raise ValueError(f"Pole: {pole} not valid")
    else:
        raise ValueError(f"grid_type: {grid_type} is not valid")
    orbit_lower, orbit_upper = get_AMSR2_orbit_range(orbit)
    orbit_dir = (
        dataroot
        / f"{satellite}_tb_orbits"
        / folder
        / f"r{orbit_lower:05d}_{orbit_upper:05d}"
    )
    orbit_dir_time = (
        dataroot / f"{satellite}_tb_orbits" / f"r{orbit_lower:05d}_{orbit_upper:05d}"
    )

    if isinstance(channel, int):
        if (channel >= 0) and (channel <= 99):
            channel_str = f"ch{channel:02d}"
        else:
            raise ValueError("Channel out of range")
    elif isinstance(channel, str):
        channel_str = channel
    else:
        raise ValueError("Channel type not valid")

    if grid_type == "equirectangular":
        if channel_str == "time":
            filename = orbit_dir_time / f"r{orbit:05d}.time.nc"
        else:
            filename = (
                orbit_dir / f"r{orbit:05d}.grid_tb.{channel_str}.{target_size:03d}km.nc"
            )
        return filename
    elif grid_type == "ease2":
        if pole in ["north", "south"]:
            if channel_str == "time":
                filename = orbit_dir_time / f"r{orbit:05d}.polar_grid_time.{pole}.nc"
            else:
                filename = (
                    orbit_dir / f"r{orbit:05d}.polar_grid_tb.{pole}."
                    f"{channel_str}.{target_size:03d}km.nc"
                )
        else:
            raise ValueError(f"Pole {pole} must be north or south")
        return filename
    else:
        raise ValueError(f"Grid type {grid_type} not valid")


def get_AMSR2_orbit_range(orbit: int) -> Tuple[int, int]:
    """Return the lower/upper bounds to an orbit.

    >>> get_orbit_range(1)
    (1, 5000)
    >>> get_orbit_range(5000)
    (1, 5000)
    >>> get_orbit_range(5001)
    (5001, 10000)
    """
    BIN_WIDTH = 5000
    j = int((orbit - 1) / BIN_WIDTH)
    orbit_lower = 1 + j * BIN_WIDTH
    orbit_upper = (j + 1) * BIN_WIDTH
    return orbit_lower, orbit_upper


def read_AMSR2_resampled_tbs(
    *,
    satellite: str,
    channel: Union[str, int],
    target_size: int,
    orbit: int,
    grid_type: str = "equirectangular",
    pole: str = "north",
    dataroot: Path = ACCESS_ROOT,
    verbose: bool = False,
) -> Tuple[Any, Path]:
    """Using xarray here to take advantage of lazy reading"""

    if satellite.lower() not in IMPLEMENTED_SATELLITES:
        raise ValueError(f"Satellite {satellite} is not implemented")

    filename = get_resampled_file_name(
        satellite=satellite,
        channel=channel,
        target_size=target_size,
        orbit=orbit,
        grid_type=grid_type,
        pole=pole,
        dataroot=dataroot,
    )
    if verbose:
        print(filename)
    ds = xr.open_dataset(filename)
    return ds.Data.values, filename
