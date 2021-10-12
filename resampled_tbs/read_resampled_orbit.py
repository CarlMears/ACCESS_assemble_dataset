import os
import typing
from pathlib import Path

import xarray as xr

if os.name == "nt":
    ACCESS_ROOT = "L:/access"
elif os.name == "posix":
    ACCESS_ROOT = "/mnt/ops1p-ren/l/access"

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


def get_orbit_range(orbit: int) -> tuple[int, int]:
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


def read_resampled_tbs(
    *,
    satellite: str,
    channel: typing.Union[str, int],
    orbit: int,
    dataroot: Path = ACCESS_ROOT,
    verbose: bool = False,
) -> tuple[typing.Any, Path]:
    if satellite not in IMPLEMENTED_SATELLITES:
        raise ValueError(f"Satellite {satellite} is not implemented")
    if isinstance(channel, int):
        if (channel < 1) or (channel > 14):
            raise ValueError(f"channel {channel} is out of range")
        channel_str = AVAILABLE_CHANNELS[channel]
    else:
        if channel in AVAILABLE_CHANNELS:
            channel_str = channel
        else:
            raise ValueError(f"Channel {channel} not valid")

    orbit_lower, orbit_upper = get_orbit_range(orbit)
    orbit_dir = dataroot.joinpath(
        f"{satellite}_tb_orbits", f"r{orbit_lower:05d}_{orbit_upper:05d}"
    )
    if channel_str == "time":
        filename = orbit_dir / f"r{orbit:05d}.time.nc"
    else:
        filename = orbit_dir / f"r{orbit:05d}.gridded_tbs.{channel_str}.nc"
    if verbose:
        print(filename)
    ds = xr.open_dataset(filename)
    return ds.Data.values, filename
