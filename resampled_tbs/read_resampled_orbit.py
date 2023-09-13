import os
from pathlib import Path
from typing import Any, Tuple, Union

import xarray as xr

if os.name == "nt":
    ACCESS_ROOT = Path("L:/access")
elif os.name == "posix":
    ACCESS_ROOT = Path("/mnt/ops1p-ren/l/access")

IMPLEMENTED_SATELLITES = ["amsr2","ssmi","smap"]

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
    ksat: str = "13"
):
    if grid_type == "equirectangular":
        folder = f"GL_{target_size}"
    elif grid_type == "ease2":
        if pole == "north":
            folder = f"NP_{target_size}"
        elif pole == "south":
            folder = f"SP_{target_size}"
        else:
            raise ValueError(f"Pole: {pole} not valid")
    else:
        raise ValueError(f"grid_type: {grid_type} is not valid")
    orbit_lower, orbit_upper = get_AMSR2_orbit_range(orbit)

    if satellite.lower() == "ssmi": 
        folder = f"f{ksat}/" + folder 
        orbit_dir = (
            dataroot
            / f"{satellite}_tb_orbits"
            / folder
            / f"r{orbit_lower:06d}_{orbit_upper:06d}"
        )
        orbit_dir_time = (
            dataroot / f"{satellite}_tb_orbits" / f"r{orbit_lower:06d}_{orbit_upper:06d}"
        )
    elif satellite.lower() == "smap":
        raise ValueError("SMAP not implemented")
    else:

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
            if satellite.lower() == "ssmi":
                filename = orbit_dir_time / f"f{ksat}_{orbit:06d}.time.nc"
            else:
                filename = orbit_dir_time / f"r{orbit:05d}.time.nc"
        else:
            if satellite.lower() == "ssmi":
                filename = (
                    orbit_dir / f"f{ksat}_{orbit:06d}.grid_tb.{channel_str}.{target_size:03d}km.nc"
                )
            else:
                filename = (
                    orbit_dir / f"r{orbit:05d}.grid_tb.{channel_str}.{target_size:03d}km.nc"
                )
        return filename
    elif grid_type == "ease2":
        if pole in ["north", "south"]:
            if channel_str == "time":
                if satellite.lower() == "ssmi":
                    filename = orbit_dir_time / f"f{ksat}_{orbit:06d}.polar_grid_time.{pole}.nc"
                else:
                    filename = orbit_dir_time / f"r{orbit:05d}.polar_grid_time.{pole}.nc"
            else:
                if satellite.lower() == "ssmi":
                    filename = (
                        orbit_dir / f"f{ksat}_{orbit:06d}.polar_grid_tb.{pole}."
                        f"{channel_str}.{target_size:03d}km.nc"
                    )
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


def read_resampled_tbs(
    *,
    satellite: str,
    ksat: str = "13",
    channel: Union[str, int],
    look: int,
    target_size: int,
    orbit: int,
    grid_type: str = "equirectangular",
    pole: str = "north",
    dataroot: Path = ACCESS_ROOT,
    verbose: bool = False,
    file_name_dict: dict = None,
) -> Tuple[Any, Path]:
    """Using xarray here to take advantage of lazy reading"""

    if satellite.lower() not in IMPLEMENTED_SATELLITES:
        raise ValueError(f"Satellite {satellite} is not implemented")

    if satellite.lower() == "smap":
        filename=file_name_dict[orbit]

        smap_var_dict = {
            'time': 'resamp_times',
            1: 'resamp_tbs',
            2: 'resamp_tbs',
            3: 'resamp_tbs',
            4: 'resamp_tbs'
        }
        
        ds = xr.open_dataset(filename)
        z = ds[smap_var_dict[channel]].values
        if channel in [1,2,3,4]:
            return z[:, :, look, channel-1], filename
        else:
            return z[:,:,look], filename
    else:
        filename = get_resampled_file_name(
            satellite=satellite,
            ksat=ksat,
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
