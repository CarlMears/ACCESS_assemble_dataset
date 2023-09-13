from contextlib import suppress
import datetime
import os
from pathlib import Path
from typing import Optional, Sequence, Any
from netCDF4 import Variable

import numpy as np
from numpy.typing import ArrayLike

from rss_lock.locked_dataset import LockedDataset

from access_io.access_attr_define import (
    common_global_attributes_access,
    resamp_tb_attributes_access,
    coord_attributes_access,
    # atm_pars_era5_attributes_access,
    anc_var_attributes_access,
)

from polar_grids import NSIDC_ease2_grids
from access_io.access_output import get_access_output_filename_daily_folder

ease2_grid_25km_north = NSIDC_ease2_grids(pole="north", resolution="25km")
ease2_grid_25km_south = NSIDC_ease2_grids(pole="south", resolution="25km")


if os.name == "nt":
    ACCESS_ROOT = Path("L:/access")
elif os.name == "posix":
    ACCESS_ROOT = Path("/mnt/ops1p-ren/l/access")

IMPLEMENTED_SATELLITES = ["amsr2"]

USED_CHANNELS = {
    "amsr2": [
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
    ]
}
AVAILABLE_CHANNELS = {
    "amsr2": [
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
}

NUM_CHANNELS = 14  # all possible AMSR2 channels
NUM_HOURS = 24
NUM_X_POLE = 720
NUM_Y_POLE = 720


class OkToSkipDay(Exception):
    pass


def set_all_attrs_polar(var: Variable, attrs: dict[str, Any]) -> None:
    """
    Set attributes for a NetCDF variable in a polar coordinate system.

    Args:
        var (Variable): The NetCDF variable to set attributes for.
        attrs (dict[str, Any]): The attributes to set.

    Returns:
        None
    """

    for name, value in attrs.items():
        if (name != "_FillValue") and (name != "global") and (value is not None):
            var.setncattr(name, value)


def write_daily_lf_netcdf_polar(
    *,
    date: datetime.date,
    satellite: str,
    ksat: str,
    target_size: int,
    pole: str = None,
    grid_type: str = None,
    version: str,
    lf_version: str,
    land_fraction: ArrayLike,
    dataroot: Path = ACCESS_ROOT,
    overwrite: bool,
    script_name: str,
    commit: str,
) -> None:
    """
    Write daily land fraction data to a NetCDF file in a polar coordinate system.

    Args:
        date (datetime.date): The date of the data.
        satellite (str): The satellite name.
        target_size (int): The target size.
        pole (str, optional): The pole ('north' or 'south'). Defaults to None.
        grid_type (str, optional): The grid type. Defaults to None.
        version (str): The version string.
        lf_version (str): The land fraction version string.
        land_fraction (ArrayLike): The land fraction data.
        dataroot (Path, optional): The data root path. Defaults to ACCESS_ROOT.
        overwrite (bool): Whether to overwrite existing files.
        script_name (str): The name of the script.
        commit (str): The commit string.

    Returns:
        None
    """

    if pole in ["north", "south"]:
        filename = get_access_output_filename_daily_folder(
            date,
            satellite,
            target_size,
            dataroot,
            "resamp_tbs",
            grid_type="ease2",
            pole=pole,
            ksat=ksat
        )
        if pole == "north":
            lats = ease2_grid_25km_north.latitude
            lons = ease2_grid_25km_north.longitude
            # crs_attrs = ease2_grid_25km_north.crs
        elif pole == "south":
            lats = ease2_grid_25km_south.latitude
            lons = ease2_grid_25km_south.longitude
            # crs_attrs = ease2_grid_25km_south.crs
        else:
            raise ValueError(f"Pole {pole} is not valid")

    else:
        raise ValueError(f"pole = {pole} is not valid")

    lf_fill = -999.0
    lf_string = f"land_frac_{lf_version}"
    filename = get_access_output_filename_daily_folder(
        date, satellite, target_size, dataroot, lf_string, grid_type, pole, ksat=ksat
    )

    if filename.is_file() and not overwrite:
        print(f"daily file for {date} exists... skipping")
        return []
    else:
        with suppress(FileNotFoundError):
            filename.unlink()

    os.makedirs(filename.parent, exist_ok=True)

    lf_attrs = anc_var_attributes_access(
        satellite, "land_fraction_" + lf_version, version=version
    )
    global_attrs = common_global_attributes_access(
        date, satellite, target_size, version=version
    )

    global_attrs.update(lf_attrs["global"])
    var_attrs = lf_attrs["var"]
    global_attrs[
        "cell_method"
    ] = f"area: resampled to {target_size}km Gaussian footprint"

    global_attrs["script_name"] = script_name
    global_attrs["commit"] = commit

    # with netcdf_dataset(filename, "w", format="NETCDF4") as nc_out:
    os.makedirs(filename.parent, exist_ok=True)
    with LockedDataset(filename, "w", 60) as nc_out:
        set_all_attrs_polar(nc_out, global_attrs)

        nc_out.createDimension("x", NUM_X_POLE)
        nc_out.createDimension("y", NUM_Y_POLE)

        x = nc_out.createVariable("x", "f4", ("x",), fill_value=-999.0)
        y = nc_out.createVariable("y", "f4", ("y",), fill_value=-999.0)

        lat_attrs = coord_attributes_access("latitude", dtype=np.float32)
        latitude = nc_out.createVariable(
            "latitude",
            "f4",
            (
                "x",
                "y",
            ),
            zlib=True,
            fill_value=lat_attrs["_FillValue"],
            least_significant_digit=2,
        )

        lon_attrs = coord_attributes_access("longitude", dtype=np.float32)
        lon_attrs["valid_min"] = -180.0
        lon_attrs["valid_max"] = 180.0
        longitude = nc_out.createVariable(
            "longitude",
            "f4",
            (
                "x",
                "y",
            ),
            zlib=True,
            fill_value=lon_attrs["_FillValue"],
            least_significant_digit=2,
        )

        lf = nc_out.createVariable(
            "land_fraction",
            "f4",
            ("x", "y"),
            zlib=True,
            fill_value=var_attrs["_FillValue"],
            least_significant_digit=3,
        )

        x_attrs = coord_attributes_access("x", dtype=np.float64)
        y_attrs = coord_attributes_access("y", dtype=np.float64)

        set_all_attrs_polar(x, x_attrs)
        set_all_attrs_polar(y, y_attrs)
        set_all_attrs_polar(latitude, lat_attrs)
        set_all_attrs_polar(longitude, lon_attrs)

        # write the coordinate data
        x_values = -8987500.0 + 25000.0 * np.arange(0, 720).astype(np.float64)
        y_values = -8987500.0 + 25000.0 * np.arange(0, 720).astype(np.float64)
        y_values = np.flip(y_values)
        x[:] = x_values
        y[:] = y_values
        latitude[:, :] = np.transpose(lats)
        longitude[:, :] = np.transpose(lons)

        lf_to_put = np.transpose(
            np.nan_to_num(
                land_fraction, nan=lf_fill, posinf=lf_fill, neginf=lf_fill
            ).astype(np.float32)
        )

        lf[:, :] = lf_to_put


def write_daily_tb_netcdf_polar(
    *,
    date: datetime.date,
    satellite: str,
    ksat: str,
    target_size: int,
    pole: str = None,
    version: str,
    tb_array_by_hour: ArrayLike,
    time_array_by_hour: ArrayLike,
    dataroot: Path = ACCESS_ROOT,
    freq_list: ArrayLike,
    file_list: Optional[Sequence[Path]],
    script_name: str = "unavailable",
    commit: str = "unavailable",
) -> None:
    """
    Writes daily brightness temperature (TB) data to a NetCDF file
    in a polar grid format.

    Parameters:
        date (datetime.date): Date of the TB data.
        satellite (str): Satellite name.
        target_size (int): Target size of the resampled grid in kilometers.
        pole (str, optional): Pole direction. Valid values are "north" or "south".
                                Defaults to None.
        version (str): Version of the TB data.
        tb_array_by_hour (ArrayLike): Array of TB values by hour.
        time_array_by_hour (ArrayLike): Array of time values by hour.
        dataroot (Path, optional): Root directory for data. Defaults to ACCESS_ROOT.
        freq_list (ArrayLike): List of frequency values.
        file_list (Optional[Sequence[Path]]): List of file paths. Defaults to None.
        script_name (str, optional): Name of the script. Defaults to "unavailable".
        commit (str, optional): Commit ID. Defaults to "unavailable".

    Returns:
        None
    """
    if pole in ["north", "south"]:
        filename = get_access_output_filename_daily_folder(
            date,
            satellite,
            target_size,
            dataroot,
            "resamp_tbs",
            grid_type="ease2",
            pole=pole,
            ksat=ksat
        )
        if pole == "north":
            lats = ease2_grid_25km_north.latitude
            lons = ease2_grid_25km_north.longitude
            crs_attrs = ease2_grid_25km_north.crs
        elif pole == "south":
            lats = ease2_grid_25km_south.latitude
            lons = ease2_grid_25km_south.longitude
            crs_attrs = ease2_grid_25km_south.crs
        else:
            raise ValueError(f"pole = {pole} is not valid")
    else:
        raise ValueError(f"pole = {pole} is not valid")

    os.makedirs(filename.parent, exist_ok=True)

    num_freq = len(freq_list)
    num_pol = 2

    # with netcdf_dataset(filename, "w", format="NETCDF4") as nc_out:
    with LockedDataset(filename, "w", 60) as nc_out:
        # set the global_attributes

        glb_attrs = common_global_attributes_access(
            date, satellite, target_size, version, tb_array_by_hour.dtype
        )
        tb_attrs = resamp_tb_attributes_access(
            satellite, version, tb_array_by_hour.dtype
        )

        glb_attrs.update(tb_attrs["global"])
        tb_attrs = tb_attrs["var"]

        glb_attrs["spatial_resolution"] = f"{target_size} km X {target_size} km"
        glb_attrs[
            "cell_method"
        ] = f"area: resampled to {target_size}km Gaussian footprint"
        glb_attrs["date_accessed"] = f"{datetime.datetime.today()}"
        glb_attrs["date_created"] = f"{datetime.datetime.today()}"

        glb_attrs["id"] = filename.name
        if file_list is not None:
            source_string = ", ".join(str(p) for p in file_list)
            glb_attrs["source"] = source_string

        glb_attrs["script_name"] = script_name
        glb_attrs["commit"] = commit
        set_all_attrs_polar(nc_out, glb_attrs)

        lat_attrs = coord_attributes_access("latitude", np.float32)
        lon_attrs = coord_attributes_access("longitude", np.float32)

        # create Dimensions
        nc_out.createDimension("x", NUM_X_POLE)
        nc_out.createDimension("y", NUM_Y_POLE)
        nc_out.createDimension("hours", NUM_HOURS)
        nc_out.createDimension("freq", num_freq)
        nc_out.createDimension("pol", num_pol)

        x = nc_out.createVariable("x", "f4", ("x",), fill_value=-999.0)
        y = nc_out.createVariable("y", "f4", ("y",), fill_value=-999.0)
        hours = nc_out.createVariable("hours", "i4", ("hours",))
        freq = nc_out.createVariable("freq", "f4", ("freq",))
        pol = nc_out.createVariable("pol", "i4", ("pol",))

        lat_attrs = coord_attributes_access("latitude", dtype=np.float32)
        latitude = nc_out.createVariable(
            "latitude",
            "f4",
            (
                "x",
                "y",
            ),
            zlib=True,
            fill_value=lat_attrs["_FillValue"],
            least_significant_digit=2,
        )

        lon_attrs = coord_attributes_access("longitude", dtype=np.float32)
        lon_attrs["valid_min"] = -180.0
        lon_attrs["valid_max"] = 180.0
        longitude = nc_out.createVariable(
            "longitude",
            "f4",
            (
                "x",
                "y",
            ),
            zlib=True,
            fill_value=lon_attrs["_FillValue"],
            least_significant_digit=2,
        )

        time = nc_out.createVariable(
            "time",
            "i8",
            (
                "x",
                "y",
                "hours",
            ),
            zlib=True,
            fill_value=-999999,
        )
        tbs = nc_out.createVariable(
            "brightness_temperature",
            "f4",
            (
                "x",
                "y",
                "hours",
                "freq",
                "pol",
            ),
            zlib=True,
            fill_value=tb_attrs["_FillValue"],
            least_significant_digit=2,
        )
        crs = nc_out.createVariable(
            "crs",
            "i4",
            (),
        )

        # set coordinate attributes from .json files

        x_attrs = coord_attributes_access("x", dtype=np.float64)
        y_attrs = coord_attributes_access("y", dtype=np.float64)
        hours_attrs = coord_attributes_access("hours", date=date, dtype=np.int32)
        freq_attrs = coord_attributes_access("frequency", dtype=np.float32)
        pol_attrs = coord_attributes_access("polarization", dtype=np.int32)
        time_attrs = coord_attributes_access("time", dtype=np.int64)

        set_all_attrs_polar(x, x_attrs)
        set_all_attrs_polar(y, y_attrs)
        set_all_attrs_polar(latitude, lat_attrs)
        set_all_attrs_polar(longitude, lon_attrs)
        set_all_attrs_polar(hours, hours_attrs)
        set_all_attrs_polar(freq, freq_attrs)
        set_all_attrs_polar(pol, pol_attrs)
        set_all_attrs_polar(time, time_attrs)
        set_all_attrs_polar(crs, crs_attrs)

        set_all_attrs_polar(tbs, tb_attrs)

        # write the coordinate data
        x_values = -8987500.0 + 25000.0 * np.arange(0, 720).astype(np.float64)
        y_values = -8987500.0 + 25000.0 * np.arange(0, 720).astype(np.float64)
        y_values = np.flip(y_values)
        x[:] = x_values
        y[:] = y_values
        latitude[:, :] = np.transpose(lats)
        longitude[:, :] = np.transpose(lons)
        hours[:] = np.arange(NUM_HOURS)
        freq[:] = freq_list
        pol[:] = np.array([0, 1], dtype=np.int32)

        # write the time layer
        time_to_put = (
            time_array_by_hour + (date - datetime.date(1900, 1, 1)).total_seconds()
        )
        time_to_put = np.nan_to_num(
            time_to_put, nan=-999998.99, posinf=-999998.99, neginf=-999998.99
        )
        time_to_put = np.floor(time_to_put).astype(np.int64)
        time_to_put = time_to_put.swapaxes(0, 1)
        time[:, :, :] = time_to_put

        fill_val = np.float32(tb_attrs["_FillValue"])
        tbs_to_put = np.nan_to_num(
            tb_array_by_hour, nan=fill_val, posinf=fill_val, neginf=fill_val
        ).astype(np.float32)
        tbs_to_put = tbs_to_put.swapaxes(0, 1)
        tbs[:, :, :, :, :] = tbs_to_put

        crs[:] = 1


def write_daily_ancillary_var_netcdf_polar(
    *,
    date: datetime.date,
    satellite: str,
    target_size: int,
    pole: str,
    grid_type: str,
    anc_data: ArrayLike,
    anc_name: str,
    anc_attrs: dict,
    global_attrs: dict,
    dataroot: Path = ACCESS_ROOT,
) -> None:
    """Write daily ancillary variable NetCDF file in polar coordinates.

    Args:
        date (datetime.date): Date of the data.
        satellite (str): Satellite name.
        target_size (int): Target size.
        pole (str): Pole name ("north" or "south").
        grid_type (str): Grid type.
        anc_data (ArrayLike): Ancillary data.
        anc_name (str): Ancillary variable name.
        anc_attrs (dict): Ancillary variable attributes.
        global_attrs (dict): Global attributes.
        dataroot (Path, optional): Root path of the data. Defaults to ACCESS_ROOT.
    """
    if pole in ["north", "south"]:
        base_filename = get_access_output_filename_daily_folder(
            date,
            satellite,
            target_size,
            dataroot,
            "resamp_tbs",
            grid_type="ease2",
            pole=pole,
        )
        var_filename = get_access_output_filename_daily_folder(
            date,
            satellite.lower(),
            target_size,
            dataroot,
            f"{anc_name}_temp",
            grid_type="ease2",
            pole=pole,
        )
        var_filename_final = get_access_output_filename_daily_folder(
            date,
            satellite.lower(),
            target_size,
            dataroot,
            anc_name,
            grid_type="ease2",
            pole=pole,
        )

        # if pole == "north":
        #     crs_attrs = ease2_grid_25km_north.crs
        # else:
        #     crs_attrs = ease2_grid_25km_south.crs

    else:
        raise ValueError(f"pole = {pole} is not valid")

    with LockedDataset(var_filename, "w", 60) as nc_out:
        with LockedDataset(base_filename, "r", 60) as root_grp:
            # Create the dimensions of the file from the base file
            for name, dim in root_grp.dimensions.items():
                if name in ["hours", "x", "y"]:
                    nc_out.createDimension(
                        name, len(dim) if not dim.isunlimited() else None
                    )

            set_all_attrs_polar(nc_out, global_attrs)

            for var_name in ["hours", "x", "y", "latitude", "longitude"]:
                # Create the time and dimension variables in the output file
                var_in = root_grp[var_name]
                # print(f'Transferring {var_name}')
                nc_out.createVariable(
                    var_name, var_in.dtype, var_in.dimensions, zlib=True
                )

                # Copy the attributes
                nc_out.variables[var_name].setncatts(
                    {a: var_in.getncattr(a) for a in var_in.ncattrs()}
                )

                # Copy the data values
                nc_out.variables[var_name][:] = root_grp.variables[var_name][:]

            # make the ancillary variable with the same dimensions as the
            # time variable in the base file
            # fill value is set from the var_attrs dictionary
            new_var = nc_out.createVariable(
                anc_name,
                np.float32,
                root_grp.variables["time"].dimensions,
                zlib=True,
                fill_value=np.float32(anc_attrs["_FillValue"]),
            )

            for key in anc_attrs.keys():
                if (
                    key != "_FillValue"
                ):  # _FillValue already set by the variable definition
                    new_var.setncattr(key, anc_attrs[key])

            new_var.coordinates = "latitude longitude"
            new_var[:, :, :] = anc_data[:, :, :]

    # everything written -- rename to final file name
    with suppress(FileNotFoundError):
        var_filename_final.unlink()

    var_filename.rename(var_filename_final)


def edit_attrs_daily_ancillary_var_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    target_size: int,
    anc_data: ArrayLike,
    anc_name: str,
    anc_attrs: dict,
    global_attrs: dict,
    dataroot: Path = ACCESS_ROOT,
) -> None:
    """Edit attributes of daily ancillary variable NetCDF file in polar coordinates.

    Args:
        date (datetime.date): Date of the data.
        satellite (str): Satellite name.
        target_size (int): Target size.
        anc_data (ArrayLike): Ancillary data.
        anc_name (str): Ancillary variable name.
        anc_attrs (dict): Ancillary variable attributes.
        global_attrs (dict): Global attributes.
        dataroot (Path, optional): Root path of the data. Defaults to ACCESS_ROOT.
    """

    var_filename_final = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, dataroot, anc_name
    )

    with LockedDataset(var_filename_final, "r+", 60) as nc_out:
        set_all_attrs_polar(nc_out, global_attrs)
        set_all_attrs_polar(nc_out.variables[anc_name], anc_attrs)

        lat_attrs = coord_attributes_access("latitude", dtype=np.float32)
        lon_attrs = coord_attributes_access("longitude", dtype=np.float32)
        hours_attrs = coord_attributes_access("hours", date=date, dtype=np.int32)

        set_all_attrs_polar(nc_out.variables["latitude"], lat_attrs)
        set_all_attrs_polar(nc_out.variables["longitude"], lon_attrs)
        set_all_attrs_polar(nc_out.variables["hours"], hours_attrs)


def write_ocean_emiss_to_daily_ACCESS_polar(
    *,
    ocean_emiss: ArrayLike,
    current_day: datetime.date,
    satellite: str,
    target_size: int,
    grid_type: str,
    pole: str,
    global_attrs: dict,
    var_attrs: dict,
    dataroot: Path,
    outputroot: Path,
    verbose: bool = False,
) -> None:
    """Write ocean emissivity data to daily ACCESS NetCDF file in polar coordinates.

    Args:
        ocean_emiss (ArrayLike): Ocean emissivity data.
        current_day (datetime.date): Current date.
        satellite (str): Satellite name.
        target_size (int): Target size.
        grid_type (str): Grid type.
        pole (str): Pole name ("north" or "south").
        global_attrs (dict): Global attributes.
        var_attrs (dict): Variable attributes.
        dataroot (Path): Root path of the input data.
        outputroot (Path): Root path for the output file.
        verbose (bool, optional): Verbosity flag. Defaults to False.
    """
    # if satellite.lower() == "amsr2":
    #     from satellite_definitions.amsr2 import REF_FREQ

    base_filename = get_access_output_filename_daily_folder(
        current_day, satellite, target_size, dataroot, "resamp_tbs", grid_type, pole
    )
    emiss_filename_final = get_access_output_filename_daily_folder(
        current_day,
        satellite,
        target_size,
        outputroot,
        "ocean_emiss_era5",
        grid_type,
        pole,
    )

    # if pole == "north":
    #     crs_attrs = ease2_grid_25km_north.crs
    # elif pole == "south":
    #     crs_attrs = ease2_grid_25km_south.crs
    # else:
    #     raise ValueError(f"Pole: {pole} is not valid")

    try:
        with LockedDataset(emiss_filename_final, "w", 60) as nc_out:
            with LockedDataset(base_filename, "r", 60) as root_grp:
                # Create the dimensions of the file from the base file
                for name, dim in root_grp.dimensions.items():
                    if name in ["hours", "x", "y", "freq", "pol"]:
                        nc_out.createDimension(
                            name, len(dim) if not dim.isunlimited() else None
                        )
                set_all_attrs_polar(nc_out, global_attrs)

                for var_name in [
                    "hours",
                    "x",
                    "y",
                    "latitude",
                    "longitude",
                    "freq",
                    "pol",
                ]:
                    # Create the time and dimension variables in the output file
                    var_in = root_grp[var_name]
                    # print(f'Transferring {var_name}')
                    nc_out.createVariable(
                        var_name, var_in.dtype, var_in.dimensions, zlib=True
                    )

                    # Copy the attributes
                    nc_out.variables[var_name].setncatts(
                        {a: var_in.getncattr(a) for a in var_in.ncattrs()}
                    )

                    # Copy the data values
                    nc_out.variables[var_name][:] = root_grp.variables[var_name][:]

                # create the new dimensions needed for the emissivity data
                # nc_out.createDimension("freq", len(REF_FREQ))
                # nc_out.createDimension("pol", 2)

                dimensions_out = ("x", "y", "hours", "freq", "pol")

                for varname, long_name, units in [
                    ("emissivity", "ocean surface emissivity", None),
                ]:
                    print(f"starting writing {varname} to {emiss_filename_final}")
                    least_significant_digit = 3
                    nc_out.createVariable(
                        varname,
                        np.float32,
                        dimensions_out,
                        fill_value=var_attrs["_FillValue"],
                        zlib=True,
                        least_significant_digit=least_significant_digit,
                    )
                    set_all_attrs_polar(nc_out, var_attrs)
                    nc_out[varname][:, :, :, :, :] = ocean_emiss

                    print(f"finished writing {varname}")

    except FileNotFoundError:
        raise OkToSkipDay
