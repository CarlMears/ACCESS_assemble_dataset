import datetime
import os
from contextlib import suppress
from pathlib import Path
from typing import Optional, Sequence, Any, Union, Literal
import numpy as np
from netCDF4 import Variable
from numpy.typing import ArrayLike, NDArray
from rss_lock.locked_dataset import LockedDataset

from access_io.access_attr_define import (  # atm_pars_era5_attributes_access,
    anc_var_attributes_access,
    common_global_attributes_access,
    coord_attributes_access,
    resamp_tb_attributes_access,
)

from polar_grids import NSIDC_ease2_grids

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
NUM_LATS = 721
NUM_LONS = 1440
NUM_HOURS = 24
NUM_X_POLE = 720
NUM_Y_POLE = 720


class OkToSkipDay(Exception):
    pass


def set_or_create_attr(var: Variable, attr_name: str, attr_value: Any) -> None:
    """seems like something like this should be part
    of the interface but I can not find it"""

    if attr_name in var.ncattrs():
        if attr_name != "_FillValue":
            var.setncattr(attr_name, attr_value)
        return
    var.UnusedNameAttribute = attr_value
    var.renameAttribute("UnusedNameAttribute", attr_name)
    return


def set_all_attrs(var: Variable, attrs: dict[str, Any]) -> None:
    for name, value in attrs.items():
        if (name != "_FillValue") and (value is not None):
            var.setncattr(name, value)


def get_access_output_filename_daily_folder(
    date: datetime.date,
    satellite: str,
    target_size: int,
    dataroot: Path,
    var: str,
    grid_type: str = "equirectangular",
    pole: str = "",
    ksat="13",
    look=0,
) -> Path:
    """
    Generates the output filename path for a daily satellite dataset.

    Args:
        date (datetime.date): The date of the dataset.
        satellite (str): The name of the satellite.
        target_size (int): The target size in kilometers.
        dataroot (Path): The root directory for the data.
        var (str): The variable name.
        grid_type (str, optional): The grid type. Defaults to "equirectangular".
        pole (str, optional): The pole type (valid only for grid_type="ease2").
                            Defaults to "".

    Returns:
        Path: The output filename path.

    Raises:
        ValueError: If the pole is not "north" or "south"
                    (valid only for grid_type="ease2").
        ValueError: If the target size is not specified
                    (valid only for grid_type="ease2").
        ValueError: If the grid type is not valid.
    """
    if satellite.lower() == "ssmi":
        dataroot = dataroot / f'f{ksat}'

    suffix = ".nc"
    if satellite.lower() == "smap":
        if look == 0:
            suffix = ".fore.nc"
        elif look == 1:
            suffix = ".aft.nc"
        else:
            raise ValueError("look must be 0 or 1")


    if grid_type == "equirectangular":
        if target_size > 0:
            return (
                dataroot
                / f"Y{date:%Y}"
                / f"M{date:%m}"
                / f"D{date:%d}"
                / f"{satellite.lower()}_{var}_{date:%Y_%m_%d}.{target_size:03d}km{suffix}"
            )
        else:
            return (
                dataroot
                / f"Y{date:%Y}"
                / f"M{date:%m}"
                / f"D{date:%d}"
                / f"{satellite.lower()}_{var}_{date:%Y_%m_%d}{suffix}"
            )
    elif grid_type == "ease2":
        if target_size > 0:
            pole = str.lower(pole)
            if pole in ["north", "south"]:
                return (
                    dataroot
                    / f"Y{date:%Y}"
                    / f"M{date:%m}"
                    / f"D{date:%d}"
                    / (
                        f"{satellite.lower()}_{var}_{date:%Y_%m_%d}."
                        f"{target_size:03d}km.{pole}{suffix}"
                    )
                )
            else:
                raise ValueError(f"Pole={pole} must be north or south")
        else:
            raise ValueError("Must specify target size")
    if target_size > 0:
        return (
            dataroot
            / f"Y{date:%Y}"
            / f"M{date:%m}"
            / f"D{date:%d}"
            / f"{satellite.lower()}_{var}_{date:%Y_%m_%d}.{target_size:03d}km{suffix}"
        )
    else:
        raise ValueError(f"Grid type {grid_type} not valid")


"""
def replace_var_in_daily_tb_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    var: ArrayLike,
    var_name: str,
    standard_name: Optional[str] = None,
    long_name: Optional[str] = None,
    valid_min: Optional[float] = None,
    valid_max: Optional[float] = None,
    units: Optional[str] = None,
    source: Optional[str] = None,
    cell_method: Optional[str] = None,
    dataroot: Path = ACCESS_ROOT,
) -> None:

    raise ValueError("Not yet updated")

    filename = get_access_output_filename_daily_folder(date, satellite, dataroot)
    with LockedDataset(filename, "a", 60) as root_grp:
        # read in existing variable
        try:
            v = root_grp.variables[var_name]
        except KeyError:
            raise KeyError(f"Variable: {var_name} in not dataset, can not replace")

        # replace any attributes that are passed
        if standard_name is not None:
            v.standard_name = standard_name
        else:
            v.standard_name = var_name

        if long_name is not None:
            v.long_name = long_name

        if valid_min is not None:
            v.valid_min = valid_min
        if valid_max is not None:
            v.valid_max = valid_max
        if units is not None:
            v.units = units
        if source is not None:
            v.source = source
        if cell_method is not None:
            v.cell_method = cell_method
        v.coordinates = "latitude longitude hours"

        v[:, :, :] = var
"""


def write_daily_lf_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    ksat: str,
    target_size: int,
    pole: str = None,
    version: str,
    lf_version: str,
    land_fraction: ArrayLike,
    dataroot: Path = ACCESS_ROOT,
    overwrite: bool,
    script_name: str,
    commit: str,
) -> None:
    """
    Writes a daily land fraction NetCDF file.

    Args:
        date (datetime.date): The date of the dataset.
        satellite (str): The name of the satellite.
        target_size (int): The target size in kilometers.
        pole (str, optional): The pole type. Defaults to None.
        version (str): The version string.
        lf_version (str): The land fraction version string.
        land_fraction (ArrayLike): The land fraction data.
        dataroot (Path, optional): The root directory for the data.
                                    Defaults to ACCESS_ROOT.
        overwrite (bool): Whether to overwrite the existing file if it exists.
        script_name (str): The name of the script that generated the file.
        commit (str): The commit string for the code that generated the file.

    Returns:
        None

    Raises:
        ValueError: If the pole is not "north" or "south".
    """

    tb_fill = -999.0
    lf_string = f"land_frac_{lf_version}"
    # filename = get_access_output_filename_daily_folder(
    #     date, satellite, target_size, dataroot, lf_string, ksat=ksat
    # )

    if pole is None:
        filename = get_access_output_filename_daily_folder(
            date, satellite, target_size, dataroot, lf_string, ksat=ksat
        )
        lats = np.arange(0, NUM_LATS) * 0.25 - 90.0
        lons = np.arange(0, NUM_LONS) * 0.25
    elif pole in ["north", "south"]:
        filename = get_access_output_filename_daily_folder(
            date,
            satellite,
            target_size,
            dataroot,
            lf_string,
            grid_type="ease2",
            pole=pole,
            ksat=ksat
        )
        if pole == "north":
            lats = ease2_grid_25km_north.latitude
            lons = ease2_grid_25km_north.longitude
        else:  # must be south
            lats = ease2_grid_25km_south.latitude
            lons = ease2_grid_25km_south.longitude
    else:
        raise ValueError(f"pole = {pole} is not valid")

    if filename.is_file() and not overwrite:
        print(f"daily file for {date} exists... skipping")
        return
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

    lat_attrs = coord_attributes_access("latitude", dtype=np.float32)
    lon_attrs = coord_attributes_access("longitude", dtype=np.float32)

    # with netcdf_dataset(filename, "w", format="NETCDF4") as nc_out:
    os.makedirs(filename.parent, exist_ok=True)
    with LockedDataset(filename, "w", 60) as nc_out:
        set_all_attrs(nc_out, global_attrs)

        nc_out.createDimension("latitude", NUM_LATS)
        nc_out.createDimension("longitude", NUM_LONS)

        latitude = nc_out.createVariable(
            "latitude", "f4", ("latitude",), fill_value=-999.0
        )
        longitude = nc_out.createVariable(
            "longitude", "f4", ("longitude",), fill_value=-999.0
        )

        lf = nc_out.createVariable(
            "land_fraction",
            "f4",
            ("latitude", "longitude"),
            zlib=True,
            fill_value=var_attrs["_FillValue"],
            least_significant_digit=3,
        )

        set_all_attrs(latitude, lat_attrs)
        set_all_attrs(longitude, lon_attrs)
        set_all_attrs(lf, var_attrs)

        latitude[:] = lats
        longitude[:] = lons

        lf_to_put = np.nan_to_num(
            land_fraction, nan=tb_fill, posinf=tb_fill, neginf=tb_fill
        ).astype(np.float32)

        lf[:, :] = lf_to_put


def write_daily_tb_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    ksat: str = '13',
    target_size: int,
    pole: str = None,
    version: str,
    tb_array_by_hour: NDArray[Any],
    time_array_by_hour: NDArray[Any],
    dataroot: Path = ACCESS_ROOT,
    freq_list: NDArray[Any],
    file_list: Optional[Sequence[Path]],
    look: int = 0,
    script_name: str = "unavailable",
    commit: str = "unavailable",
) -> None:
    """
    Writes a daily brightness temperature NetCDF file.

    Args:
        date (datetime.date): The date of the dataset.
        satellite (str): The name of the satellite.
        target_size (int): The target size in kilometers.
        pole (str, optional): The pole type. Defaults to None.
        version (str): The version string.
        tb_array_by_hour (NDArray[Any]): The brightness temperature array by hour.
        time_array_by_hour (NDArray[Any]): The time array by hour.
        dataroot (Path, optional): The root directory for the data.
            Defaults to ACCESS_ROOT.
        freq_list (NDArray[Any]): The frequency list.
        file_list (Optional[Sequence[Path]]): The list of file paths. Defaults to None.
        script_name (str, optional): The name of the script. Defaults to "unavailable".
        commit (str, optional): The commit string. Defaults to "unavailable".

    Returns:
        None

    Raises:
        ValueError: If the pole is not "north" or "south".
    """

    if pole is None:
        filename = get_access_output_filename_daily_folder(
            date, satellite, target_size, dataroot, "resamp_tbs",ksat=ksat,look=look
        )
        lats = np.arange(0, NUM_LATS) * 0.25 - 90.0
        lons = np.arange(0, NUM_LONS) * 0.25
    elif pole in ["north", "south"]:
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
        else:  # must be south
            lats = ease2_grid_25km_south.latitude
            lons = ease2_grid_25km_south.longitude
    else:
        raise ValueError(f"pole = {pole} is not valid")

    os.makedirs(filename.parent, exist_ok=True)

    num_freq = len(freq_list)
    if satellite.lower() == "smap":
        num_pol = 4
    else:
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
        set_all_attrs(nc_out, glb_attrs)

        # create Dimensions
        nc_out.createDimension("latitude", NUM_LATS)
        nc_out.createDimension("longitude", NUM_LONS)
        nc_out.createDimension("hours", NUM_HOURS)
        nc_out.createDimension("freq", num_freq)
        nc_out.createDimension("pol", num_pol)

        latitude = nc_out.createVariable(
            "latitude", "f4", ("latitude",), fill_value=-999.0
        )
        longitude = nc_out.createVariable(
            "longitude", "f4", ("longitude",), fill_value=-999.0
        )
        hours = nc_out.createVariable("hours", "i4", ("hours",))
        freq = nc_out.createVariable("freq", "f4", ("freq",))
        pol = nc_out.createVariable("pol", "i4", ("pol",))

        time = nc_out.createVariable(
            "time",
            "i8",
            (
                "latitude",
                "longitude",
                "hours",
            ),
            zlib=True,
            fill_value=-999999,
        )
        tbs = nc_out.createVariable(
            "brightness_temperature",
            "f4",
            (
                "latitude",
                "longitude",
                "hours",
                "freq",
                "pol",
            ),
            zlib=True,
            fill_value=tb_attrs["_FillValue"],
            least_significant_digit=2,
        )

        # set coordinate attributes from .json files

        lat_attrs = coord_attributes_access("latitude", dtype=np.float32)
        lon_attrs = coord_attributes_access("longitude", dtype=np.float32)
        hours_attrs = coord_attributes_access("hours", date=date, dtype=np.int32)
        freq_attrs = coord_attributes_access("frequency", dtype=np.float32)
        pol_attrs = coord_attributes_access("polarization", dtype=np.int32)
        time_attrs = coord_attributes_access("time", dtype=np.int64)

        set_all_attrs(latitude, lat_attrs)
        set_all_attrs(longitude, lon_attrs)
        set_all_attrs(hours, hours_attrs)
        set_all_attrs(freq, freq_attrs)
        set_all_attrs(pol, pol_attrs)
        set_all_attrs(time, time_attrs)

        set_all_attrs(tbs, tb_attrs)

        # write the coordinate data
        latitude[:] = lats
        longitude[:] = lons
        hours[:] = np.arange(NUM_HOURS)
        freq[:] = freq_list
        
        pol[:] = np.arange(num_pol, dtype=np.int32)

        # write the time layer
        time_to_put = (
            time_array_by_hour + (date - datetime.date(1900, 1, 1)).total_seconds()
        )
        time_to_put = np.nan_to_num(
            time_to_put, nan=-999998.99, posinf=-999998.99, neginf=-999998.99
        )
        time_to_put = np.floor(time_to_put).astype(np.int64)
        time[:, :, :] = time_to_put

        fill_val = np.float32(tb_attrs["_FillValue"])
        tbs_to_put = np.nan_to_num(  # type: ignore
            tb_array_by_hour,
            nan=fill_val,
            posinf=fill_val,
            neginf=fill_val,
        ).astype(np.float32)
        tbs[:, :, :, :, :] = tbs_to_put


def edit_attrs_daily_tb_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    target_size: int,
    version: str,
    tb_array_by_hour: ArrayLike,
    time_array_by_hour: ArrayLike,
    dataroot: Path = ACCESS_ROOT,
    freq_list: ArrayLike,
    file_list: Optional[Sequence[Path]],
    script_name: str = "unavailable",
    commit: str = "unavailable",
) -> None:
    filename = get_access_output_filename_daily_folder(
        date, satellite, target_size, dataroot, "resamp_tbs"
    )

    # with netcdf_dataset(filename, "w", format="NETCDF4") as nc_out:
    with LockedDataset(filename, "r+", 60) as nc_out:
        # set the global_attributes
        try:
            data_type = tb_array_by_hour.dtype
        except AttributeError:
            data_type = np.float32

        glb_attrs = common_global_attributes_access(
            date, satellite, target_size, version, data_type
        )
        tb_attrs = resamp_tb_attributes_access(satellite, version, data_type)

        glb_attrs.update(tb_attrs["global"])
        tb_attrs = tb_attrs["var"]

        glb_attrs["spatial_resolution"] = f"{target_size} km X {target_size} km"
        glb_attrs[
            "cell_method"
        ] = f"area: resampled to {target_size}km Gaussian footprint"
        set_all_attrs(nc_out, glb_attrs)
        set_all_attrs(nc_out["brightness_temperature"], tb_attrs)

        lat_attrs = coord_attributes_access("latitude", dtype=np.float32)
        lon_attrs = coord_attributes_access("longitude", dtype=np.float32)
        hours_attrs = coord_attributes_access("hours", date=date, dtype=np.int32)
        freq_attrs = coord_attributes_access("frequency", dtype=np.float32)
        pol_attrs = coord_attributes_access("polarization", dtype=np.int32)
        time_attrs = coord_attributes_access("time", dtype=np.int64)

        set_all_attrs(nc_out["latitude"], lat_attrs)
        set_all_attrs(nc_out["longitude"], lon_attrs)
        set_all_attrs(nc_out["hours"], hours_attrs)
        set_all_attrs(nc_out["freq"], freq_attrs)
        set_all_attrs(nc_out["pol"], pol_attrs)
        set_all_attrs(nc_out["time"], time_attrs)


def write_daily_ancillary_var_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    target_size: int,
    anc_data: NDArray[Any],
    anc_name: str,
    anc_attrs: dict[str, Any],
    global_attrs: Union[dict[str, Any], Literal["copy"]],
    dataroot: Path = ACCESS_ROOT,
) -> None:
    """
    Writes a daily NetCDF file for an ancillary variable.

    Args:
        date (datetime.date): The date of the dataset.
        satellite (str): The name of the satellite.
        target_size (int): The target size in kilometers.
        anc_data (NDArray[Any]): The ancillary data array.
        anc_name (str): The name of the ancillary variable.
        anc_attrs (dict[str, Any]): The attributes for the ancillary variable.
        global_attrs (Union[dict[str, Any], Literal["copy"]]):
            the global attributes for the NetCDF file.
            It can be a dictionary or the string "copy" to copy attributes
                from another file.
        dataroot (Path, optional): The root directory for the data.
            Defaults to ACCESS_ROOT.

    Returns:
        None
    """

    base_filename = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, dataroot, "resamp_tbs"
    )
    var_filename = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, dataroot, f"{anc_name}_temp"
    )
    var_filename_final = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, dataroot, anc_name
    )
    with suppress(FileNotFoundError):
        var_filename.unlink()

    with LockedDataset(var_filename, "w", 60) as nc_out:
        with LockedDataset(base_filename, "r", 60) as root_grp:
            # Create the dimensions of the file
            for name, dim in root_grp.dimensions.items():
                if name in ["hours", "latitude", "longitude"]:
                    nc_out.createDimension(
                        name, len(dim) if not dim.isunlimited() else None
                    )

            set_all_attrs(nc_out, global_attrs)

            for var_name in ["hours", "longitude", "latitude"]:
                # Create the time and dimension variables in the output file
                var_in = root_grp[var_name]
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
    var_filename_final = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, dataroot, anc_name
    )

    with LockedDataset(var_filename_final, "r+", 60) as nc_out:
        set_all_attrs(nc_out, global_attrs)
        set_all_attrs(nc_out.variables[anc_name], anc_attrs)

        lat_attrs = coord_attributes_access("latitude", dtype=np.float32)
        lon_attrs = coord_attributes_access("longitude", dtype=np.float32)
        hours_attrs = coord_attributes_access("hours", date=date, dtype=np.int32)

        set_all_attrs(nc_out.variables["latitude"], lat_attrs)
        set_all_attrs(nc_out.variables["longitude"], lon_attrs)
        set_all_attrs(nc_out.variables["hours"], hours_attrs)


def write_ocean_emiss_to_daily_ACCESS(
    *,
    ocean_emiss: ArrayLike,
    current_day: datetime.date,
    satellite: str,
    ksat: str = "13",
    target_size: int,
    glb_attrs: dict[str, Any],
    var_attrs: dict[str, Any],
    dataroot: Path,
    outputroot: Path,
    verbose: bool = False,
) -> None:
    """
    Writes the ocean emissivity data to a daily ACCESS NetCDF file.

    Args:
        ocean_emiss (ArrayLike): The ocean emissivity data array.
        current_day (datetime.date): The current day.
        satellite (str): The name of the satellite.
        target_size (int): The target size in kilometers.
        glb_attrs (dict[str, Any]): The global attributes for the NetCDF file.
        var_attrs (dict[str, Any]): The variable attributes for the emissivity variable.
        dataroot (Path): The root directory for the input data.
        outputroot (Path): The root directory for the output data.
        verbose (bool, optional): Whether to print verbose output. Defaults to False.

    Returns:
        None

    Raises:
        OkToSkipDay: If the base file is not found, indicating it is OK to skip the day.
    """

    if verbose:
        print(f"Opening base file for {satellite} on {current_day} in {dataroot}")

    if satellite.lower() == "amsr2":
        from satellite_definitions.amsr2 import REF_FREQ

    base_filename = get_access_output_filename_daily_folder(
        current_day, satellite, target_size, dataroot, "resamp_tbs",ksat=ksat
    )
    emiss_filename_final = get_access_output_filename_daily_folder(
        current_day, satellite, target_size, outputroot, "ocean_emiss_era5",ksat=ksat
    )

    try:
        with LockedDataset(base_filename, "r") as root_grp:
            os.makedirs(emiss_filename_final.parent, exist_ok=True)

            with LockedDataset(emiss_filename_final, mode="w") as trg:
                for name, dim in root_grp.dimensions.items():
                    if name in ["hours", "latitude", "longitude","freq","pol"]:
                        trg.createDimension(
                            name, len(dim) if not dim.isunlimited() else None
                        )

                # create the new dimensions needed for the emissivity data
                #trg.createDimension("freq", len(REF_FREQ))
                #trg.createDimension("pol", 2)

                set_all_attrs(trg, glb_attrs)

                for var_name in ["hours", "longitude", "latitude", "freq", "pol"]:
                    # Create the time and dimension variables in the output file
                    var_in = root_grp[var_name]
                    trg.createVariable(
                        var_name, var_in.dtype, var_in.dimensions, zlib=True
                    )

                    # Copy the attributes from the base file
                    trg.variables[var_name].setncatts(
                        {a: var_in.getncattr(a) for a in var_in.ncattrs()}
                    )
                    trg[var_name][:] = var_in[:]

                dimensions_out = ("latitude", "longitude", "hours", "freq", "pol")

                for varname, long_name, units in [
                    ("emissivity", "ocean surface emissivity", None),
                ]:
                    print(f"starting writing {varname} to {emiss_filename_final}")
                    least_significant_digit = 3
                    trg.createVariable(
                        varname,
                        np.float32,
                        dimensions_out,
                        fill_value=var_attrs["_FillValue"],
                        zlib=True,
                        least_significant_digit=least_significant_digit,
                    )
                    set_all_attrs(trg, var_attrs)
                    trg[varname][:, :, :, :, :] = ocean_emiss

                    print(f"finished writing {varname}")

    except FileNotFoundError:
        raise OkToSkipDay
