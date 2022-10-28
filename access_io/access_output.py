import datetime
import os
from pathlib import Path
from typing import Any, Optional, Sequence, Union

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


class OkToSkipDay(Exception):
    pass


def set_or_create_attr(var, attr_name, attr_value):
    """seems like something like this should be part
    of the interface but I can not find it"""

    if attr_name in var.ncattrs():
        if attr_name != "_FillValue":
            var.setncattr(attr_name, attr_value)
        return
    var.UnusedNameAttribute = attr_value
    var.renameAttribute("UnusedNameAttribute", attr_name)
    return


def get_access_output_filename_daily_folder(
    date: datetime.date, satellite: str, target_size: int, dataroot: Path, var: str
) -> Path:

    if target_size > 0:
        return (
            dataroot
            / f"Y{date:%Y}"
            / f"M{date:%m}"
            / f"D{date:%d}"
            / f"{satellite.lower()}_{var}_{date:%Y_%m_%d}.{target_size:03d}km.nc"
        )
    else:
        return (
            dataroot
            / f"Y{date:%Y}"
            / f"M{date:%m}"
            / f"D{date:%d}"
            / f"{satellite.lower()}_{var}_{date:%Y_%m_%d}.nc"
        )


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


def write_daily_lf_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    target_size: int,
    version: str,
    lf_version: str,
    land_fraction: ArrayLike,
    dataroot: Path = ACCESS_ROOT,
    overwrite: bool,
) -> None:

    tb_fill = -999.0
    lf_string = f"land_frac_{lf_version}"
    filename = get_access_output_filename_daily_folder(
        date, satellite, target_size, dataroot, lf_string
    )

    if filename.is_file() and not overwrite:
        print(f"daily file for {date} exists... skipping")
        return []
    else:
        try:
            filename.unlink()
        except FileNotFoundError:
            pass

    os.makedirs(filename.parent, exist_ok=True)

    lats = np.arange(0, NUM_LATS) * 0.25 - 90.0
    lons = np.arange(0, NUM_LONS) * 0.25

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
    ] = f"area: resampled to {target_size}km guassian footprint"

    # with netcdf_dataset(filename, "w", format="NETCDF4") as nc_out:
    with LockedDataset(filename, "w", 60) as nc_out:

        for key in global_attrs.keys():
            value = global_attrs[key]
            set_or_create_attr(nc_out, key, value)

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

        lat_attrs = coord_attributes_access("latitude")
        for key in lat_attrs.keys():
            value = lat_attrs[key]
            set_or_create_attr(latitude, key, value)

        lon_attrs = coord_attributes_access("longitude")
        for key in lon_attrs.keys():
            value = lon_attrs[key]
            set_or_create_attr(longitude, key, value)

        for key in var_attrs.keys():
            value = var_attrs[key]
            set_or_create_attr(lf, key, value)

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
    target_size: int,
    version: str,
    tb_array_by_hour: ArrayLike,
    time_array_by_hour: ArrayLike,
    dataroot: Path = ACCESS_ROOT,
    freq_list: ArrayLike,
    file_list: Optional[Sequence[Path]],
) -> None:

    filename = get_access_output_filename_daily_folder(
        date, satellite, target_size, dataroot, "resamp_tbs"
    )
    os.makedirs(filename.parent, exist_ok=True)

    lats = np.arange(0, NUM_LATS) * 0.25 - 90.0
    lons = np.arange(0, NUM_LONS) * 0.25

    num_freq = len(freq_list)
    num_pol = 2

    # with netcdf_dataset(filename, "w", format="NETCDF4") as nc_out:
    with LockedDataset(filename, "w", 60) as nc_out:

        # set the global_attributes

        glb_attrs = common_global_attributes_access(
            date, satellite, target_size, version
        )
        tb_attrs = resamp_tb_attributes_access(satellite, version)

        glb_attrs.update(tb_attrs["global"])
        tb_attrs = tb_attrs["var"]

        glb_attrs["spatial_resolution"] = f"{target_size} km X {target_size} km"
        glb_attrs[
            "cell_method"
        ] = f"area: resampled to {target_size}km guassian footprint"
        glb_attrs["date_accessed"] = f"{datetime.datetime.today()}"
        glb_attrs["creation_date"] = f"{datetime.datetime.today()}"

        for key in glb_attrs.keys():
            value = glb_attrs[key]
            set_or_create_attr(nc_out, key, value)

        nc_out.id = filename.name

        if file_list is not None:
            source_string = ", ".join(str(p) for p in file_list)
            nc_out.source = source_string

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

        lat_attrs = coord_attributes_access("latitude")
        for key in lat_attrs.keys():
            value = lat_attrs[key]
            set_or_create_attr(latitude, key, value)

        lon_attrs = coord_attributes_access("longitude")
        for key in lon_attrs.keys():
            value = lon_attrs[key]
            set_or_create_attr(longitude, key, value)

        hours_attrs = coord_attributes_access("hours", date=date)
        for key in hours_attrs.keys():
            value = hours_attrs[key]
            set_or_create_attr(hours, key, value)

        freq_attrs = coord_attributes_access("frequency")
        for key in freq_attrs.keys():
            value = freq_attrs[key]
            set_or_create_attr(freq, key, value)

        pol_attrs = coord_attributes_access("polarization")
        for key in pol_attrs.keys():
            value = pol_attrs[key]
            set_or_create_attr(pol, key, value)

        time_attrs = coord_attributes_access("time")
        for key in time_attrs.keys():
            value = time_attrs[key]
            set_or_create_attr(time, key, value)

        # set tb attributes from .json file
        for key in tb_attrs.keys():
            value = tb_attrs[key]
            set_or_create_attr(tbs, key, value)

        # write the coordinate data
        latitude[:] = lats
        longitude[:] = lons
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
        time[:, :, :] = time_to_put

        fill_val = np.float32(tb_attrs["_FillValue"])
        tbs_to_put = np.nan_to_num(
            tb_array_by_hour, nan=fill_val, posinf=fill_val, neginf=fill_val
        ).astype(np.float32)
        tbs[:, :, :, :, :] = tbs_to_put


def write_daily_ancillary_var_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    target_size: int,
    anc_data: ArrayLike,
    anc_name: str,
    anc_attrs: dict[str, Any],
    global_attrs: Union[dict[str, Any], str],
    dataroot: Path = ACCESS_ROOT,
) -> None:

    base_filename = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, dataroot, "resamp_tbs"
    )
    var_filename = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, dataroot, f"{anc_name}_temp"
    )
    var_filename_final = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, dataroot, anc_name
    )
    try:
        var_filename.unlink()
    except FileNotFoundError:
        pass

    with LockedDataset(var_filename, "w", 60) as nc_out:
        with LockedDataset(base_filename, "r", 60) as root_grp:
            # Create the dimensions of the file
            for name, dim in root_grp.dimensions.items():
                if name in ["hours", "latitude", "longitude"]:
                    nc_out.createDimension(
                        name, len(dim) if not dim.isunlimited() else None
                    )

            if global_attrs == "copy":
                # Copy the global attributes
                nc_out.setncatts({a: root_grp.getncattr(a) for a in root_grp.ncattrs()})
            else:
                for key in global_attrs.keys():
                    value = global_attrs[key]
                    set_or_create_attr(nc_out, key, value)

            for var_name in ["time", "hours", "longitude", "latitude"]:
                # Create the time and dimension variables in the output file
                var_in = root_grp[var_name]
                nc_out.createVariable(
                    var_name, var_in.dtype, var_in.dimensions, zlib=True
                )

                # Copy the attributes
                nc_out.variables[var_name].setncatts(
                    {a: var_in.getncattr(a) for a in var_in.ncattrs()}
                )

                # Copy the time values
                nc_out.variables[var_name][:] = root_grp.variables[var_name][:]

            # make the ancillary variable with the same dimensions as the
            # time variable in the base file
            # fill value is set from the var_attrs dictionary
            new_var = nc_out.createVariable(
                anc_name,
                np.float32,
                nc_out.variables["time"].dimensions,
                zlib=True,
                fill_value=np.float32(anc_attrs["_FillValue"]),
            )

            # convert numerical attributes to numbers
            attrs_to_convert_to_floats = ["valid_max", "valid_min", "v_fill"]
            for attr_to_convert in attrs_to_convert_to_floats:
                if attr_to_convert in anc_attrs.keys():
                    anc_attrs[attr_to_convert] = np.float32(anc_attrs[attr_to_convert])

            for key in anc_attrs.keys():
                if (
                    key != "_FillValue"
                ):  # _FillValue already set by the variable definition
                    new_var.setncattr(key, anc_attrs[key])

            new_var.coordinates = "latitude longitude"
            new_var[:, :, :] = anc_data[:, :, :]

    # everything written -- rename to final file name
    try:
        # delete the final file to make room for the new one
        var_filename_final.unlink()
    except FileNotFoundError:
        pass
    var_filename.rename(var_filename_final)


def write_ocean_emiss_to_daily_ACCESS(
    *,
    ocean_emiss: ArrayLike,
    current_day: datetime.date,
    satellite: str,
    target_size: int,
    glb_attrs: dict[str, Any],
    var_attrs: dict[str, Any],
    dataroot: Path,
    outputroot: Path,
    verbose: bool = False,
) -> None:

    if verbose:
        print(f"Opening base file for {satellite} on {current_day} in {dataroot}")

    if satellite.lower() == "amsr2":
        from satellite_definitions.amsr2 import REF_FREQ

    base_filename = get_access_output_filename_daily_folder(
        current_day, satellite, target_size, dataroot, "resamp_tbs"
    )
    try:
        with LockedDataset(base_filename, "r") as root_grp:

            emiss_filename_final = get_access_output_filename_daily_folder(
                current_day, satellite, target_size, outputroot, "ocean_emiss_era5"
            )
            os.makedirs(emiss_filename_final.parent, exist_ok=True)

            with LockedDataset(emiss_filename_final, mode="w") as trg:
                for name, dim in root_grp.dimensions.items():
                    if name in ["hours", "latitude", "longitude"]:
                        trg.createDimension(
                            name, len(dim) if not dim.isunlimited() else None
                        )

                # create the new dimensions needed for the emissivity data
                trg.createDimension("freq", len(REF_FREQ))
                trg.createDimension("pol", 2)

                for key in glb_attrs.keys():
                    value = glb_attrs[key]
                    set_or_create_attr(trg, key, value)

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
                    print(f"starting writing {varname}")
                    least_significant_digit = 3
                    trg.createVariable(
                        varname,
                        np.float32,
                        dimensions_out,
                        fill_value=var_attrs["_FillValue"],
                        zlib=True,
                        least_significant_digit=least_significant_digit,
                    )

                    for key in var_attrs.keys():
                        value = var_attrs[key]
                        set_or_create_attr(trg[varname], key, value)

                    trg[varname][:, :, :, :, :] = ocean_emiss

                    print(f"finished writing {varname}")
                    print()

    except FileNotFoundError:
        raise OkToSkipDay
