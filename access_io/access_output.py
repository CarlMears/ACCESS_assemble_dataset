import datetime
import os
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
from numpy.typing import ArrayLike

from rss_lock.locked_dataset import LockedDataset

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
    "asmr2": [
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


def set_or_create_attr(var, attr_name, attr_value):
    # seems like something like this should be part of the interface but I can not find it
    if attr_name in var.ncattrs():
        var.setncattr(attr_name, attr_value)
        return
    var.UnusedNameAttribute = attr_value
    var.renameAttribute("UnusedNameAttribute", attr_name)
    return


def get_access_output_filename(
    date: datetime.date, satellite: str, dataroot: Path, var: str
) -> Path:
    return (
        dataroot
        / f"Y{date:%Y}"
        / f"M{date:%m}"
        / f"{satellite}_{var}_{date:%Y_%m_%d}.nc"
    )


def get_access_output_filename_daily_folder(
    date: datetime.date, satellite: str, dataroot: Path, var: str
) -> Path:
    return (
        dataroot
        / f"Y{date:%Y}"
        / f"M{date:%m}"
        / f"D{date:%d}"
        / f"{satellite}_{var}_{date:%Y_%m_%d}.nc"
    )


def append_var_to_daily_tb_netcdf(
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
    v_fill: float = -999.0,
    dataroot: Path = ACCESS_ROOT,
    overwrite: bool = False,
) -> None:

    filename = get_access_output_filename(date, satellite, dataroot)
    with LockedDataset(filename, "a", 60) as root_grp:
        # with netcdf_dataset(filename, "a", format="NETCDF4") as root_grp:
        try:
            v = root_grp.createVariable(
                var_name,
                "f4",
                (
                    "latitude",
                    "longitude",
                    "hours",
                ),
                zlib=True,
                fill_value=v_fill,
            )
        except RuntimeError:
            if overwrite:
                print(f"Warning: variable {var_name} already exists - overwriting")
                v = root_grp.variables[var_name]
            else:
                raise RuntimeError(f"Variable {var_name} already exists")

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


def write_var_corresponding_to_daily_tb_netcdf(
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
    v_fill: float = -999.0,
    dataroot: Path = ACCESS_ROOT,
    overwrite: bool = False,
) -> None:

    filename = get_access_output_filename(date, satellite, dataroot)
    with LockedDataset(filename, "a", 60) as root_grp:
        # with netcdf_dataset(filename, "a", format="NETCDF4") as root_grp:
        try:
            v = root_grp.createVariable(
                var_name,
                "f4",
                (
                    "latitude",
                    "longitude",
                    "hours",
                ),
                zlib=True,
                fill_value=v_fill,
            )
        except RuntimeError:
            if overwrite:
                print(f"Warning: variable {var_name} already exists - overwriting")
                v = root_grp.variables[var_name]
            else:
                raise RuntimeError(f"Variable {var_name} already exists")

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

    filename = get_access_output_filename(date, satellite, dataroot)
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


def append_const_var_to_daily_tb_netcdf(
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
    v_fill: float = -999.0,
    dataroot: Path = ACCESS_ROOT,
    overwrite: bool = False,
    verbose: bool = False,
    script_name: Optional[str] = None,
    commit: Optional[str] = None,
    lock_stale_time: float = 86400.0,
) -> None:
    filename = get_access_output_filename(date, satellite, dataroot)
    with LockedDataset(
        filename, "a", 60, lock_stale_time=lock_stale_time, verbose=verbose
    ) as root_grp:
        # with netcdf_dataset(filename, "a", format="NETCDF4") as root_grp:
        try:
            v = root_grp.createVariable(
                var_name,
                "f4",
                (
                    "latitude",
                    "longitude",
                ),
                zlib=True,
                fill_value=v_fill,
            )
        except RuntimeError:
            if overwrite:
                print(f"Warning: variable {var_name} already exists - overwriting")
                v = root_grp.variables[var_name]
            else:
                raise RuntimeError(f"Variable {var_name} already exists")

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
        if script_name is not None:
            v.script_name = script_name
        if commit is not None:
            v.commit = commit

        v.coordinates = "latitude longitude"

        v[:, :] = var


def append_lf_daily_tb_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    land_fraction: ArrayLike,
    dataroot: Path = ACCESS_ROOT,
) -> None:
    lf_fill = -999.0
    filename = get_access_output_filename(date, satellite, dataroot, "land_frac")
    os.makedirs(filename.parent, exist_ok=True)
    with LockedDataset(filename, "a", 60) as root_grp:
        # with netcdf_dataset(filename, "a", format="NETCDF4") as root_grp:
        lf = root_grp.createVariable(
            "land_fraction",
            "f4",
            (
                "latitude",
                "longitude",
            ),
            zlib=True,
            fill_value=lf_fill,
        )

        lf.standard_name = "land_fraction"
        lf.long_name = "land_fraction"
        lf.missing = lf_fill
        lf.valid_min = 0.0
        lf.valid_max = 1.0
        lf.coordinates = "latitude longitude"

        lf[:, :] = land_fraction


def write_daily_lf_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    land_fraction: ArrayLike,
    dataroot: Path = ACCESS_ROOT,
) -> None:

    tb_fill = -999.0
    filename = get_access_output_filename(date, satellite, dataroot, "land_frac")

    os.makedirs(filename.parent, exist_ok=True)

    lats = np.arange(0, NUM_LATS) * 0.25 - 90.0
    lons = np.arange(0, NUM_LONS) * 0.25

    # with netcdf_dataset(filename, "w", format="NETCDF4") as root_grp:
    with LockedDataset(filename, "w", 60) as root_grp:
        root_grp.Conventions = "CF-1.8"
        root_grp.standard_name_vocabulary = (
            "CF Standard Name Table (v78, 21 September 2021)"
        )
        root_grp.id = filename.name
        root_grp.title = f"Land Fraction Corresponding to resampled {satellite} brightness temperatures"
        root_grp.product_version = "v00r00"
        root_grp.date_issued = "2021-10-01"
        root_grp.summary = (
            "Remote Sensing Systems (RSS) Land Fraction for Resampled brightness temperature; "
            "intercalibrated and homogenized brightness temperature "
            "polar-orbiting resampled to a regular Earth grid"
        )

        root_grp.keywords_vocabulary = (
            "NASA Global Change Master Directory (GCMD) "
            "Earth Science Keywords, Version 6.0"
        )

        root_grp.cdm_data_type = "Grid"
        root_grp.program = (
            "NASA ACCESS-0031 > Machine Learning Datasets "
            "for the Earth’s Natural Microwave Emission"
        )

        root_grp.date_created = datetime.datetime.now().isoformat()
        root_grp.creator_name = "Carl Mears"
        root_grp.creator_url = "http://www.remss.com/"
        root_grp.creator_email = "mears@remss.com"
        root_grp.institution = "Remote Sensing Systems"
        root_grp.processing_level = "NASA Level 4"
        root_grp.references = "None"
        root_grp.history = (
            datetime.datetime.now().isoformat()
            + " Created Resampled Brightness Temperature from RSS AMSR2 L1A data"
        )
        root_grp.geospatial_lat_min = -90.0
        root_grp.geospatial_lat_max = 90.0
        root_grp.geospatial_lon_min = 0.0
        root_grp.geospatial_lon_max = 359.9999
        root_grp.geospatial_lat_units = "degrees_north"
        root_grp.geospatial_lon_units = "degrees_east"
        root_grp.spatial_resolution = "30 km X 30 km"

        root_grp.license = "No restrictions on access or use"
        root_grp.contributor_name = "Frank Wentz, Carl Mears"
        root_grp.contributor_role = (
            "Principal investigator and originator of "
            "input/source or antenna temperature data, "
            "Processor and author of entire driver routine "
            "which resamples RSS native brightness temperature to a fixed Earth grid"
        )

        root_grp.createDimension("latitude", NUM_LATS)
        root_grp.createDimension("longitude", NUM_LONS)
        root_grp.createDimension("hours", NUM_HOURS)
        root_grp.createDimension("channels", NUM_CHANNELS)

        latitude = root_grp.createVariable(
            "latitude", "f4", ("latitude",), fill_value=-999.0
        )
        longitude = root_grp.createVariable(
            "longitude", "f4", ("longitude",), fill_value=-999.0
        )
        hours = root_grp.createVariable("hours", "i4", ("hours",))
        channels = root_grp.createVariable("channels", "i4", ("channels",))

        time = root_grp.createVariable(
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
        tbs = root_grp.createVariable(
            "brightness_temperature",
            "f4",
            (
                "latitude",
                "longitude",
                "hours",
                "channels",
            ),
            zlib=True,
            fill_value=tb_fill,
            least_significant_digit=2,
        )

        latitude.standard_name = "latitude"
        latitude.long_name = "latitude"
        latitude.units = "degrees_north"
        latitude.valid_range = (-90.0, 90.0)

        longitude.standard_name = "longitude"
        longitude.long_name = "longitude"
        longitude.units = "degrees_east"
        longitude.valid_range = (0.0, 360.0)

        hours.standard_name = "time"
        hours.units = f"hours since {date.isoformat()} 00:00:00.0"
        hours.valid_min = 0
        hours.valid_max = 23

        channels.standard_name = "channel"
        channels.long_name = f"{satellite} channel index"
        channel_names = ", ".join(USED_CHANNELS["amsr2"])
        channels.channel_names = channel_names

        # Each day's file is offset by a half-hour into the previous day
        time.standard_name = "time"
        time.long_name = "time of satellite observation"
        time.units = "seconds since 1900-01-01 00:00:00.0"
        time.missing = -999999
        time.valid_range = 0
        time.valid_max = 200 * 366 * 24 * 3600
        time.coordinates = "latitude longitude"

        tbs.standard_name = "brightness_temperature"
        tbs.units = "degrees kelvin"
        tbs.missing = tb_fill
        tbs.valid_min = 50.0
        tbs.valid_max = 350.0
        tbs.long_name = f"resampled {satellite} brightness temperature on lat/lon grid"
        tbs.channel_names = channel_names
        tbs.coordinates = "latitude longitude"

        latitude[:] = lats
        longitude[:] = lons
        hours[:] = np.arange(NUM_HOURS)
        channels[:] = np.arange(NUM_CHANNELS) + 1

        time_to_put = (
            time_array_by_hour + (date - datetime.date(1900, 1, 1)).total_seconds()
        )
        time_to_put = np.nan_to_num(
            time_to_put, nan=-999998.99, posinf=-999998.99, neginf=-999998.99
        )
        time_to_put = np.floor(time_to_put).astype(np.int64)
        time[:, :, :] = time_to_put

        tbs_to_put = np.nan_to_num(
            tb_array_by_hour, nan=tb_fill, posinf=tb_fill, neginf=tb_fill
        ).astype(np.float32)
        tbs[:, :, :, :] = tbs_to_put


def write_daily_tb_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    tb_array_by_hour: ArrayLike,
    time_array_by_hour: ArrayLike,
    dataroot: Path = ACCESS_ROOT,
    file_list: Optional[Sequence[Path]],
) -> None:
    tb_fill = -999.0
    filename = get_access_output_filename_daily_folder(
        date, satellite, dataroot, "resamp_tbs"
    )
    os.makedirs(filename.parent, exist_ok=True)

    lats = np.arange(0, NUM_LATS) * 0.25 - 90.0
    lons = np.arange(0, NUM_LONS) * 0.25

    # with netcdf_dataset(filename, "w", format="NETCDF4") as root_grp:
    with LockedDataset(filename, "w", 60) as root_grp:
        root_grp.Conventions = "CF-1.8"
        root_grp.standard_name_vocabulary = (
            "CF Standard Name Table (v78, 21 September 2021)"
        )
        root_grp.id = filename.name
        root_grp.title = f"Resampled {satellite} brightness temperatures"
        root_grp.product_version = "v00r00"
        root_grp.date_issued = "2021-10-01"
        root_grp.summary = (
            "Remote Sensing Systems (RSS) Resampled brightness temperature; "
            "intercalibrated and homogenized brightness temperature "
            "polar-orbiting resampled to a regular Earth grid"
        )
        root_grp.keywords = (
            "EARTH SCIENCE > SPECTRAL/ENGINEERING > MICROWAVE > BRIGHTNESS TEMPERATURE"
        )
        root_grp.keywords_vocabulary = (
            "NASA Global Change Master Directory (GCMD) "
            "Earth Science Keywords, Version 6.0"
        )
        root_grp.platform = "GCOM-W1, JAXA"
        root_grp.sensor = "AMSR2 > Advanced Microwave Scanning Radiometer 2"
        root_grp.cdm_data_type = "Grid"
        root_grp.program = (
            "NASA ACCESS-0031 > Machine Learning Datasets "
            "for the Earth’s Natural Microwave Emission"
        )
        if file_list is not None:
            source_string = ", ".join(str(p) for p in file_list)
            root_grp.source = source_string
        root_grp.date_created = datetime.datetime.now().isoformat()
        root_grp.creator_name = "Carl Mears"
        root_grp.creator_url = "http://www.remss.com/"
        root_grp.creator_email = "mears@remss.com"
        root_grp.institution = "Remote Sensing Systems"
        root_grp.processing_level = "NASA Level 4"
        root_grp.references = "None"
        root_grp.history = (
            datetime.datetime.now().isoformat()
            + " Created Resampled Brightness Temperature from RSS AMSR2 L1A data"
        )
        root_grp.geospatial_lat_min = -90.0
        root_grp.geospatial_lat_max = 90.0
        root_grp.geospatial_lon_min = 0.0
        root_grp.geospatial_lon_max = 359.9999
        root_grp.geospatial_lat_units = "degrees_north"
        root_grp.geospatial_lon_units = "degrees_east"
        root_grp.spatial_resolution = "30 km X 30 km"
        day_boundary = datetime.datetime.combine(date, datetime.time())
        start_date = day_boundary - datetime.timedelta(minutes=30)
        end_date = day_boundary + datetime.timedelta(minutes=1410.0)
        root_grp.time_coverage_start = start_date.isoformat()
        root_grp.time_coverage_end = end_date.isoformat()
        root_grp.time_coverage_duration = "P24H"
        root_grp.license = "No restrictions on access or use"
        root_grp.contributor_name = "Frank Wentz, Carl Mears"
        root_grp.contributor_role = (
            "Principal investigator and originator of "
            "input/source or antenna temperature data, "
            "Processor and author of entire driver routine "
            "which resamples RSS native brightness temperature to a fixed Earth grid"
        )

        root_grp.createDimension("latitude", NUM_LATS)
        root_grp.createDimension("longitude", NUM_LONS)
        root_grp.createDimension("hours", NUM_HOURS)
        root_grp.createDimension("channels", NUM_CHANNELS)

        latitude = root_grp.createVariable(
            "latitude", "f4", ("latitude",), fill_value=-999.0
        )
        longitude = root_grp.createVariable(
            "longitude", "f4", ("longitude",), fill_value=-999.0
        )
        hours = root_grp.createVariable("hours", "i4", ("hours",))
        channels = root_grp.createVariable("channels", "i4", ("channels",))

        time = root_grp.createVariable(
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
        tbs = root_grp.createVariable(
            "brightness_temperature",
            "f4",
            (
                "latitude",
                "longitude",
                "hours",
                "channels",
            ),
            zlib=True,
            fill_value=tb_fill,
            least_significant_digit=2,
        )

        latitude.standard_name = "latitude"
        latitude.long_name = "latitude"
        latitude.units = "degrees_north"
        latitude.valid_range = (-90.0, 90.0)

        longitude.standard_name = "longitude"
        longitude.long_name = "longitude"
        longitude.units = "degrees_east"
        longitude.valid_range = (0.0, 360.0)

        hours.standard_name = "time"
        hours.units = f"hours since {date.isoformat()} 00:00:00.0"
        hours.valid_min = 0
        hours.valid_max = 23

        channels.standard_name = "channel"
        channels.long_name = f"{satellite} channel index"
        channel_names = ", ".join(USED_CHANNELS["amsr2"])
        channels.channel_names = channel_names

        # Each day's file is offset by a half-hour into the previous day
        time.standard_name = "time"
        time.long_name = "time of satellite observation"
        time.units = "seconds since 1900-01-01 00:00:00.0"
        time.missing = -999999
        time.valid_range = 0
        time.valid_max = 200 * 366 * 24 * 3600
        time.coordinates = "latitude longitude"

        tbs.standard_name = "brightness_temperature"
        tbs.units = "degrees kelvin"
        tbs.missing = tb_fill
        tbs.valid_min = 50.0
        tbs.valid_max = 350.0
        tbs.long_name = f"resampled {satellite} brightness temperature on lat/lon grid"
        tbs.channel_names = channel_names
        tbs.coordinates = "latitude longitude"

        latitude[:] = lats
        longitude[:] = lons
        hours[:] = np.arange(NUM_HOURS)
        channels[:] = np.arange(NUM_CHANNELS) + 1

        time_to_put = (
            time_array_by_hour + (date - datetime.date(1900, 1, 1)).total_seconds()
        )
        time_to_put = np.nan_to_num(
            time_to_put, nan=-999998.99, posinf=-999998.99, neginf=-999998.99
        )
        time_to_put = np.floor(time_to_put).astype(np.int64)
        time[:, :, :] = time_to_put

        tbs_to_put = np.nan_to_num(
            tb_array_by_hour, nan=tb_fill, posinf=tb_fill, neginf=tb_fill
        ).astype(np.float32)
        tbs[:, :, :, :] = tbs_to_put


def write_daily_ancillary_var_netcdf(
    *,
    date: datetime.date,
    satellite: str,
    anc_data: ArrayLike,
    anc_name: str,
    anc_attrs: dict,
    global_attrs: Union[dict, str],
    dataroot: Path = ACCESS_ROOT,
) -> None:

    base_filename = get_access_output_filename_daily_folder(
        date, satellite.lower(), dataroot, "resamp_tbs"
    )
    var_filename = get_access_output_filename_daily_folder(
        date, satellite.lower(), dataroot, f"{anc_name}_temp"
    )
    var_filename_final = get_access_output_filename_daily_folder(
        date, satellite.lower(), dataroot, anc_name
    )

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

            # make the rain rate variable with the same dimensions as the time variable in the base file
            new_var = nc_out.createVariable(
                anc_name, np.float32, nc_out.variables["time"].dimensions, zlib=True
            )
            for key in anc_attrs.keys():
                val = anc_attrs[key]
                set_or_create_attr(new_var, key, val)

            new_var.coordinates = "latitude longitude"

            new_var[:, :, :] = anc_data[:, :, :]

    # everything written -- rename to final file name
    try:
        # delete the file is
        var_filename_final.unlink()
    except FileNotFoundError:
        pass
    var_filename.rename(var_filename_final)
