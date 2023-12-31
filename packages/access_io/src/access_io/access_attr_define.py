import datetime
import json
import os
from pathlib import Path
from typing import Optional, Any, cast
import numpy as np

"""
The routines in this file define the attributes for the ACCESS project output file
Constant attributes are read in from the .json files in attr_define_root
"""
if os.name == "nt":
    _attr_define_root = Path(
        "M:/job_access/python/dataset_assembly/packages/access_io/src/access_io/attr_define_json"
    )
elif os.name == "posix":
    _attr_define_root = Path(
        "/mnt/m/job_access/python/dataset_assembly/packages/access_io/src/access_io/attr_define_json"
    )

_allowed_numeric_types = ["np.float32", "np.float64", "np.int32", "np.int64"]

# the types of these attrs are set to match the type of the variable
matched_attributes = ["valid_min", "valid_max", "_FillValue", "missing_value"]

# the types of these attrs are set to np.float32
float32_attributes = [
    "geospatial_lat_min",
    "geospatial_lat_max",
    "geospatial_lon_min",
    "geospatial_lon_max",
]


def load_attrs(
    *, project: Optional[str] = None, satellite: Optional[str] = None, var: str
) -> dict[str, Any]:
    """
    Loads attributes from a JSON file based on the specified project,
    satellite, and variable.

    Parameters:
        project (Optional[str]): The project name. Default is None.
        satellite (Optional[str]): The satellite name. Default is None.
        var (str): The variable name.

    Returns:
        dict[str, Any]: The loaded attributes.
    """

    if len(var) == 0:
        raise ValueError("var must be specified")

    filename = "attr."
    if project is not None:
        filename = f"{filename}{project}."
    if satellite is not None:
        filename = f"{filename}{satellite.upper()}."

    filename = f"{filename}{var}.json"
    path_to_file = _attr_define_root / filename
    with open(path_to_file) as json_file:
        attrs = json.load(json_file)

    return cast(dict[str, Any], attrs)


def load_access_attrs(*, satellite: Optional[str] = None, var: str) -> dict:
    """
    Loads access attributes from a JSON file based on the specified
    satellite and variable.

    Parameters:
        satellite (Optional[str]): The satellite name. Default is None.
        var (str): The variable name.

    Returns:
        dict: The loaded access attributes.
    """

    project = "access"

    return load_attrs(project=project, satellite=satellite, var=var)


def fix_attr_types(attrs: dict, var_dtype):
    """
    Fixes the attribute types in the given dictionary based on
    the specified variable data type.

    Parameters:
        attrs (dict): The attributes dictionary.
        var_dtype: The variable data type.

    Returns:
        dict: The fixed attributes dictionary.
    """

    for att_name, value in attrs.items():
        if isinstance(value, dict):
            fix_attr_types(value, var_dtype)
        else:
            if att_name in matched_attributes:
                attrs[att_name] = np.asarray(value, dtype=var_dtype)
            elif att_name in float32_attributes:
                attrs[att_name] = np.asarray(value, dtype=np.float32)

    return attrs


def common_global_attributes_access(
    date: datetime.date,
    satellite: str,
    target_size: int,
    version: str = "v00r00",
    dtype=np.float32,
) -> dict:
    """
    Retrieves common global attributes for a specific date,
    satellite, target size, and version.

    Parameters:
        date (datetime.date): The date.
        satellite (str): The satellite name.
        target_size (int): The target size.
        version (str): The version. Default is "v00r00".
        dtype: The data type. Default is np.float32.

    Returns:
        dict: The common global attributes.
    """

    attrs = load_access_attrs(var="common")

    day_boundary = datetime.datetime.combine(date, datetime.time())
    start_date = day_boundary - datetime.timedelta(minutes=30)
    end_date = day_boundary + datetime.timedelta(minutes=1410.0)

    attrs["title"] = f"Resampled {satellite} brightness temperatures"
    attrs["date_created"] = datetime.datetime.now().isoformat()
    attrs["spatial_resolution"] = f"{target_size} km X {target_size} km"
    attrs["time_coverage_start"] = start_date.isoformat()
    attrs["time_coverage_end"] = end_date.isoformat()
    if version is not None:
        attrs["version"] = version
    else:
        try:
            del attrs["version"]
        except KeyError:
            pass

    attrs = fix_attr_types(attrs, dtype)
    return attrs


def resamp_tb_attributes_access(satellite: str, version="v01r00", dtype=np.float32):
    """
    Retrieves resampled brightness temperature attributes for a
    specific satellite and version.

    Parameters:
        satellite (str): The satellite name.
        version (str): The version. Default is "v01r00".
        dtype: The data type. Default is np.float32.

    Returns:
        dict: The resampled brightness temperature attributes.
    """

    attrs = load_access_attrs(satellite=satellite, var="resamp_tbs")
    if version is not None:
        attrs["global"]["version"] = version

    attrs = fix_attr_types(attrs, dtype)
    return attrs


def atm_pars_era5_attributes_access(
    satellite: str, target_size: int, version="v00r00", dtype=np.float32
):
    """
    Retrieves atmospheric parameters (ERA5) attributes for a
    specific satellite, target size, and version.

    Parameters:
        satellite (str): The satellite name.
        target_size (int): The target size.
        version (str): The version. Default is "v00r00".
        dtype: The data type. Default is np.float32.

    Returns:
        dict: The atmospheric parameters (ERA5) attributes.
    """

    attrs = load_access_attrs(satellite=satellite, var="atm_pars_era5")
    attrs["global"]["version"] = version
    attrs["global"]["date_accessed"] = f"{datetime.datetime.now()}"
    attrs = fix_attr_types(attrs, dtype)
    return attrs


def anc_var_attributes_access(
    satellite: str, var: str, version="v00r00", dtype=np.float32
):
    """
    Retrieves ancillary variable attributes for a specific
    satellite, variable, and version.

    Parameters:
        satellite (str): The satellite name.
        var (str): The variable name.
        version (str): The version. Default is "v00r00".
        dtype: The data type. Default is np.float32.

    Returns:
        dict: The ancillary variable attributes.
    """

    attrs = load_access_attrs(satellite=satellite, var=var)

    if "global" in attrs.keys():
        attrs["global"]["version"] = version
        attrs["global"]["creation-date"] = f"{datetime.datetime.now()}"
    else:
        attrs["version"] = version
        attrs["creation-date"] = f"{datetime.datetime.now()}"
    attrs = fix_attr_types(attrs, dtype)
    return attrs


def coord_attributes_access(
    coord: str, date: Optional[datetime.date] = None, dtype=np.float32
) -> dict[str, Any]:
    """
    Retrieves coordinate attributes for a specific coordinate and date.

    Parameters:
        coord (str): The coordinate name.
        date (Optional[datetime.date]): The date. Default is None.
        dtype: The data type. Default is np.float32.

    Returns:
        dict[str, Any]: The coordinate attributes.
    """

    attrs = load_access_attrs(var=coord)

    if coord == "hours" and date is not None:
        attrs["units"] = f"hours since {date.isoformat()} 00:00:00.0"

    attrs = fix_attr_types(attrs, dtype)
    return attrs


def attrs_as_string(attrs, prefix=""):
    """
    Converts attributes dictionary to a string representation.

    Parameters:
        attrs (dict): The attributes dictionary.
        prefix (str): The prefix string. Default is an empty string.

    Returns:
        str: The string representation of the attributes.
    """

    out_str = ""
    for key in attrs.keys():
        if isinstance(attrs[key], dict):
            out_str += f"---{key}---\n"
            s = attrs_as_string(attrs[key], prefix=f"{prefix} {key}:")
            out_str += s
        else:
            attr_as_truncated_string = str(attrs[key])[0 : 50 - len(prefix)]
            type_str = str(type(attrs[key]))
            if "numpy" in type_str:
                type_str += f": {attrs[key].dtype}"
            out_str += f"{prefix}{key:25}:  {attr_as_truncated_string}: {type_str}\n"
    return out_str


if __name__ == "__main__":
    # exercise the code

    project = "access"
    satellite = "AMSR2"
    version = "v01r00"

    AMSR2_vars = [
        "atm_pars_era5",
        "land_fraction_modis",
        "ocean_emiss_era5",
        "rain_rate_imerg",
        "resamp_tbs",
        "skt_era5",
        "tclw_era5",
        "tcwv_era5",
        "u10n_era5",
        "v10n_era5",
    ]
    access_vars = [
        "common",
        "frequency",
        "hours",
        "latitude",
        "longitude",
        "polarization",
        "time",
    ]

    print("Testing loading json files")
    for var in access_vars:
        try:
            attrs = load_attrs(project=project, var=var)
        except Exception as e:
            print(e)
            print(f"Test Failed: Could not read attr.{project}.{var}.json")
            raise

    for var in AMSR2_vars:
        try:
            attrs = load_attrs(project=project, satellite=satellite, var=var)
        except Exception as e:
            print(e)
            print(f"Test Failed: Could not read attr.{project}.{var}.json")
            raise

    print("PASSED: Testing loading json files")

    print("Testing Attribute Routines and dumping attributes")

    filename = (
        _attr_define_root / f"attributes_access_{satellite}.test.{version}.dump.txt"
    )
    print(f"dump filename: {filename}")

    with open(filename, "w") as f:
        date = datetime.date(2013, 4, 12)
        target_size = 30
        try:
            attrs = common_global_attributes_access(
                date=date, satellite=satellite, target_size=target_size, version=version
            )
            f.write(attrs_as_string(attrs, prefix="common:"))

        except Exception as e:
            print(e)
            print(
                "Test Failed: Could not load common attrs:."
                + f"{project}.{satellite}.{var}.json"
            )
            raise

        try:
            attrs = resamp_tb_attributes_access(satellite=satellite, version=version)
            f.write(attrs_as_string(attrs, prefix=attrs["name"]))
        except Exception as e:
            print(e)
            print(
                "Test Failed: Could not load resamp_tbs attrs:."
                + f"{project}.{satellite}.resamp_tbs.json"
            )
            raise

        try:
            attrs = atm_pars_era5_attributes_access(
                satellite=satellite, target_size=target_size, version=version
            )
            f.write(attrs_as_string(attrs, prefix=attrs["name"]))

        except Exception as e:
            print(e)
            print(
                "Test Failed: Could not load atm_par_era5 attrs:"
                + f".{project}.{satellite}.atm_pars_era5.json"
            )
            raise

        anc_var_names = [
            "rain_rate_imerg",
            "ocean_emiss_era5",
            "skt_era5",
            "tclw_era5",
            "tcwv_era5",
            "u10n_era5",
            "v10n_era5",
        ]
        for var in anc_var_names:
            try:
                attrs = anc_var_attributes_access(
                    satellite=satellite, var=var, version=version
                )
                f.write(attrs_as_string(attrs, prefix=attrs["name"] + ": " + var))
            except Exception as e:
                print(e)
                print(
                    "Test Failed: Could not load anc attrs:"
                    + f".{project}.{satellite}.{var}.json"
                )
                raise

        coord_var_names = [
            "frequency",
            "hours",
            "latitude",
            "longitude",
            "polarization",
            "time",
        ]

        for var in coord_var_names:
            try:
                attrs = coord_attributes_access(coord=var, date=date)
                f.write(attrs_as_string(attrs, prefix="coord: " + var))
            except Exception as e:
                print(e)
                print(
                    "Test Failed: Could not load coord attrs:"
                    + f"attr.{project}.{var}.json"
                )
                raise
    print("All Tests PASSED")
