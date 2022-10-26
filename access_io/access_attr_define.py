import datetime
import json
import numpy as np
import os
from pathlib import Path
from typing import Callable

"""
The routines in this file define the attributes for the ACCESS project output file
Constant attributes are read in from the .json files in attr_define_root
"""
if os.name == "nt":
    attr_define_root = Path(
        "M:/job_access/python/dataset_assembly/access_io/attr_define_json"
    )
elif os.name == "posix":
    attr_define_root = Path(
        "/mnt/ops1p-ren/m/job_access/python/dataset_assembly/access_io/attr_define_json"
    )


def _convert_attrs_to_numbers(
    attrs: dict,
    int_converter: Callable = np.int32,
    float_converter: Callable = np.float32,
    keys_to_exclude: list[str] = ["key_to_exclude"],
):
    for key in attrs.keys():
        if key in keys_to_exclude:
            continue
        value = attrs[key]
        try:
            _ = attrs[key].keys()  # determine if dictionary-like, e.g. has keys()
            attrs[key] = _convert_attrs_to_numbers(value)
        except AttributeError:
            try:
                value_int = int_converter(value)
                attrs[key] = value_int
            except ValueError:  # maybe too big for int32 - try int64
                try:
                    value_int = np.int64(value)
                    attrs[key] = value_int
                except ValueError:
                    try:
                        value_flt = float_converter(value)
                        attrs[key] = value_flt
                    except ValueError:
                        pass

    return attrs


def load_attrs(project: str = "", satellite: str = "", var: str = "") -> dict:

    if len(var) == 0:
        raise ValueError("var must be specified")

    filename = "attr."
    if len(project) > 0:
        filename = f"{filename}{project}."
    if len(satellite) > 0:
        filename = f"{filename}{satellite.upper()}."
    filename = f"{filename}{var}.json"
    path_to_file = attr_define_root / filename
    with open(path_to_file) as json_file:
        attrs = json.load(json_file)

    return attrs


def load_access_attrs(satellite: str = "", var: str = "") -> dict:

    project = "access"

    return load_attrs(project=project, satellite=satellite, var=var)


def common_global_attributes_access(
    date: datetime.datetime, satellite: str, target_size: int, version: str = "v00r00"
) -> dict:

    attrs = load_access_attrs(var="common")

    day_boundary = datetime.datetime.combine(date, datetime.time())
    start_date = day_boundary - datetime.timedelta(minutes=30)
    end_date = day_boundary + datetime.timedelta(minutes=1410.0)

    attrs["title"] = f"Resampled {satellite} brightness temperatures"
    attrs["date_created"] = datetime.datetime.now().isoformat()
    attrs["spatial_resolution"] = f"{target_size} km X {target_size} km"
    attrs["time_coverage_start"] = start_date.isoformat()
    attrs["time_coverage_end"] = end_date.isoformat()

    attrs = _convert_attrs_to_numbers(attrs)

    return attrs


def resamp_tb_attributes_access(satellite: str, version="v00r00"):
    attrs = load_access_attrs(satellite=satellite, var="resamp_tbs")
    attrs["global"]["date_created"] = datetime.datetime.now().isoformat()
    attrs["global"]["version"] = version

    attrs = _convert_attrs_to_numbers(attrs)

    return attrs


def atm_pars_era5_attributes_access(satellite: str, target_size: int, version="v00r00"):

    attrs = load_access_attrs(satellite=satellite, var="atm_pars_era5")

    attrs["global"]["version"] = version
    attrs["global"]["date_accessed"] = f"{datetime.datetime.now()}"

    attrs = _convert_attrs_to_numbers(attrs)

    return attrs


def anc_var_attributes_access(satellite: str, var: str, version="v00r00"):

    attrs = load_access_attrs(satellite=satellite, var=var)

    if "global" in attrs.keys():
        attrs["global"]["version"] = version
        attrs["global"]["creation-date"] = f"{datetime.datetime.now()}"
    else:
        attrs["version"] = version
        attrs["creation-date"] = f"{datetime.datetime.now()}"

    attrs = _convert_attrs_to_numbers(attrs)

    return attrs


def coord_attributes_access(coord: str, date=None):

    attrs = load_access_attrs(var=coord)

    if coord == "hours":
        attrs["units"] = f"hours since {date.isoformat()} 00:00:00.0"

    if coord == 'time':
        attrs = _convert_attrs_to_numbers(attrs,int_converter=np.int64)
    else:
        attrs = _convert_attrs_to_numbers(attrs)

    return attrs


def write_attrs(file, attrs, prefix=""):
    for key in attrs.keys():
        if isinstance(attrs[key], dict):
            if file == "screen":
                print(f"{prefix}---{key}---")
            else:
                file.write(f"---{key}---\n")
            write_attrs(file, attrs[key], prefix=f"{prefix} {key}:")
        else:
            attr_as_truncated_string = str(attrs[key])[0 : 50 - len(prefix)]
            if file == "screen":
                print(f"{prefix}{key:25}:  {attr_as_truncated_string}")
            else:
                file.write(f"{prefix}{key:25}:  {attr_as_truncated_string}\n")


if __name__ == "__main__":
    # exercise the code

    satellite = "AMSR2"
    version = "v0r01"

    filename = attr_define_root / f"attributes_access_{satellite}.test.{version}.txt"
    with open(filename, "w") as file:
        date = datetime.date(2021, 3, 12)
        target_size = 70
        version = "v100r01"
        attrs = common_global_attributes_access(
            date, satellite, target_size, version=version
        )
        file.write("--------------------------\n")
        file.write("ACCESS Common Attributes\n")
        file.write("--------------------------\n")
        write_attrs(file, attrs)

        satellite = "AMSR2"
        attrs = resamp_tb_attributes_access(satellite, version=version)
        file.write("--------------------------\n")
        file.write("ACCESS Resamp Tbs Attributes\n")
        file.write("--------------------------\n")
        write_attrs(file, attrs)

        attrs = atm_pars_era5_attributes_access(satellite, version=version)
        file.write("--------------------------\n")
        file.write("Atm Pars Attributes\n")
        file.write("--------------------------\n")
        write_attrs(file, attrs)

        var_list = [
            "v10n_era5",
            "ocean_emiss_era5",
            "skt_era5",
            "imerg_rr",
            "tcwv_era5",
            "tclw_era5",
        ]

        glb_required_attributes = ["summary", "keywords", "title", "history", "source"]
        var_required_attributes = [
            "standard_name",
            "long_name",
            "valid_min",
            "valid_max",
            "units",
            "_FillValue",
            "missing_value",
        ]

        for var in var_list:
            attrs = anc_var_attributes_access(satellite, var, version=version)
            file.write("--------------------------")
            file.write(f"ACCESS {var} Attributes")
            file.write("--------------------------")
            write_attrs(file, attrs)

        print("-----------------------------------------------------")
        print("-- Testing for Required Attributes in .json files ---")
        print("-----------------------------------------------------")
        for var in var_list:
            attrs = anc_var_attributes_access(satellite, var, version=version)
            for required_attr in glb_required_attributes:
                if required_attr not in attrs["global"].keys():
                    print(f"required attribite: {required_attr} missing for {var}")

            for required_attr in var_required_attributes:
                if required_attr not in attrs["var"].keys():
                    print(f"required attribite: {required_attr} missing for {var}")
