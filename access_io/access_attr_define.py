import datetime
import json
import os
from pathlib import Path

"""
The routines in this file define the attributes for the ACCESS project output file
Constant attributes are read in from the .json files in attr_define_root
"""
if os.name == "nt":
    attr_define_root = Path("M:/job_access/python/dataset_assembly/access_io/attr_define_json")
elif os.name == "posix":
    attr_define_root = Path("/mnt/ops1p-ren/m/job_access/python/dataset_assembly/access_io/attr_define_json")

def common_global_attributes_access(
    date: datetime.datetime, satellite: str, target_size: int, version: str = "v00r00"
) -> dict:

    day_boundary = datetime.datetime.combine(date, datetime.time())
    start_date = day_boundary - datetime.timedelta(minutes=30)
    end_date = day_boundary + datetime.timedelta(minutes=1410.0)

    file = attr_define_root / f"common_attributes_access.json"

    with open(file) as json_file:
        attrs_glob = json.load(json_file)

    attrs_glob["title"] = f"Resampled {satellite} brightness temperatures"
    attrs_glob["date_created"] = datetime.datetime.now().isoformat()
    attrs_glob["spatial_resolution"] = f"{target_size} km X {target_size} km"
    attrs_glob["time_coverage_start"] = start_date.isoformat()
    attrs_glob["time_coverage_end"] = end_date.isoformat()

    return attrs_glob

def resamp_tb_attributes_access(satellite: str, version="v00r00"):
    if satellite == "AMSR2":
        file = attr_define_root / "resamp_tb_attributes_access_AMSR2.json"
    else:
        raise ValueError(
            f"Satellite {satellite} not valid in resamp_tb_attributes_access"
        )

    with open(file) as json_file:
        attrs_glob_tb = json.load(json_file)

    attrs_glob_tb["date_created"] = datetime.datetime.now().isoformat()
    attrs_glob_tb["version"] = version

    return attrs_glob_tb

def atm_pars_era5_attributes_access(satellite: str, version="v00r00"):

    file = attr_define_root / f"atm_pars_era5_attributes_access_{satellite}.json"

    with open(file) as json_file:
        attrs = json.load(json_file)

    attrs["version"] = version
    attrs["creation-date"] = f"{datetime.datetime.now()}"
    return attrs

def anc_var_attributes_access(satellite: str, var: str, version="v00r00"):

    file = attr_define_root / f"{var}_attributes_access_{satellite}.json"

    with open(file) as json_file:
        attrs_glob = json.load(json_file)

    attrs_glob["version"] = version
    attrs_glob["creation-date"] = f"{datetime.datetime.now()}"

    return attrs_glob

if __name__ == "__main__":
    # exercise the code

    def write_attrs(file,attrs):
        for key in attrs.keys():
            file.write(f"{key:25}:  {attrs[key]}\n")

    satellite = "AMSR2"
    version = "v0r01"

    filename = attr_define_root / f"attributes_access_{satellite}.test.{version}.txt"
    with open(filename,"w") as file:
        date = datetime.date(2021, 3, 12)
        target_size = 70
        version = "v100r01"
        attrs = common_global_attributes_access(date, satellite, target_size, version=version)
        file.write("--------------------------\n")
        file.write("ACCESS Common Attributes\n")
        file.write("--------------------------\n")
        write_attrs(file,attrs)

        satellite = "AMSR2"
        attrs = resamp_tb_attributes_access(satellite, version=version)
        file.write("--------------------------\n")
        file.write("ACCESS Resamp Tbs Attributes\n")
        file.write("--------------------------\n")
        write_attrs(file,attrs)

        attrs = atm_pars_era5_attributes_access(satellite, version=version)
        file.write("--------------------------\n")
        file.write("Atm Pars Attributes\n")
        file.write("--------------------------\n")
        write_attrs(file,attrs)

        var_list = [
            "v10n_era5",
            "ocean_emiss_era5",
            "skt_era5",
            "imerg_rr",
            "tcwv_era5",
            "tclw_era5",
        ]

        for var in var_list:
            attrs = anc_var_attributes_access(satellite, var, version=version)
            file.write("--------------------------")
            file.write(f"ACCESS {var} Attributes")
            file.write("--------------------------")
            write_attrs(file,attrs)

