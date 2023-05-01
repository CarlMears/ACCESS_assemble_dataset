import argparse
from contextlib import suppress
import datetime
import git
import os
from pathlib import Path
from typing import Tuple, Union
from netCDF4 import Dataset as netcdf_dataset
import numpy as np

from era5_request.era5_requests import era5_hourly_single_level_request
from access_io.access_output import get_access_output_filename_daily_folder
from access_io.access_output import write_daily_ancillary_var_netcdf


output_root = Path("N:/data/model/ERA5/hourly")
var_list = ["tcwv", "u10n", "v10n", "tclw", "skt"]
year = 2017

for var in var_list:
    for month in range(1, 13):

        var_dict = {
            "skt": "Skin temperature",
            "tcwv": "Total column water vapour",
            "tclw": "total_column_cloud_liquid_water",
            "u10n": "10m_u_component_of_neutral_wind",
            "v10n": "10m_v_component_of_neutral_wind",
        }
        try:
            variable = (var, var_dict[var])
        except KeyError:
            print(f"Variable {var} not defined - skipping")
            continue

        os.makedirs(output_root, exist_ok=True)

        current_day = datetime.datetime(year, month, 15)
        file1 = era5_hourly_single_level_request(
            date=current_day,
            variable=variable[1],
            target_path=output_root,
            full_day=True,
            full_month=True,
            simpler_path=True,
        )
