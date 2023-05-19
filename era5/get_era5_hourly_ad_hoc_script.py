import datetime
from Era5_requests.era5_requests import era5_hourly_single_level_request
import os
from pathlib import Path

# variable = ("total_column_cloud_liquid_water", "tclw")
variable = ("Total column water vapour", "tcwv")
temproot = Path("L:/era5_vapor/")

current_day = datetime.date(2022, 9, 24)
os.makedirs(temproot, exist_ok=True)
for iday in range(0, 14):
    print(current_day)
    file1 = era5_hourly_single_level_request(
        date=current_day,
        variable=variable[0],
        target_path=temproot,
        full_day=True,
        full_month=False,
    )
    current_day = current_day + datetime.timedelta(days=1)
