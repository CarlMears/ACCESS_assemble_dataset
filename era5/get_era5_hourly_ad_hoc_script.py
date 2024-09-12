import datetime
from era5_requests.era5_requests import era5_hourly_single_level_request
import os
from pathlib import Path

# variable = ("total_column_cloud_liquid_water", "tclw")

variable = ("Skin temperature","skt")
temproot = Path("//flux/cm/ERA5/1hr/skt")


start_day = datetime.date(2024,6,17)
end_day = datetime.date(2024,6,30)
os.makedirs(temproot, exist_ok=True)
day_to_do = start_day
while day_to_do <= end_day:
    print(day_to_do)
    try:
        file1 = era5_hourly_single_level_request(
            date=day_to_do,
            variable=variable[0],
            target_path=temproot,
            full_day=True,
            full_month=False,
            verbose=True
        )
    except FileExistsError:
        print(f'Download for {day_to_do} failed ... skipping')
    day_to_do = day_to_do + datetime.timedelta(days=1)
