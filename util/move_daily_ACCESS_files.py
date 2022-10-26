from access_io.access_output import (
    get_access_output_filename,
    get_access_output_filename_daily_folder,
)
import os
from pathlib import Path
import datetime
import xarray as xr


start_date = datetime.datetime(2012, 1, 1)
end_date = datetime.datetime(2012, 12, 31)

if os.name == "nt":
    access_root = Path("L:/access/amsr2_out")
elif os.name == "posix":
    access_root = Path("/mnt/ops1p-ren/laccess/amsr2_out")

satellite = "amsr2"
possible_extra_vars = ["rainfall_rate", "skt"]
date_to_do = start_date
while date_to_do <= end_date:

    old_file = get_access_output_filename(
        date_to_do, satellite, access_root, "resamp_tbs"
    )
    new_file = get_access_output_filename_daily_folder(
        date_to_do, satellite, access_root, "resamp_tbs"
    )

    print(f"Old File {old_file}")
    print(f"New File {new_file}")
    try:
        with xr.open_dataset(old_file) as ds:
            extra_vars_not_present = True
            for test_var in possible_extra_vars:
                if test_var in ds.keys():
                    extra_vars_not_present = False
                    break
    except FileNotFoundError:
        date_to_do += datetime.timedelta(days=1)
        continue

    if extra_vars_not_present:
        os.makedirs(new_file.parent, exist_ok=True)
        if new_file.exists():
            continue
        print(f"Renaming {old_file} to")
        print(f"         {new_file}")
        os.rename(old_file, new_file)
    else:
        print(f"Need to remake {old_file}")

    print

    date_to_do += datetime.timedelta(days=1)
