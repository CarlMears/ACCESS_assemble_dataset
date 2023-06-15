import datetime
import os
from pathlib import Path

import xarray as xr

from access_io.access_output import get_access_output_filename_daily_folder

start_date = datetime.datetime(2012, 1, 1)
end_date = datetime.datetime(2021, 12, 31)

if os.name == "nt":
    access_root = Path("L:/access/")
elif os.name == "posix":
    access_root = Path("/mnt/ops1p-ren/l/access")


satellite = "amsr2"
date_to_do = start_date
while date_to_do <= end_date:
    old_file = get_access_output_filename_daily_folder(
        date_to_do, satellite, 70, access_root / "amsr2_out_NP_70", "resamp_tbs"
    )
    new_file = get_access_output_filename_daily_folder(
        date_to_do, satellite, 70, access_root / "amsr2_out_GL_70", "resamp_tbs"
    )

    # print(f"Old File {old_file}")
    # print(f"New File {new_file}")
    
    
    if os.path.isfile(old_file):
        os.makedirs(new_file.parent, exist_ok=True)
        if new_file.exists():
            continue
        print(f"Renaming {old_file} to")
        print(f"         {new_file}")
        os.rename(old_file, new_file)
    

    date_to_do += datetime.timedelta(days=1)
