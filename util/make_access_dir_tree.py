from pathlib import Path
import os
import datetime

def make_access_dir_tree(base_dir, start_date, end_date):
    date = start_date
    while date <= end_date:
        dir_to_make = base_dir / f'Y{date.year}' / f'M{date.month:02d}' / f'D{date.day:02d}'
        os.makedirs(dir_to_make, exist_ok=True)
        date += datetime.timedelta(days=1)


if __name__ == "__main__":
    start_date = datetime.datetime(2012, 7, 1)
    end_date = datetime.datetime(2021, 12, 31)
    if os.name == 'nt':
        base_dir = Path("L:/access/amsr2_out_SP_30")
    else:
        base_dir = Path("/mnt/l/access/amsr2_out_SP_30")

    make_access_dir_tree(base_dir, start_date, end_date)

    if os.name == 'nt':
        base_dir = Path("L:/access/amsr2_out_SP_70")
    else:
        base_dir = Path("/mnt/l/access/amsr2_out_SP_70")
    make_access_dir_tree(base_dir, start_date, end_date)