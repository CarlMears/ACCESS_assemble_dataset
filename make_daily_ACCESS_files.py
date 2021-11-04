import os
from datetime import date
from pathlib import Path
from typing import Collection, List

import matplotlib.pyplot as plt
import numpy as np

# must be installed from rss_plotting package
from rss_plotting.global_map import global_map

from access_io.access_output import write_daily_tb_netcdf
from resampled_tbs.read_resampled_orbit import read_resampled_tbs
from util.numpy_date_utils import convert_to_sec_in_day
from util.orbit_times_amsr2 import find_orbits_in_day, read_amsr2_orbit_times

NUM_LATS = 721
NUM_LONS = 1440
NUM_HOURS = 24
NUM_CHANNELS = 14  # all possible AMSR2 channels
AVAILABLE_CHANNELS = [
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


def make_daily_ACCESS_tb_file(
    *,
    current_day: date,
    satellite: str,
    dataroot: Path,
    channels: Collection[int],
    verbose: bool = False,
    plot_example_map: bool = True,
) -> List[Path]:
    if satellite == "amsr2":
        orbit_times = read_amsr2_orbit_times()
    else:
        raise ValueError(f"Orbit Times for {satellite} not implemented yet")

    # initialize arrays for daily data
    at_least_one_orbit = False
    tb_array_by_hour = np.full(
        (NUM_LATS, NUM_LONS, NUM_HOURS, NUM_CHANNELS), np.nan, dtype=np.float32
    )
    time_array_by_hour = np.full((NUM_LATS, NUM_LONS, NUM_HOURS), np.nan)

    file_list = []
    orbits_to_do = find_orbits_in_day(times_np64=orbit_times, date=current_day)
    print(f"Processing {current_day:%Y/%m/%d}, orbit: ", end="")
    for orbit in orbits_to_do:
        print(f"{orbit} ", end="")
        # get the times for this orbit, and convert to time within the day...
        try:
            ob_time, filename = read_resampled_tbs(
                satellite=satellite, channel="time", orbit=orbit, verbose=False
            )
        except FileNotFoundError:
            print(f"No time file found for orbit: {orbit}")
            continue

        file_list.append(filename)
        obtime_in_day = convert_to_sec_in_day(ob_time, current_day)
        for hour in range(0, 24):
            time_slice = time_array_by_hour[:, :, hour]
            start_time_sec = hour * 3600.0
            end_time_sec = start_time_sec + 3600.0
            ok = np.all(
                [(obtime_in_day > start_time_sec), (obtime_in_day <= end_time_sec)],
                axis=0,
            )
            if np.any(ok):
                time_slice[ok] = obtime_in_day[ok]

        if verbose:
            print("reading resampled tb files")
        for channel in channels:
            try:
                tbs, filename = read_resampled_tbs(
                    satellite=satellite, channel=channel, orbit=orbit
                )
            except FileNotFoundError:
                print(f"No file found for orbit: {orbit}, channel: {channel}")
                continue
            file_list.append(filename)
            at_least_one_orbit = True
            # loop through the hours, saving the tbs and times in hour-long slices
            for hour in range(0, 24):
                tb_slice = tb_array_by_hour[:, :, hour, channel - 1]
                start_time_sec = hour * 3600.0
                end_time_sec = start_time_sec + 3600.0
                ok = np.all(
                    [(obtime_in_day > start_time_sec), (obtime_in_day <= end_time_sec)],
                    axis=0,
                )
                num_ok = np.sum(ok)
                if num_ok > 0:
                    tb_slice[ok] = tbs[ok]
                    if verbose:
                        print(
                            f"orbit:{orbit}, "
                            f"channel:{AVAILABLE_CHANNELS[channel]}, "
                            f"time_range({start_time_sec}:{end_time_sec}) "
                            f"Num Obs: {num_ok}"
                        )
    print()
    if at_least_one_orbit:
        if plot_example_map:
            global_map(tb_array_by_hour[:, :, 0, 5], vmin=0.0, vmax=330)

        write_daily_tb_netcdf(
            date=current_day,
            satellite=satellite,
            tb_array_by_hour=tb_array_by_hour,
            time_array_by_hour=time_array_by_hour,
            file_list=file_list,
            dataroot=dataroot,
        )

    return file_list


if __name__ == "__main__":
    year = 2012
    month = 7
    channels = list(range(5, 13))
    satellite = "amsr2"
    if os.name == "nt":
        dataroot = Path("L:/access")
    elif os.name == "posix":
        dataroot = Path("/mnt/ops1p-ren/l/access")

    for day in range(12, 32):
        make_daily_ACCESS_tb_file(
            current_day=date(year, month, day),
            satellite=satellite,
            dataroot=dataroot,
            channels=channels,
            verbose=False,
            plot_example_map=True,
        )

    plt.show()
