from datetime import date
from pathlib import Path
from typing import Collection, List

import matplotlib.pyplot as plt
import numpy as np

# must be installed from rss_plotting package
from rss_plotting.global_map import plot_global_map

# these packages are located in folders in the local path
from access_io.access_output import write_daily_tb_netcdf, get_access_output_filename
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
    overwrite: bool = False,
) -> List[Path]:
    if satellite.lower() == "amsr2":
        orbit_times = read_amsr2_orbit_times()
    else:
        raise ValueError(f"Orbit Times for {satellite} not implemented yet")

    filename = get_access_output_filename(current_day, satellite, dataroot)
    if filename.is_file() and not overwrite:
        print(f"daily file for {current_day} exists... skipping")
        return []

    # initialize arrays for daily data
    at_least_one_orbit = False
    tb_array_by_hour = np.full(
        (NUM_LATS, NUM_LONS, NUM_HOURS, NUM_CHANNELS), np.nan, dtype=np.float32
    )
    time_array_by_hour = np.full((NUM_LATS, NUM_LONS, NUM_HOURS), np.nan)

    file_list = []
    try:
        orbits_to_do = find_orbits_in_day(times_np64=orbit_times, date=current_day)
    # TODO: what is the actual exception expected? AssertionError?
    except Exception:
        print(f"No orbits found for {current_day}")
        return []
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
            plot_global_map(tb_array_by_hour[:, :, 0, 5], vmin=0.0, vmax=330)

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

    import argparse
    import datetime

    parser = argparse.ArgumentParser(
        description=("Arrange resampled TBs into an ACCESS output file. ")
    )
    parser.add_argument(
        "access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument(
        "temp_root", type=Path, help="Root directory store temporary files"
    )
    parser.add_argument(
        "start_date",
        type=datetime.date.fromisoformat,
        help="First Day to process, as YYYY-MM-DD",
    )
    parser.add_argument(
        "end_date",
        type=datetime.date.fromisoformat,
        help="Last Day to process, as YYYY-MM-DD",
    )
    parser.add_argument("sensor", choices=["amsr2"], help="Microwave sensor to use")
    parser.add_argument(
        "--overwrite", help="force overwrite if file exists", action="store_true"
    )
    parser.add_argument("--plot_map", help="plot an example map", action="store_true")
    parser.add_argument(
        "--verbose", help="enable more verbose screen output", action="store_true"
    )

    args = parser.parse_args()

    access_root: Path = args.access_root
    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()
    channels = list(range(5, 13))

    day_to_do = START_DAY
    while day_to_do <= END_DAY:
        print(f"{day_to_do}")
        make_daily_ACCESS_tb_file(
            current_day=day_to_do,
            satellite=satellite,
            dataroot=access_root,
            channels=channels,
            verbose=args.verbose,
            plot_example_map=args.plot_map,
            overwrite=args.overwrite,
        )
        if args.plot_map:
            plt.show()

        day_to_do += datetime.timedelta(days=1)
