from contextlib import suppress
from datetime import date
from pathlib import Path
from typing import Collection, List
from copy import copy

import git
import matplotlib.pyplot as plt
import numpy as np

# must be installed from rss_plotting package
from rss_plotting.global_map import plot_global_map

# these packages are located in folders in the local path
from access_io.access_output import (
    write_daily_tb_netcdf,
    edit_attrs_daily_tb_netcdf,
    get_access_output_filename_daily_folder,
    write_daily_tb_netcdf,
)
from access_io.access_output_polar import write_daily_tb_netcdf_polar
from resampled_tbs.read_resampled_orbit import (
    read_AMSR2_resampled_tbs,
    get_resampled_file_name,
)
from util.numpy_date_utils import convert_to_sec_in_day
from util.orbit_times_amsr2 import find_orbits_in_day, read_amsr2_orbit_times
from util.file_times import need_to_process_base_file

NUM_LATS = 721
NUM_LONS = 1440
NUM_X = 720
NUM_Y = 720
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


def get_mtime_multi_try(path_to_file, max_num_trys=10):

    num_so_far = 0
    while num_so_far < max_num_trys:
        try:
            mtime = path_to_file.stat().st_mtime
            return mtime
        except Exception:
            # catch all errors
            num_so_far += 1

    raise IOError(f"Problem getting mtime for {path_to_file}")


def decide_not_to_process(
    *,
    filename: Path,
    satellite: str,
    channels: list,
    target_size: int,
    orbits_to_do: list[int],
    overwrite: bool,
    update: bool,
):

    """Logic to determine if a daily Tb file needs to be reprocessed
    If the daily file does not exist:
         process
    If the daily file exists:
         If overwrite -> process
         If update:
             If any of the Tb orbit files (or the time file) are newer than
             the daily file -> process
         Else:
             don't process
         If neither overwrite or update are True:
             don't process
    """

    if filename.is_file():
        if overwrite:
            with suppress(FileNotFoundError):
                filename.unlink()
            return False  # force overwrite
        else:
            if update:
                tb_day_file_time = filename.stat().st_mtime

                # check to see if any of the included tb orbit files are newer than
                # the daily tb file

                need_to_update = False
                for orbit in orbits_to_do:
                    channels_to_do = copy(channels)
                    channels_to_do.append("time")
                    for channel in channels_to_do:
                        tb_orbit_file = get_resampled_file_name(
                            satellite=satellite,
                            channel=channel,
                            target_size=target_size,
                            orbit=orbit,
                        )
                        if tb_orbit_file.is_file():
                            tb_file_time = tb_orbit_file.stat().st_mtime
                            if tb_file_time > tb_day_file_time:
                                need_to_update = True
                        else:
                            print(f"Warning - Tb file {tb_orbit_file} is missing")

                if need_to_update:
                    with suppress(FileNotFoundError):
                        filename.unlink()
                    return (
                        False  # At least 1 Tb orbit files are newer than the daily file
                    )
                else:
                    return True  # Up to date
            else:
                return True  # don't update because update not set
    else:
        return False  # daily file does not exist


def redo_attrs_daily_ACCESS_tb_file(
    *,
    current_day: date,
    satellite: str,
    target_size: int,
    version: str,
    dataroot: Path,
    channels: Collection[int],
    verbose: bool = False,
    plot_example_map: bool = True,
    overwrite: bool = False,
    script_name: str = "unavailable",
    commit: str = "unavailable",
) -> None:

    filename = get_access_output_filename_daily_folder(
        current_day, satellite, target_size, dataroot, "resamp_tbs"
    )
    if filename.is_file():
        edit_attrs_daily_tb_netcdf(
            date=current_day,
            satellite=satellite,
            target_size=target_size,
            version=None,
            tb_array_by_hour=None,
            time_array_by_hour=None,
            dataroot=dataroot,
            freq_list=None,
            file_list=None,
            script_name=None,
            commit=None,
        )


def make_daily_ACCESS_tb_file(
    *,
    current_day: date,
    satellite: str,
    target_size: int,
    region: str = "global",
    version: str,
    dataroot: Path,
    channels: Collection[int],
    verbose: bool = False,
    plot_example_map: bool = True,
    overwrite: bool = False,
    update: bool = False,
    script_name: str = "unavailable",
    commit: str = "unavailable",
) -> List[Path]:
    if satellite.lower() == "amsr2":
        orbit_times = read_amsr2_orbit_times()
        from satellite_definitions.amsr2 import (
            CHANNEL_TO_FREQ_MAP,
            CHANNEL_TO_POL_MAP,
            REF_FREQ,
            SAT_NAME,
        )

        assert SAT_NAME.lower() == satellite.lower()
        NUM_FREQS = len(REF_FREQ)
        NUM_POLS = 2
    else:
        raise ValueError(f"No satellite definition file for {satellite}")

    if region == "global":
        filename = get_access_output_filename_daily_folder(
            current_day, satellite, target_size, dataroot, "resamp_tbs"
        )
        grid_type = "equirectangular"
        pole = None
    elif region in ["north", "south"]:
        filename = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            target_size,
            dataroot,
            "resamp_tbs",
            grid_type="ease2",
            pole="north",
        )
        grid_type = "ease2"
        pole = region
    else:
        raise ValueError(f"region {region} not valid")

    orbits_to_do = find_orbits_in_day(times_np64=orbit_times, date=current_day)
    # TODO: what is the actual exception expected? AssertionError?

    if len(orbits_to_do) == 0:
        print(f"No orbits found for {current_day}")
        return []

    # if decide_not_to_process(
    #     filename=filename,
    #     satellite=satellite,
    #     channels=channels,
    #     target_size=target_size,
    #     orbits_to_do=orbits_to_do,
    #     overwrite=overwrite,
    #     update=update,
    # ):
    #     print(f"No processing needed for {satellite} base file on {current_day}")
    #     return []

    need_to_process = need_to_process_base_file(
        date=current_day,
        satellite=satellite,
        target_size=target_size,
        grid_type=grid_type,
        pole=pole,
        orbits_to_do=orbits_to_do,
        channels_to_do=channels,
        dataroot=dataroot,
        outputroot=dataroot,
        var="resamp_tbs",
        overwrite=overwrite,
        update=update,
    )

    if need_to_process is False:
        print(f"No processing needed for {satellite} base file on {current_day}")
        return []

    # initialize arrays for daily data
    at_least_one_orbit = False
    file_list = []
    if region == "global":
        tb_array_by_hour = np.full(
            (NUM_LATS, NUM_LONS, NUM_HOURS, NUM_FREQS, NUM_POLS),
            np.nan,
            dtype=np.float32,
        )
        time_array_by_hour = np.full((NUM_LATS, NUM_LONS, NUM_HOURS), np.nan)

    elif region in ["north", "south"]:
        tb_array_by_hour = np.full(
            (NUM_X, NUM_Y, NUM_HOURS, NUM_FREQS, NUM_POLS), np.nan, dtype=np.float32
        )
        time_array_by_hour = np.full((NUM_X, NUM_Y, NUM_HOURS), np.nan)
    else:
        raise ValueError(f"Region {region} is not valid")

    print(f"Processing {current_day:%Y/%m/%d}, orbit: ", end="")
    for orbit in orbits_to_do:
        print(f"{orbit} ", end="")

        # get the times for this orbit, and convert to time within the day...
        try:
            ob_time, filename = read_AMSR2_resampled_tbs(
                satellite=satellite,
                channel="time",
                target_size=target_size,
                orbit=orbit,
                grid_type=grid_type,
                pole=pole,
                verbose=False,
            )
        except FileNotFoundError:
            print(f"Time file corrupt or missing for orbit: {orbit}, - skipping orbit")
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
                tbs, filename = read_AMSR2_resampled_tbs(
                    satellite=satellite,
                    channel=channel,
                    target_size=target_size,
                    orbit=orbit,
                    grid_type=grid_type,
                    pole=pole,
                )
            # There are a lot of possible errors for corrupted files
            except:
                print(
                    f"File corrupt or missing for orbit: {orbit}, channel: {channel} - skipping channel"
                )
                continue

            file_list.append(filename)
            at_least_one_orbit = True
            # loop through the hours, saving the tbs and times in hour-long slices
            for hour in range(0, 24):
                freq_index = CHANNEL_TO_FREQ_MAP[channel - 1]
                pol_index = CHANNEL_TO_POL_MAP[channel - 1]
                tb_slice = tb_array_by_hour[:, :, hour, freq_index, pol_index]
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
    print(" ")  # force next line

    if at_least_one_orbit:
        if plot_example_map:
            plot_global_map(tb_array_by_hour[:, :, 0, 5], vmin=0.0, vmax=330)

        if region == "global":
            write_daily_tb_netcdf(
                date=current_day,
                satellite=satellite,
                target_size=target_size,
                version=version,
                tb_array_by_hour=tb_array_by_hour,
                time_array_by_hour=time_array_by_hour,
                freq_list=REF_FREQ,
                file_list=file_list,
                dataroot=dataroot,
                script_name=script_name,
                commit=commit,
            )
        elif region in ["north", "south"]:
            write_daily_tb_netcdf_polar(
                date=current_day,
                satellite=satellite,
                target_size=target_size,
                pole=pole,
                version=version,
                tb_array_by_hour=tb_array_by_hour,
                time_array_by_hour=time_array_by_hour,
                freq_list=REF_FREQ,
                file_list=file_list,
                dataroot=dataroot,
                script_name=script_name,
                commit=commit,
            )
        else:
            raise ValueError(f"Region {region} is not valid")
    return file_list


if __name__ == "__main__":
    import argparse
    import datetime

    parser = argparse.ArgumentParser(
        description=("Arrange resampled TBs into an ACCESS output file. ")
    )
    parser.add_argument(
        "--access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument(
        "--temp_root", type=Path, help="Root directory store temporary files"
    )
    parser.add_argument(
        "--start_date",
        type=datetime.date.fromisoformat,
        help="First Day to process, as YYYY-MM-DD",
    )
    parser.add_argument(
        "--end_date",
        type=datetime.date.fromisoformat,
        help="Last Day to process, as YYYY-MM-DD",
    )
    parser.add_argument("--sensor", choices=["amsr2"], help="Microwave sensor to use")
    parser.add_argument(
        "--target_size", choices=["30", "70"], help="Size of target footprint in km"
    )
    parser.add_argument(
        "--region", choices=["global", "north", "south"], default="global"
    )
    parser.add_argument("--version", help="version sting - e.g. v01r00")
    parser.add_argument(
        "--overwrite", help="force overwrite if file exists", action="store_true"
    )
    parser.add_argument(
        "--update",
        help="force overwrite if file is older than orbit files",
        action="store_true",
    )
    parser.add_argument("--plot_map", help="plot an example map", action="store_true")
    parser.add_argument(
        "--verbose", help="enable more verbose screen output", action="store_true"
    )
    parser.add_argument(
        "--redo_attrs",
        help="rewrite attributes using the current .json files",
        action="store_true",
    )

    args = parser.parse_args()

    access_root: Path = args.access_root
    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()
    target_size = int(args.target_size)

    if target_size == 30:
        channels = list(range(5, 13))
    else:
        channels = list(range(1, 13))

    script_name = parser.prog
    repo = git.Repo(search_parent_directories=True)
    commit = repo.head.object.hexsha

    day_to_do = START_DAY
    while day_to_do <= END_DAY:
        if args.redo_attrs:
            redo_attrs_daily_ACCESS_tb_file(
                current_day=day_to_do,
                satellite=satellite,
                target_size=target_size,
                region=args.region,
                version=args.version,
                dataroot=access_root,
                channels=channels,
                verbose=args.verbose,
                plot_example_map=args.plot_map,
                overwrite=args.overwrite,
                update=args.update,
                script_name=script_name,
                commit=commit,
            )
        else:
            make_daily_ACCESS_tb_file(
                current_day=day_to_do,
                satellite=satellite,
                target_size=target_size,
                region=args.region,
                version=args.version,
                dataroot=access_root,
                channels=channels,
                verbose=args.verbose,
                plot_example_map=args.plot_map,
                overwrite=args.overwrite,
                update=args.update,
                script_name=script_name,
                commit=commit,
            )
            if args.plot_map:
                plt.show()
                
        day_to_do += datetime.timedelta(days=1)
