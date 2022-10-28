from datetime import date
from pathlib import Path
from typing import Any, Collection

import numpy as np
import pandas as pd
import xarray as xr

# these packages are located in folders in the local path
from access_io.access_output import get_access_output_filename_daily_folder

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


def inventory_daily_ACCESS_tb_file(
    *,
    current_day: date,
    satellite: str,
    dataroot: Path,
    channels: Collection[int],
    verbose: bool = False,
) -> list[Any]:

    filename = get_access_output_filename_daily_folder(
        current_day, satellite, 0, dataroot, "UNKNOWN"
    )
    if filename.is_file():
        ds = xr.open_dataset(filename)
        return list(ds.data_vars.keys())
    else:
        return []


def plot_inventory(access_inventory):

    from rss_plotting.plot_2d_array import plot_2d_array

    inv_np = np.transpose(access_inventory.to_numpy())

    print(inv_np.shape)

    time = access_inventory.index

    vars = list(access_inventory.columns)

    fig, ax = plot_2d_array(
        inv_np,
        time,
        np.arange(0, len(vars)),
        zrange=[0.0, 1.2],
        cmap="Greens",
        plt_colorbar=False,
        figsize=(13, 5),
    )

    for y in np.arange(0, len(vars) - 1):
        ax.plot(time, np.full((len(time)), y + 0.5), color="grey")
    ax.set_yticks(np.arange(0, len(vars)))
    ax.set_yticklabels(vars)
    fig.subplots_adjust(left=0.2)
    print

    return fig, ax


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

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()
    channels = list(range(5, 13))

    day_to_do = START_DAY

    date_rng = pd.date_range(start=START_DAY, end=END_DAY, freq="D")
    df = pd.DataFrame(index=date_rng)

    while day_to_do <= END_DAY:

        list_of_vars = inventory_daily_ACCESS_tb_file(
            current_day=day_to_do,
            satellite=satellite,
            dataroot=access_root,
            channels=channels,
            verbose=args.verbose,
        )
        print(f"{day_to_do}: {len(list_of_vars)}")
        for var in list_of_vars:
            timestamp = pd.Timestamp(day_to_do)
            z = df.loc[timestamp]
            if var in df.columns:
                df.loc[timestamp][var] = 1.0
            else:
                # create new column
                df[var] = np.zeros((len(date_rng)))
                df.loc[timestamp][var] = 1.0

        day_to_do += datetime.timedelta(days=1)

    fig, ax = plot_inventory(df)
    png_file = args.access_root / "inventory.png"
    fig.savefig(png_file)
    print
