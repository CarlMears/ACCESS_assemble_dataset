from pathlib import Path
import numpy as np
import datetime
import matplotlib.pyplot as plt
from rss_plotting.plot_2d_array import plot_2d_array


def is_file_multiple(path_to_test: Path, max_tries: int = 10):
    num_tries = 0
    is_file_multiple = False
    while num_tries < max_tries:
        try:
            is_file_multiple = path_to_test.is_file()
            break
        except OSError:
            print(f"{num_tries}: path_to_test")
            pass
    return is_file_multiple


def radiance_file_exists(date: datetime.date, output_root: Path):
    file = output_root / f"era5_tbs_{date.year:04d}-{date.month:02d}-{date.day:02d}.nc"

    return is_file_multiple(file)


def plot_summary(
    start_date: datetime.date, end_date: datetime.date, file_exists_all: np.ndarray
):
    import matplotlib.dates as mdates

    yvals = np.arange(2)
    xvals = np.arange((end_date - start_date).days + 1)
    days = mdates.drange(
        start_date, end_date + datetime.timedelta(days=1), datetime.timedelta(days=1)
    )
    xvals = days

    fig = plt.figure(figsize=(10, 5))
    fig, ax = plot_2d_array(
        file_exists_all, xvals, yvals, zrange=(0.0, 1.0), fig_in=fig, cmap="Greens"
    )
    ax.set_title("Simulated Radiance Inventory")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax.tick_params(axis="x", labelrotation=90)
    plt.subplots_adjust(left=0.3, bottom=0.25)

    return fig, ax


if __name__ == "__main__":
    import os

    if os.name == "nt":
        output_root = Path("A:/_access_temp/rtm/tbs")
    elif os.name == "posix":
        output_root = Path("/mnt/A/_access_temp/rtm/tbs)")
    else:
        raise ValueError

    for year in range(2021, 2022):
        start_date = datetime.date(year, 1, 1)
        end_date = datetime.date(year, 12, 31)

        num_days = end_date - start_date
        file_exist_all = np.zeros((2, num_days.days + 1), dtype=np.int32)

        date_to_do = start_date
        date_index = 0
        while date_to_do <= end_date:
            exists = radiance_file_exists(date_to_do, output_root)
            if exists:
                file_exist_all[0, date_index] = 1
                file_exist_all[1, date_index] = 1
            date_to_do += datetime.timedelta(days=1)
            date_index += 1

        fig, ax = plot_summary(start_date, end_date, file_exist_all)
        summary_path = output_root / "summaries"
        png_file = summary_path / f"dataset_summary_{year:04d}.png"
        os.makedirs(png_file.parent, exist_ok=True)
        fig.savefig(png_file)
        plt.show()
        print()
