from pathlib import Path
import numpy as np
import datetime
import matplotlib.pyplot as plt
from rss_plotting.plot_2d_array import plot_2d_array


def inventory_access_output_files(
    output_root: Path, date_to_do: datetime.date, footprint_size: int, base_template:Path, templates: list[Path]
) -> tuple:

    debug = True
    output_root_this_day = (
        output_root
        / f"Y{date_to_do.year:04d}"
        / f"M{date_to_do.month:02d}"
        / f"D{date_to_do.day:02d}"
    )
    list_of_missing_files = []
    exists = np.zeros(len(templates)+1, np.int32)

    date_str = f"{date_to_do.year:04d}_{date_to_do.month:02d}_{date_to_do.day:02d}"
    file_name = f"{base_template}_{date_str}.{footprint_size:03d}km.nc"
    path_to_file = output_root_this_day / file_name

    if not path_to_file.is_file():
        list_of_missing_files.append(file_name)
        exists[0] = 0
        if debug:
            print(f"Base File Missing: {path_to_file}")
    else:
        exists[0] = 2
        base_file_time = path_to_file.stat().st_mtime


    if exists[0] > 0:
        for file_index, file_template in enumerate(templates):
            date_str = f"{date_to_do.year:04d}_{date_to_do.month:02d}_{date_to_do.day:02d}"
            file_name = f"{file_template}_{date_str}.{footprint_size:03d}km.nc"
            path_to_file = output_root_this_day / file_name

            if not path_to_file.is_file():
                list_of_missing_files.append(file_name)
                exists[file_index+1] = 0
                if debug:
                    print(f"Missing: {path_to_file}")
            else:
                file_time = path_to_file.stat().st_mtime
                if file_time > base_file_time:
                    exists[file_index+1] = 2
                else:
                    exists[file_index+1] = 1

    return exists, list_of_missing_files


def plot_summary(
    start_date: datetime.date,
    end_date: datetime.date,
    file_exists_all: np.ndarray,
    base_template: Path,
    templates: list[Path],
):

    import matplotlib.dates as mdates

    yvals = np.arange(len(templates)+1)
    xvals = np.arange((end_date - start_date).days + 1)
    days = mdates.drange(
        start_date, end_date + datetime.timedelta(days=1), datetime.timedelta(days=1)
    )
    xvals = days
    templates.insert(0,base_template)

    fig = plt.figure(figsize=(10, 5))
    fig, ax = plot_2d_array(
        file_exists_all, xvals, yvals, zrange=(0.0, 2.0), fig_in=fig, cmap="Greens"
    )
    ax.set_yticks(yvals)
    ax.set_yticklabels(templates)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax.tick_params(axis="x", labelrotation=90)
    plt.subplots_adjust(left=0.3, bottom=0.25)

    return fig, ax


if __name__ == "__main__":

    import os

    if os.name == "nt":
        output_root = Path("L:/access/amsr2_out_test3")
    elif os.name == "posix":
        output_root = Path("/mnt/ops1p-ren/l/access/amsr2_out_test3")
    else:
        raise ValueError

    for year in range(2020, 2022):
        start_date = datetime.date(year, 1, 1)
        end_date = datetime.date(year, 12, 31)
        footprint_size = 30

        base_template = "amsr2_resamp_tbs"

        template_list = [
            "amsr2_atm_par_era5",
            "amsr2_land_frac_modis",
            "amsr2_ocean_emiss_era5",
            "amsr2_rainfall_rate",
            "amsr2_skt_era5",
            "amsr2_tclw_era5",
            "amsr2_tcwv_era5",
            "amsr2_u10n_era5",
            "amsr2_v10n_era5",
        ]

        num_days = end_date - start_date
        num_files_tested = len(template_list)+1
        file_exist_all = np.zeros((num_files_tested, num_days.days + 1), dtype=np.int32)

        date_to_do = start_date
        date_index = 0
        while date_to_do <= end_date:
            exists, missing_files = inventory_access_output_files(
                output_root, date_to_do, footprint_size, base_template,template_list
            )
            file_exist_all[:, date_index] = exists

            if np.sum(exists) < num_files_tested:
                print(missing_files)
                print()

            date_to_do += datetime.timedelta(days=1)
            date_index += 1

        fig, ax = plot_summary(start_date, end_date, file_exist_all, base_template,template_list)
        summary_path = output_root / "summaries"
        png_file = summary_path / f"dataset_summary_{year:04d}.png"
        os.makedirs(png_file.parent, exist_ok=True)
        fig.savefig(png_file)

