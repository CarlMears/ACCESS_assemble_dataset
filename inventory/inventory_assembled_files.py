from pathlib import Path
import numpy as np
import datetime
import matplotlib.pyplot as plt
from rss_plotting.plot_2d_array import plot_2d_array

from access_io.access_output import get_access_output_filename_daily_folder


def inventory_access_output_files(
    *,
    output_root: Path,
    date_to_do: datetime.date,
    footprint_size: int,
    satellite="AMSR2",
    ksat="15",
    grid_type="equirectangular",
    pole="",
    base_template: str,
    templates: list[str],
) -> tuple:
    debug = True
    list_of_missing_files = []

    if satellite.upper() == "SMAP":
        exists = np.zeros((len(templates) + 1, 2),np.int32)
    else:
        exists = np.zeros(len(templates) + 1, np.int32)

    if satellite.upper() == "SMAP":
        for look in [0,1]:
            base_filename = get_access_output_filename_daily_folder(
                date_to_do,
                satellite,
                footprint_size,
                output_root,
                base_template,
                grid_type=grid_type,
                pole=pole,
                ksat=ksat,
                look=look
            )

            if not base_filename.is_file():
                list_of_missing_files.append(base_filename.name)
                exists[0,look] = 0
                if debug:
                    print(f"Base File Missing: {base_filename}")
            else:
                exists[0,look] = 2
                base_file_time = base_filename.stat().st_mtime

            if exists[0,look] > 0:
                for file_index, file_template in enumerate(templates):
                    var_filename = get_access_output_filename_daily_folder(
                        date_to_do,
                        satellite,
                        footprint_size,
                        output_root,
                        file_template,
                        grid_type=grid_type,
                        pole=pole,
                        ksat=ksat,
                        look=look
                    )

                    if not var_filename.is_file():
                        list_of_missing_files.append(var_filename.name)
                        exists[file_index + 1,look] = 0
                        if debug:
                            print(f"Missing: {var_filename}")
                    else:
                        file_time = var_filename.stat().st_mtime
                        if file_time > base_file_time:
                            exists[file_index + 1,look] = 2
                        else:
                            exists[file_index + 1,look] = 1
    else:
        base_filename = get_access_output_filename_daily_folder(
                date_to_do,
                satellite,
                footprint_size,
                output_root,
                base_template,
                grid_type=grid_type,
                pole=pole,
                ksat=ksat,
            )

        if not base_filename.is_file():
            list_of_missing_files.append(base_filename.name)
            exists[0] = 0
            if debug:
                print(f"Base File Missing: {base_filename}")
        else:
            exists[0] = 2
            base_file_time = base_filename.stat().st_mtime

        if exists[0] > 0:
            for file_index, file_template in enumerate(templates):
                var_filename = get_access_output_filename_daily_folder(
                    date_to_do,
                    satellite,
                    footprint_size,
                    output_root,
                    file_template,
                    grid_type=grid_type,
                    pole=pole,
                    ksat=ksat,
                )

                if not var_filename.is_file():
                    list_of_missing_files.append(var_filename.name)
                    exists[file_index + 1] = 0
                    if debug:
                        print(f"Missing: {var_filename}")
                else:
                    file_time = var_filename.stat().st_mtime
                    if file_time > base_file_time:
                        exists[file_index + 1] = 2
                    else:
                        exists[file_index + 1] = 1

    return exists, list_of_missing_files


def plot_summary(
    start_date: datetime.date,
    end_date: datetime.date,
    file_exists_all: np.ndarray,
    base_template: Path,
    templates: list[Path],
):
    import matplotlib.dates as mdates

    yvals = np.arange(len(templates) + 1)
    xvals = np.arange((end_date - start_date).days + 1)
    days = mdates.drange(
        start_date, end_date + datetime.timedelta(days=1), datetime.timedelta(days=1)
    )
    xvals = days
    templates.insert(0, base_template)

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

def plot_summary_two_looks(
    start_date: datetime.date,
    end_date: datetime.date,
    file_exists_all: np.ndarray,
    base_template: Path,
    templates: list[Path],
        ):
    import matplotlib.dates as mdates

    yvals = np.arange(2*(len(templates) + 1))
    xvals = np.arange((end_date - start_date).days + 1)
    days = mdates.drange(
        start_date, end_date + datetime.timedelta(days=1), datetime.timedelta(days=1)
    )
    xvals = days
    templates.insert(0, base_template)

    fig = plt.figure(figsize=(10, 5))

    file_exists_all_to_plot = np.zeros((2*len(templates), file_exists_all.shape[2]),dtype=np.int32)

    for i in range(len(templates)):
        file_exists_all_to_plot[2*i,:] = file_exists_all[i,0,:]
        file_exists_all_to_plot[2*i+1,:] = file_exists_all[i,1,:]

    fig, ax = plot_2d_array(
        file_exists_all_to_plot, xvals, yvals, zrange=(0.0, 2.0), fig_in=fig, cmap="Greens"
    )
    ytick_vals = np.arange(len(templates))*2+1
    ax.set_yticks(ytick_vals)
    ax.set_yticklabels(templates)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax.tick_params(axis="x", labelrotation=90)
    plt.subplots_adjust(left=0.3, bottom=0.25)

    return fig, ax



if __name__ == "__main__":
    import os

    if os.name == "nt":
        output_root = Path("L:/access")
    elif os.name == "posix":
        output_root = Path("/mnt/l/access/")
    else:
        raise ValueError

    satellite = "smap"
    ksat="15"
    footprint_size = 70
    region = "global"

    if region == "global":
        grid_type = "equirectangular"
        output_root = output_root / f"{satellite}_out_GL_{footprint_size}"
        pole = None
    elif region in ["north", "south"]:
        pole = region
        grid_type = "ease2"
        if pole == "north":
            output_root = output_root / f"{satellite}_out_NP_{footprint_size}"
        elif pole == "south":
            output_root = output_root / f"{satellite}_out_SP_{footprint_size}"
        else:
            raise ValueError(f"pole {pole} not valid")
    else:
        raise ValueError(f"region {region} not valid")

    for year in range(2015, 2023):
        start_date = datetime.date(year, 1, 1)
        end_date = datetime.date(year, 12, 31)

        base_template = "resamp_tbs"

        template_list = [
            "atm_par_era5",
            "land_frac_modis",
            "ocean_emiss_era5",
            "rainfall_rate",
            "skt_era5",
            "tclw_era5",
            "tcwv_era5",
            "u10n_era5",
            "v10n_era5",
        ]

        num_days = end_date - start_date
        num_files_tested = len(template_list) + 1

        if satellite.lower()=='smap':
            file_exist_all=np.zeros((num_files_tested,2,num_days.days+1),dtype=np.int32)
        else:
            file_exist_all = np.zeros((num_files_tested, num_days.days + 1), dtype=np.int32)

        date_to_do = start_date
        date_index = 0
        while date_to_do <= end_date:
            exists, missing_files = inventory_access_output_files(
                output_root=output_root,
                date_to_do=date_to_do,
                satellite=satellite,
                ksat=ksat,
                footprint_size=footprint_size,
                grid_type=grid_type,
                pole=pole,
                base_template=base_template,
                templates=template_list,
            )

            if np.sum(exists) < num_files_tested:
                print(missing_files)
                print()

            if satellite=='smap':
                file_exist_all[:, :, date_index] = exists
            else:
                file_exist_all[:, date_index] = exists

            date_to_do += datetime.timedelta(days=1)
            date_index += 1

        if satellite.lower()=='smap':
            fig,axis = plot_summary_two_looks(
                start_date, end_date, file_exist_all, base_template, template_list
            )
        else:
            fig, ax = plot_summary(
                start_date, end_date, file_exist_all, base_template, template_list
            )
        if satellite=='ssmi':
            summary_path = output_root / f'f{int(ksat):02d}'/ "summaries"
        else:
            summary_path = output_root / "summaries"
        png_file = summary_path / f"dataset_summary_{year:04d}.png"
        os.makedirs(png_file.parent, exist_ok=True)
        print(png_file)
        fig.savefig(png_file)
