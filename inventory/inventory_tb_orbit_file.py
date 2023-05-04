from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from rss_plotting.plot_2d_array import plot_2d_array


def is_file_multiple_try(path_to_test: Path, max_tries: int = 10):

    num_tries = 0
    is_file_multiple = False
    while num_tries < max_tries:
        try:
            is_file_multiple = path_to_test.is_file()
            break
        except OSError:
            print(f"{num_tries}: {path_to_test}")
            pass
    return is_file_multiple


def inventory_access_tb_orbit_files(
    *,
    tb_orbit_root: Path,
    start_orbit: int,
    end_orbit: int,
    channel_list: list,
    footprint_size: int,
    file_template: str,
) -> np.ndarray:

    debug = False
    if os.name == "nt":
        AMSR2_L2B_root = Path("J:/AMSR2/L2B_V08.2/")
    elif os.name == "posix":
        AMSR2_L2B_root = Path("/mnt/ops1p-ren/j/AMSR2/L2B_V08.2")

    num_orbits = end_orbit - start_orbit + 1
    exists = np.zeros((len(channel_list), num_orbits), dtype=np.int32)

    for orbit_num in range(start_orbit, end_orbit + 1):
        if orbit_num % 500 == 0:
            print(orbit_num)
        lower_range = (5000 * int(np.floor((orbit_num - 1) / 5000.0))) + 1
        upper_range = lower_range + 4999
        tb_orbit_root_this_orbit = (
            tb_orbit_root / f"r{lower_range:05d}_{upper_range:05d}"
        )
        l2b_root_this_orbit = AMSR2_L2B_root / f"r{lower_range:05d}_{upper_range:05d}"
        for ich, channel in enumerate(channel_list):
            filename = (
                tb_orbit_root_this_orbit / f"r{orbit_num:05d}.{file_template}."
                f"ch{channel:02d}.{footprint_size:03d}km.nc"
            )
            if debug:
                print(filename)
            isFile = is_file_multiple_try(filename)
            if isFile:
                exists[channel - channel_list[0], orbit_num - start_orbit] = 2
            else:
                l2b_filename = l2b_root_this_orbit / f"r{orbit_num:05d}.dat"
                isFileL2B = is_file_multiple_try(l2b_filename)
                if isFileL2B:
                    exists[channel - channel_list[0], orbit_num - start_orbit] = 0
                    print(f"{filename} is missing but L2B file exists")
                else:
                    exists[channel - channel_list[0], orbit_num - start_orbit] = 1
    return exists


def plot_orbit_summary(
    start_orbit: int, end_orbit: int, file_exists_all: np.ndarray, channel_list: list
):

    channel_names = []
    for channel in channel_list:
        channel_names.append(f"ch{channel:02d}")

    yvals = np.arange(len(channel_list))
    xvals = np.arange(start_orbit, end_orbit + 1)

    fig = plt.figure(figsize=(10, 5))
    fig, ax = plot_2d_array(
        file_exists_all, xvals, yvals, zrange=(0.0, 2.0), fig_in=fig, cmap="Greens"
    )
    ax.set_yticks(yvals)
    ax.set_yticklabels(channel_names)
    plt.subplots_adjust(left=0.3, bottom=0.25)

    return fig, ax


if __name__ == "__main__":

    import os

    if os.name == "nt":
        tb_orbit_root = Path("L:/access/amsr2_tb_orbits")
    elif os.name == "posix":
        tb_orbit_root = Path("/mnt/ops1p-ren/l/access/amsr2_tb_orbits")
    else:
        raise ValueError

    channel_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    channel_name = [
        "6.9V",
        "6.9H",
        "7.3V",
        "7.3H",
        "11V",
        "11H",
        "19V",
        "19H",
        "24V",
        "24H",
        "37V",
        "37H",
    ]

    footprint_size = 30
    region = "north"

    if region == "global":
        template = "grid_tb"
        tb_orbit_root = tb_orbit_root / f"GL_{footprint_size:02d}"
    elif region == "north":
        template = "polar_grid_tb.north"
        tb_orbit_root = tb_orbit_root / f"NP_{footprint_size:02d}"
    elif region == "south":
        template = "polar_grid_tb.south"
        tb_orbit_root = tb_orbit_root / f"SP_{footprint_size:02d}"
    else:
        raise ValueError(f"Region: {region} is not valid")

    for orbit_group in range(0, 12):
        start_orbit = 5000 * orbit_group + 1
        end_orbit = start_orbit + 4999

        exists = inventory_access_tb_orbit_files(
            tb_orbit_root=tb_orbit_root,
            start_orbit=start_orbit,
            end_orbit=end_orbit,
            channel_list=channel_list,
            footprint_size=footprint_size,
            file_template=template,
        )

        print(np.sum(exists))
        fig, ax = plot_orbit_summary(start_orbit, end_orbit, exists, channel_list)
        ax.set_title(
            (
                f"AMSR2 Tb Orbit Inventory {start_orbit:05d}-{end_orbit:05d}"
                f" {region} {footprint_size:03d}km"
            )
        )
        print()
        png_file = (
            tb_orbit_root
            / "inventory"
            / (
                f"tb_orbit_inv_r{start_orbit:05d}_{end_orbit:05d}."
                f"{footprint_size:03d}km.{region}.png"
            )
        )
        os.makedirs(png_file.parent, exist_ok=True)
        fig.savefig(png_file)
        print(f"Saving: {png_file}")
