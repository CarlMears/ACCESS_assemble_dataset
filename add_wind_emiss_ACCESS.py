import argparse
import datetime
from contextlib import suppress
from pathlib import Path

import git
import numpy as np
from geomod10 import wind_emiss  # python wrapper for geomod10b and geomod10c
from rss_lock.locked_dataset import LockedDataset

from access_io.access_output import write_ocean_emiss_to_daily_ACCESS
from access_io.access_output_polar import write_ocean_emiss_to_daily_ACCESS_polar
from access_io.access_attr_define import common_global_attributes_access
from access_io.access_attr_define import (
    anc_var_attributes_access,
)  # old version

from util.file_times import need_to_process

from geomod10 import wind_emiss  # python wrapper for geomod10b and geomod10c
from access_io.access_attr_define import (  # old version
    anc_var_attributes_access,
    common_global_attributes_access,
)
from access_io.access_output import (
    get_access_output_filename_daily_folder,
    write_ocean_emiss_to_daily_ACCESS,
)
from satellite_definitions.amsr2 import REF_EIA, REF_FREQ

# these are for debugging only
# from rss_plotting.global_map import plot_global_map
# import matplotlib.pyplot as plt


class OkToSkipDay(Exception):
    pass


def calc_emissivity_maps(
    *, date, wind_source, sst_source, target_size, grid_type, pole
):

    print(f"{date}")

    # Get u wind info from data repository
    anc_name = f"u10n_{wind_source}"
    var_filename_final = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, access_root, anc_name, grid_type, pole
    )
    # ds = xr.open_dataset(var_filename_final)
    with LockedDataset(var_filename_final, "r", 60) as root_grp:
        # read in existing u10n
        try:
            u10n = root_grp.variables[anc_name][:, :, :]
        except KeyError:
            raise KeyError(
                f"Variable: {anc_name} in not dataset, can not find emissivity"
            )

    # Get v wind info from data repository
    anc_name = f"v10n_{wind_source}"
    var_filename_final = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, access_root, anc_name, grid_type, pole
    )

    with LockedDataset(var_filename_final, "r", 60) as root_grp:
        # read in existing v10n
        try:
            v10n = root_grp.variables[anc_name][:, :, :]
        except KeyError:
            raise KeyError(
                f"Variable: {anc_name} in not dataset, can not find emissivity"
            )

    w10n = np.sqrt(np.square(u10n) + np.square(v10n)).filled(fill_value=np.nan)

    anc_name = f"skt_{sst_source}"
    var_filename_final = get_access_output_filename_daily_folder(
        date, satellite.lower(), target_size, access_root, anc_name, grid_type, pole
    )
    with LockedDataset(var_filename_final, "r", 60) as root_grp:
        try:
            sst = root_grp.variables[anc_name][:, :, :]
        except KeyError:
            raise KeyError(
                f"Variable: {anc_name} in not dataset, can not find emissivity"
            )

    wind_ok = np.isfinite(w10n)
    
    w10n_ok = w10n[wind_ok]

    sst = sst - 273.15
    sst[sst < -3.0] = 3.0
    sst[sst > 34.0] = 34.0
    sst_ok = sst[wind_ok]

    phir = -999.0
    phir = np.full_like(sst_ok, -999.0)

    num_freq = len(REF_FREQ)
    num_pol = 2
    sz = w10n.shape
    num_lat = sz[0]
    num_lon = sz[1]
    num_time = sz[2]
    shape_5d = (num_lat, num_lon, num_time, num_freq, num_pol)
    shape_3d = shape_5d[0:3]

    emiss_maps = np.full(shape_5d, np.nan, dtype=np.float32)

    for ifreq, freq in enumerate(REF_FREQ):
        print(f"Processing {freq:.3f}")
        tht = REF_EIA[ifreq]
        if np.sum(wind_ok) > 0:
            emiss_this_freq = wind_emiss(freq, tht, sst_ok, w10n_ok, phir)

            for ipol in [0, 1]:
                maps_temp = np.full(shape_3d, np.nan, dtype=np.float32)
                maps_temp[wind_ok] = emiss_this_freq[:, ipol]
                emiss_maps[:, :, :, ifreq, ipol] = maps_temp

    return emiss_maps


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            (
                """Calculate and Add Wind emissivity based on
                 wind variables already present in the ACCESS
                 dataset."""
            )
        ),
    )
    parser.add_argument(
        "--access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument(
        "--output_root",
        type=Path,
        help="directory to output data",
    )

    parser.add_argument("--wind", choices=["era5"], help="Source of wind info")
    parser.add_argument("--sst", choices=["era5"], help="Source of wind info")

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
    parser.add_argument("--target_size", choices=["30", "70"], help="Target size in km")
    parser.add_argument(
        "--region", help="region to process", choices=["global", "north", "south"]
    )
    parser.add_argument("--version", help="String describing version to be created")
    parser.add_argument(
        "--verbose", help="enable more verbose screen output", action="store_true"
    )
    parser.add_argument(
        "--overwrite", help="overwrite output file no matter what", action="store_true"
    )
    parser.add_argument(
        "--update",
        help="overwrite output file if older than base file",
        action="store_true",
    )

    
    args = parser.parse_args()
    access_root: Path = args.access_root
    output_root: Path = args.output_root
    wind_source = args.wind
    sst_source = args.sst
    version = args.version
    satellite = args.sensor
    target_size = int(args.target_size)
    region = args.region
    overwrite = args.overwrite
    update = args.update

    script_name = parser.prog
    repo = git.Repo(search_parent_directories=True)
    commit = repo.head.object.hexsha

    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor.upper()
    date = START_DAY
while date <= END_DAY:
    emiss_filename_final = get_access_output_filename_daily_folder(
        date, satellite, target_size, output_root, "ocean_emiss_era5"
    )

    if region in ["north", "south"]:
        grid_type = "ease2"
        pole = region
    elif region == "global":
        grid_type = "equirectangular"
        pole = None
    else:
        raise ValueError(f"Region: {region} not valid")

    while date <= END_DAY:

        emiss_filename_final = get_access_output_filename_daily_folder(
            date,
            satellite,
            target_size,
            output_root,
            "ocean_emiss_era5",
            grid_type=grid_type,
            pole=pole,
        )

        base_filename = get_access_output_filename_daily_folder(
            date,
            satellite,
            target_size,
            access_root,
            "resamp_tbs",
            grid_type=grid_type,
            pole=pole,
        )
        var = "ocean_emiss_era5"
        if need_to_process(
            date=date,
            satellite=satellite,
            target_size=target_size,
            grid_type=grid_type,
            pole=pole,
            dataroot=access_root,
            outputroot=output_root,
            var=var,
            overwrite=overwrite,
            update=update,
        ):
            with suppress(FileNotFoundError):
                emiss_filename_final.unlink()

            ocean_emiss = calc_emissivity_maps(
                date=date, wind_source="era5", sst_source="era5", target_size=target_size,grid_type=grid_type, pole=pole
            )

            # common global_attributes for the project
            glb_attrs = common_global_attributes_access(
                date, satellite, target_size, version
            )

            # variable-specific attributes
            var = "ocean_emiss"
            var_attrs = anc_var_attributes_access(satellite, var + "_era5", version)

            # add the global part of these to the global_attrs
            glb_attrs.update(var_attrs["global"])
            glb_attrs["corresponding_resampled_tb_file"] = base_filename.name
            glb_attrs["script_name"] = script_name
            glb_attrs["commit"] = commit

            # keep the variable decription parts in var_attrs

            var_attrs = var_attrs["var"]

            if grid_type == 'equirectangular':
                write_ocean_emiss_to_daily_ACCESS(
                    ocean_emiss=ocean_emiss,
                    current_day=date,
                    satellite=satellite,
                    target_size=target_size,
                    glb_attrs=glb_attrs,
                    var_attrs=var_attrs,
                    dataroot=access_root,
                    outputroot=output_root,
                    verbose=True
                )
            elif grid_type == 'ease2':
                write_ocean_emiss_to_daily_ACCESS_polar(
                    ocean_emiss=ocean_emiss,
                    current_day=date,
                    satellite=satellite,
                    target_size=target_size,
                    grid_type=grid_type,
                    pole=pole,
                    global_attrs=glb_attrs,
                    var_attrs=var_attrs,
                    dataroot=access_root,
                    outputroot=output_root,
                    verbose=True,
                )
            else:
                raise ValueError(f'grid_type: {grid_type} is not valid')
        else:
            print(f"No Processing Needed for {var} on {date}")

        date = date + datetime.timedelta(days=1.0)