"""Run the RTM and add the atmospheric terms to the resampled TB file.

ERA5 data is downloaded if missing.
"""

import argparse
from contextlib import suppress
import git
from datetime import date
from pathlib import Path

import numpy as np

from netCDF4 import Dataset
from rss_lock.locked_dataset import LockedDataset
import os
import datetime

from access_io.access_output import get_access_output_filename_daily_folder
from access_io.access_output import set_all_attrs
from access_io.access_attr_define import common_global_attributes_access
from access_io.access_attr_define import atm_pars_era5_attributes_access

from util.access_interpolators import time_interpolate_synoptic_maps_ACCESS
from util.file_times import need_to_process
from satellite_definitions.amsr2 import REF_FREQ_mapping

from era5.resample_ERA5 import ResampleERA5

# Reference frequencies (in GHz) to use
# REF_FREQ = np.array([6.9, 7.3, 10.7, 18.7, 23.8, 37.0], np.float32)
# REF_FREQ_mapping = np.array([1,1,2,3,4,5],np.int32)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
# REF_EIA = np.array([53.0, 53.0, 53.0, 53.0, 53.0, 53.0], np.float32)


class OkToSkipDay(Exception):
    pass


class DailyRtm:
    """RTM results for the entire day."""

    def __init__(self, date_to_load: datetime.datetime, data_root: Path):

        filename = (
            f"era5_tbs_{date_to_load.year}-"
            + f"{date_to_load.month:02d}-"
            + f"{date_to_load.day:02d}.nc"
        )
        path_to_data = data_root / filename
        with Dataset(path_to_data, "r") as f:
            tb_down = f["tb_down"][:, :, :, :]
            tb_up = f["tb_up"][:, :, :, :]
            trans = f["tran"][:, :, :, :]

        dims = tb_down.shape
        num_lats = dims[1]
        num_lons = dims[2]
        num_hours = dims[0]
        num_freqs = dims[3]

        shape_rtm_4d = (num_lats, num_lons, num_hours + 1, num_freqs)

        self.downwelling_tb = np.full(shape_rtm_4d, np.nan, dtype=np.float32)
        self.upwelling_tb = np.full(shape_rtm_4d, np.nan, dtype=np.float32)
        self.transmissivity = np.full(shape_rtm_4d, np.nan, dtype=np.float32)

        for hour in range(24):
            for freq in range(num_freqs - 1):
                self.downwelling_tb[:, :, hour, freq] = tb_down[hour, :, :, freq]
                self.upwelling_tb[:, :, hour, freq] = tb_up[hour, :, :, freq]
                self.transmissivity[:, :, hour, freq] = trans[hour, :, :, freq]

        date_to_load_plus_one = date_to_load + datetime.timedelta(days=1)
        filename = (
            f"era5_tbs_{date_to_load_plus_one.year}-"
            + f"{date_to_load_plus_one.month:02d}-"
            + f"{date_to_load_plus_one.day:02d}.nc"
        )
        path_to_data_plus_one = data_root / filename

        with Dataset(path_to_data_plus_one, "r") as f:
            tb_down = f["tb_down"][:, :, :, :]
            tb_up = f["tb_up"][:, :, :, :]
            trans = f["tran"][:, :, :, :]

        for hour in [0]:
            for freq in range(num_freqs - 1):
                self.downwelling_tb[:, :, hour + 24, freq] = tb_down[hour, :, :, freq]
                self.upwelling_tb[:, :, hour + 24, freq] = tb_up[hour, :, :, freq]
                self.transmissivity[:, :, hour + 24, freq] = trans[hour, :, :, freq]

        self.time_in_day = np.arange(0, 25) * 3600.0


def write_atmosphere_to_daily_ACCESS(
    current_day,
    satellite: str,
    target_size: int,
    region: str,
    dataroot: Path,
    outputroot: Path,
    temproot: Path,
    version: str,
    script_name: str,
    commit: str,
    verbose: bool = False,
    overwrite: bool = False,
    update: bool = False,
) -> None:
    if satellite.lower() == "amsr2":
        from satellite_definitions.amsr2 import REF_FREQ


        # Do some logic about the region and file names

    if region == "global":
        grid_type = "equirectangular"
        base_filename = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            target_size,
            dataroot,
            "resamp_tbs",
            grid_type=grid_type,
        )

        atm_filename = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            target_size,
            outputroot,
            "atm_par_era5_temp",
            grid_type=grid_type,
        )
        atm_filename_final = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            target_size,
            outputroot,
            "atm_par_era5",
            grid_type=grid_type,
        )
        grid_type = "equirectangular"
        pole = None
    elif region in ["north", "south"]:
        pole = region
        grid_type = "ease2"
        base_filename = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            target_size,
            dataroot,
            "resamp_tbs",
            grid_type=grid_type,
            pole=pole,
        )
        atm_filename = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            target_size,
            outputroot,
            "atm_par_era5_temp",
            grid_type=grid_type,
            pole=pole,
        )
        atm_filename_final = get_access_output_filename_daily_folder(
            current_day,
            satellite,
            target_size,
            outputroot,
            "atm_par_era5",
            grid_type=grid_type,
            pole=region,
        )
    else:
        raise ValueError(f"region {region} not valid")

    if need_to_process(
        date=current_day,
        satellite=satellite,
        target_size=target_size,
        dataroot=dataroot,
        outputroot=outputroot,
        var="atm_par_era5",
        overwrite=overwrite,
        update=update,
        grid_type=grid_type,
        pole=pole
    ):

        if atm_filename_final.is_file():
            if overwrite or update:
                atm_filename_final.unlink()
            else:
                print(f"File {atm_filename_final} exists, skipping")
                return
        with suppress(FileNotFoundError):
            atm_filename.unlink()

        if verbose:
            print(f"Opening data for {satellite} on {current_day} in {dataroot}")

        try:
            with LockedDataset(base_filename, "r",lock_stale_time = 0.1) as root_grp:

                if verbose:
                    print(
                        f"Reading ERA5 computed RTM data {satellite} "
                        + f"on {current_day} in {temproot}"
                    )

                # num_lats = root_grp["latitude"].shape[0]
                # num_lons = root_grp["longitude"].shape[0]
                # num_hours = root_grp["hours"].shape[0]
                # num_chan = len(REF_FREQ)

                rtm_data = DailyRtm(current_day, temproot)

                resample_required = True
                if (target_size == 30) and (grid_type == "equirectangular"):
                    resample_required = False
                resampler = None

                # WORKING HERE!
                if resample_required:
                    if resampler is None:
                        resampler = ResampleERA5(target_size=target_size, region=region)

                    # time_begin = datetime.datetime.now()
                    # var = resampler.resample_fortran(var)
                    # time = datetime.datetime.now()-time_begin
                    # print(time)time = dateti

                glb_attrs = common_global_attributes_access(
                    current_day,
                    satellite,
                    target_size,
                    version=version,
                    dtype=np.float32,
                )

                atm_attrs = atm_pars_era5_attributes_access(
                    satellite,
                    target_size=target_size,
                    version=version,
                    dtype=np.float32,
                )

                glb_attrs.update(atm_attrs["global"])
                glb_attrs["corresponding_resampled_tb_file"] = base_filename.name
                glb_attrs["script_name"] = script_name
                glb_attrs["commit"] = commit

                with LockedDataset(atm_filename, mode="w",lock_stale_time = 0.1) as trg:
                    for name, dim in root_grp.dimensions.items():
                        if grid_type == 'equirectangular':
                            if name in ["hours", "latitude", "longitude"]:
                                trg.createDimension(
                                    name, len(dim) if not dim.isunlimited() else None
                                )
                        elif grid_type == 'ease2':
                            if name in ["hours", "x", "y"]:
                                trg.createDimension(
                                    name, len(dim) if not dim.isunlimited() else None
                                )
                        else:
                            raise ValueError(f'Grid Type {grid_type} is not valid')

                    trg.createDimension("freq", len(REF_FREQ))

                    # Set global attributes
                    set_all_attrs(trg, glb_attrs)
                    if grid_type == 'equirectangular':
                        for var_name in ["time", "hours", "longitude", "latitude", "freq"]:
                            # Create the time and dimension variables in the output file
                            var_in = root_grp[var_name]
                            trg.createVariable(
                                var_name, var_in.dtype, var_in.dimensions, zlib=True
                            )

                            # Copy the attributes
                            trg.variables[var_name].setncatts(
                                {a: var_in.getncattr(a) for a in var_in.ncattrs()}
                            )
                            trg[var_name][:] = var_in[:]
                        dimensions_out = ("latitude", "longitude", "hours", "freq")
                    elif grid_type == 'ease2':
                        for var_name in ["time", "hours", "x", "y", "freq"]:
                            # Create the time and dimension variables in the output file
                            var_in = root_grp[var_name]
                            trg.createVariable(
                                var_name, var_in.dtype, var_in.dimensions, zlib=True
                            )

                            # Copy the attributes
                            trg.variables[var_name].setncatts(
                                {a: var_in.getncattr(a) for a in var_in.ncattrs()}
                            )
                            trg[var_name][:] = var_in[:]

                        for var_name in ["longitude", "latitude"]:
                            # Create the time and dimension variables in the output file
                            var_in = root_grp[var_name]
                            trg.createVariable(
                                var_name, var_in.dtype, var_in.dimensions, zlib=True
                            )

                            # Copy the attributes
                            trg.variables[var_name].setncatts(
                                {a: var_in.getncattr(a) for a in var_in.ncattrs()}
                            )
                            trg[var_name][:,:] = var_in[:,:]
                        dimensions_out = ("y", "x", "hours", "freq")
                    else:
                        raise ValueError(f'Grid Type {grid_type} is not valid')

                    

                    for varname, long_name, units in [
                        ("transmissivity", "atmospheric transmissivity", None),
                        ("upwelling_tb", "upwelling brightness temperature", "kelvin"),
                        (
                            "downwelling_tb",
                            "downwelling brightness temperature",
                            "kelvin",
                        ),
                    ]:
                        var_attrs = atm_attrs[varname]
                        print(f"starting writing {varname}")
                        if varname == "transmissivity":
                            least_significant_digit = 3
                        else:
                            least_significant_digit = 2
                        trg.createVariable(
                            varname,
                            np.float32,
                            dimensions_out,
                            fill_value=var_attrs["_FillValue"],
                            zlib=True,
                            least_significant_digit=least_significant_digit,
                        )

                        set_all_attrs(trg[varname], var_attrs)
                        # for key in var_attrs.keys():
                        #     if key != "_FillValue":
                        #         value = var_attrs[key]
                        #         set_or_create_attr(trg[varname], key, value)

                        for freq_index, freq in enumerate(REF_FREQ):
                            var = getattr(rtm_data, varname)[
                                :, :, :, REF_FREQ_mapping[freq_index]
                            ]
                            var = np.moveaxis(var, -1, 0)

                            if resample_required:

                                print(f'Resampling {varname} to polar map for freq = {freq}')
                                time_begin = datetime.datetime.now()
                                var = resampler.resample_fortran(var)
                                time = datetime.datetime.now()-time_begin
                                
                            var_times = rtm_data.time_in_day

                            print('Interpolating in time',end="")
                            for hour_index in range(len(root_grp["hours"][:])):
                                time_map = root_grp["time"][:, :, hour_index]
                                time_map = (
                                    time_map
                                    - (
                                        current_day - datetime.date(1900, 1, 1)
                                    ).total_seconds()
                                )
                                var_at_time_map = time_interpolate_synoptic_maps_ACCESS(
                                    var, var_times, time_map
                                )
                                trg[varname][
                                    :, :, hour_index, freq_index
                                ] = var_at_time_map
                                print(".", end="")
                            print()
                        print()
                        print(f"finished writing {varname}")
        except FileNotFoundError:
            raise OkToSkipDay

        with suppress(FileNotFoundError):
            atm_filename_final.unlink()

        os.rename(atm_filename, atm_filename_final)
    else:
        print(f"No Processing needed for atmosperic parameters on {current_day}")


if __name__ == "__main__":
    cds_help = (
        "For downloading ERA5 data from CDS, the UID and API key "
        "must be set as arguments or in the 'CDS_UID' and 'CDS_API_KEY` "
        "environment variables"
    )
    parser = argparse.ArgumentParser(
        description=(
            "Compute and append atmospheric RTM terms to ACCESS output file. "
            "ERA5 data is downloaded if required."
        ),
        epilog=cds_help,
    )
    parser.add_argument(
        "--access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument(
        "--output_root", type=Path, help="Root directory to ACCESS project"
    )

    parser.add_argument(
        "--temp_root", type=Path, help="Root directory to ACCESS project"
    )

    parser.add_argument(
        "--start_date", type=date.fromisoformat, help="Day to process, as YYYY-MM-DD"
    )
    parser.add_argument(
        "--end_date", type=date.fromisoformat, help="Day to process, as YYYY-MM-DD"
    )
    parser.add_argument("--sensor", choices=["amsr2"], help="Microwave sensor to use")
    parser.add_argument(
        "--target_size", choices=["30", "70"], help="target footprint size in km"
    )
    parser.add_argument("--version")
    parser.add_argument(
        "--region", help="region to process", choices=["global", "north", "south"]
    )
    parser.add_argument(
        "--overwrite", help="force overwrite if file exists", action="store_true"
    )
    parser.add_argument(
        "--update", help="force overwrite file older than base", action="store_true"
    )

    args = parser.parse_args()

    access_root: Path = args.access_root
    temp_root: Path = args.temp_root
    output_root: Path = args.output_root
    target_size: int = int(args.target_size)
    region: str = args.region
    version: str = args.version
    rtm_dir = temp_root / "rtm" / "tbs"
    START_DAY = args.start_date
    END_DAY = args.end_date
    satellite = args.sensor
    version = args.version
    overwrite = args.overwrite
    update = args.update

    script_name = parser.prog
    repo = git.Repo(search_parent_directories=True)
    commit = repo.head.object.hexsha

    print("Adding ERA5 drrived RTM parameters")
    print(f"ACCESS Root: {access_root}")
    print(f"Output Root: {output_root}")
    print(f"RTM Root: {rtm_dir}")
    print(f"Date Range:  {START_DAY} - {END_DAY}")
    print(f"Satellite:   {satellite}")
    print(f"Target Size: {target_size}")
    print(f"Region: {region}")
    print(f"Version:     {version}")
    if overwrite:
        print("Overwriting old files")
    else:
        if update:
            print("Updating older files by overwriting")
    print()
    print(f"Script Name: {script_name}")
    print(f"git Commit: {commit}")

    date_to_do = args.start_date
    while date_to_do <= args.end_date:
        try:
            write_atmosphere_to_daily_ACCESS(
                date_to_do,
                satellite,
                target_size,
                region,
                access_root,
                output_root,
                rtm_dir,
                version=version,
                verbose=True,
                overwrite=overwrite,
                update=update,
                script_name=script_name,
                commit=commit,
            )
        except OkToSkipDay:
            print(f"Problem finding file...skipping {date_to_do}")
        date_to_do = date_to_do + datetime.timedelta(days=1)
