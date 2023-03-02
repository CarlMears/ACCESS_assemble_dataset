import datetime
from pathlib import Path
from access_io.access_output import get_access_output_filename_daily_folder
from resampled_tbs.read_resampled_orbit import get_resampled_file_name

def get_mtime_multi_try(path_to_file,max_num_trys=10):

    num_so_far = 0
    while num_so_far < max_num_trys:
        try:
            mtime = path_to_file.stat().st_mtime
            return mtime
        except:
            #catch all errors
            num_so_far += 1
    
    raise IOError(f'Problem getting mtime for {path_to_file}')


def need_to_process(*, 
                    date: datetime.date, 
                    satellite: str, 
                    target_size: int, 
                    dataroot: Path, 
                    outputroot: Path,
                    var:str,
                    overwrite: bool,
                    update: bool):



    base_filename = get_access_output_filename_daily_folder(
        date, satellite, target_size, dataroot, "resamp_tbs"
        )

    var_filename = get_access_output_filename_daily_folder(
                date, satellite, target_size, outputroot, var)

    if base_filename.is_file():
        if var_filename.is_file():
            if overwrite:
                return True
            else:
                if update:
                    base_file_time = get_mtime_multi_try(base_filename)
                    var_file_time = get_mtime_multi_try(var_filename)
                    if base_file_time > var_file_time:
                        return True
                    else:
                        return False
                else:
                    return False
        else:
            return True
    else:
        return False

def need_to_process_base_file(*, 
                    date: datetime.date, 
                    satellite: str, 
                    target_size: int, 
                    grid_type: str = 'equirectangular',
                    pole: str,
                    orbits_to_do: list[int],
                    channels_to_do: list[int],
                    dataroot: Path, 
                    outputroot: Path,
                    var:str,
                    overwrite: bool,
                    update: bool):

    if len(orbits_to_do) <= 0:
        return False

    base_filename = get_access_output_filename_daily_folder(
        date, satellite, target_size, dataroot, 
        "resamp_tbs",grid_type=grid_type,pole=pole
    )

    found_tb_orbit_file = False
    tb_orbit_file_newer_than_base_file = False

    if base_filename.is_file():
        base_file_time = get_mtime_multi_try(base_filename)
    else:
        base_file_time = 0.0

    for orbit in orbits_to_do:
        for channel in channels_to_do:
            tb_orbit_file = get_resampled_file_name(
                satellite=satellite,
                channel=channel,
                target_size=target_size,
                orbit=orbit,
                grid_type=grid_type,
                pole=pole
            )
            if tb_orbit_file.is_file():
                found_tb_orbit_file = True
                tb_file_time = get_mtime_multi_try(tb_orbit_file)
                if tb_file_time > base_file_time:
                    tb_orbit_file_newer_than_base_file

    if base_filename.is_file():
        if found_tb_orbit_file:
            if overwrite:
                return True
            else:
                if update:
                    if tb_orbit_file_newer_than_base_file:
                        return True
                    else:
                        return False
                else:
                    return False 
        else:
            return False #no orbit files
    else:
        if found_tb_orbit_file:
            return True
        else:
            return False  