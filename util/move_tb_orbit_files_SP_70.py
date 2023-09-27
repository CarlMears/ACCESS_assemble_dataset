import os
from pathlib import Path
from time import sleep

# move tb_orbit files into subdirectory

for orbit_set in range(0, 12):
    num_moved = 0
    start_orbit = 1 + orbit_set * 5000
    end_orbit = start_orbit + 4999

    if os.name=="posix":
        source_dir = Path(f"/mnt/l/access/amsr2_tb_orbits/r{start_orbit:05d}_{end_orbit:05d}/")
        dest_dir = Path(
            f"/mnt/l/access/amsr2_tb_orbits/SP_70/r{start_orbit:05d}_{end_orbit:05d}/"
        )
    elif os.name == 'nt':
        source_dir = Path(f"L:/access/amsr2_tb_orbits/r{start_orbit:05d}_{end_orbit:05d}/")
        dest_dir = Path(
            f"L:/access/amsr2_tb_orbits/SP_70/r{start_orbit:05d}_{end_orbit:05d}/"
        )
    else:
        raise ValueError(f"OS {os.name} not supported")
    
    os.makedirs(dest_dir, exist_ok=True)
    files = os.listdir(source_dir)
    files_to_move = [
        file
        for file in files
        if ((".polar_grid_tb.south." in file) and ("070km" in file))
    ]

    for f in files_to_move:
        source_file = source_dir / f
        dest_file = dest_dir / f
        if num_moved < 10:
            print(f"Moving {source_file} to {dest_file}")
        elif num_moved == 10:
            print("...", end="")
        num_moved += 1
        if num_moved % 100 == 0:
            print('.', end="")
        
        num_tries = 0
        success = False
        while num_tries < 10:
            try:
                os.rename(source_file, dest_file)
                success = True
                break
            except OSError:
                sleep(1.0)
                num_tries += 1
        if ~success:
            f"Could not move {source_file} to {dest_file} in {num_tries} tries"
             
    print(' ')
    print()
