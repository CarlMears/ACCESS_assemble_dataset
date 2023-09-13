import os
from pathlib import Path

# move tb_orbit files into subdirectory

for orbit_set in range(0, 12):
    num_moved = 0
    start_orbit = 1 + orbit_set * 5000
    end_orbit = start_orbit + 4999

    source_dir = Path(f"L:/access/amsr2_tb_orbits/r{start_orbit:05d}_{end_orbit:05d}/")
    dest_dir = Path(
        f"L:/access/amsr2_tb_orbits/SP_30/r{start_orbit:05d}_{end_orbit:05d}/"
    )
    os.makedirs(dest_dir, exist_ok=True)
    files = os.listdir(source_dir)
    files_to_move = [
        file
        for file in files
        if ((".polar_grid_tb.south." in file) and ("030km" in file))
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
        
        os.rename(source_file, dest_file)
    print(' ')
    print()
