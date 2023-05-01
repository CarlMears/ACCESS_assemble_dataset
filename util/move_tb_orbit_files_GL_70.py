import os
from pathlib import Path

# move tb_orbit files into subdirectory


for orbit_set in range(0,12):
    start_orbit = 1 + orbit_set*5000
    end_orbit = start_orbit+4999


    source_dir = Path(f'L:/access/amsr2_tb_orbits/r{start_orbit:05d}_{end_orbit:05d}/')
    dest_dir = Path(f'L:/access/amsr2_tb_orbits/GL_70/r{start_orbit:05d}_{end_orbit:05d}/')

    files = os.listdir(source_dir)

    files_to_move = [file for file in files if (('.grid_tb.' in file) and ('070km' in file))]

    for f in files_to_move:
        source_file = source_dir / f
        dest_file = dest_dir / f
        print(f'Moving {source_file} to {dest_file}')
        os.rename(source_file,dest_file)
    print()

