import os
from pathlib import Path

# move tb_orbit files into subdirectory
if os.name == "nt":
    LDRIVE = Path("L:/")
elif os.name == "posix":
    LDRIVE = Path("/mnt/l/")    
else:
    raise ValueError(f'FILE SYSTEM {os.name} not SUPPORTED')

for orbit_set in range(0, 14):
    start_orbit = 1 + orbit_set * 5000
    end_orbit = start_orbit + 4999

    source_dir = Path(f"{LDRIVE}/access/amsr2_tb_orbits/r{start_orbit:05d}_{end_orbit:05d}/")
    dest_dir = Path(
        f"{LDRIVE}/access/amsr2_tb_orbits/GL_30/r{start_orbit:05d}_{end_orbit:05d}/"
    )

    files = os.listdir(source_dir)

    files_to_move = [
        file for file in files if ((".grid_tb." in file) and ("030km" in file))
    ]

    for f in files_to_move:
        source_file = source_dir / f
        dest_file = dest_dir / f
        print(f"Moving {source_file} to {dest_file}")
        os.rename(source_file, dest_file)
    print()

    files_to_move = [
        file for file in files if ((".time." in file) and (".nc" in file))
    ]

    for f in files_to_move:
        source_file = source_dir / f
        dest_file = dest_dir / f
        print(f"Moving {source_file} to {dest_file}")
        os.rename(source_file, dest_file)
    print()
