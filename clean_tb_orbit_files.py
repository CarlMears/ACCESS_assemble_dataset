from pathlib import Path
import numpy as np


def range_str(iorbit):
    irange = int(np.floor((iorbit - 1) / 5000))
    lower_orbit = 1 + irange * 5000
    upper_orbit = lower_orbit + 4999

    return f"r{lower_orbit:05d}_{upper_orbit:05d}"


if __name__ == "__main__":
    import os

    if os.name == "nt":
        tb_root = Path("L:/access/amsr2_tb_orbits")
    elif os.name == "posix":
        tb_root = Path("/mnt/l/access/amsr2_tb_orbits")
    else:
        raise ValueError

    start_orbit = 35001
    end_orbit = 55000
    channel_list = [5, 6, 7, 8, 9, 10, 11, 12]

    footprint_size = 30

    template_list = ["grid", "dist"]

    for iorbit in range(start_orbit, end_orbit + 1):
        rng_str = range_str(iorbit)
        if (iorbit % 100) == 0:
            print(iorbit)
        for file_template in template_list:
            for channel in channel_list:
                filename = (
                    tb_root
                    / rng_str
                    / f"r{iorbit:05d}.{file_template}.{channel:02d}.nc"
                )
                try:
                    isFile = filename.is_file()
                except Exception:
                    isFile = True

                if isFile:
                    if ((iorbit % 100) == 0) and (channel == 5):
                        print(filename)
                    max_tries = 10
                    num_tries = 0
                    while num_tries < max_tries:
                        try:
                            filename.unlink()
                            break
                        except Exception:
                            pass
                            num_tries += 1
                            print(num_tries)
