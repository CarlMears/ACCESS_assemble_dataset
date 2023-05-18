import datetime
import os
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from util.numpy_date_utils import convert_to_np_datetime64

if os.name == "nt":
    AMSR2_ORBIT_FILE = Path("j:/amsr2/tables/orbit_times.dat")
elif os.name == "posix":
    AMSR2_ORBIT_FILE = Path("/mnt/ops1p-ren/j/amsr2/tables/orbit_times.dat")


def read_amsr2_orbit_times(
    filename_orb: Path = AMSR2_ORBIT_FILE,
) -> NDArray[Any]:
    times = np.fromfile(filename_orb, dtype=np.float64)
    return convert_to_np_datetime64(times, ref_year=1993)


def find_orbits_in_day(
    *, times_np64: NDArray[Any], date: datetime.date, verbose: bool = False
) -> NDArray[Any]:
    target_day_begin_np64 = np.datetime64(f"{date:%Y-%m-%d}T00:00")
    target_day_end_np64 = target_day_begin_np64 + np.timedelta64(24, "h")

    try:
        orbits = np.where(
            np.all(
                [
                    (times_np64 > target_day_begin_np64),
                    (times_np64 < target_day_end_np64),
                ],
                axis=0,
            )
        )[0]
        orbits = np.insert(orbits, 0, orbits[0] - 1)
    except IndexError:
        orbits = np.empty((0),dtype=np.int32)
        return orbits

    orbits = np.append(orbits, orbits[-1] + 1)

    assert len(orbits) < 19  # This should always be true

    if verbose:
        print(times_np64[orbits])
    orbits += 1  # fix indexing offset in fortran file
    if verbose:
        print(orbits)
    return orbits


if __name__ == "__main__":
    times = read_amsr2_orbit_times()
    date = datetime.date(2012, 7, 11)
    orbits = find_orbits_in_day(times_np64=times, date=date)
    print(orbits)
