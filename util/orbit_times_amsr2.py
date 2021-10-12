def read_amsr2_orbit_times(filename_orb="j:/amsr2/tables/orbit_times.dat"):
    import numpy as np
    from util.numpy_date_utils import convert_to_np_datetime64

    times = np.fromfile(filename_orb, dtype=np.float64)
    times_np64 = convert_to_np_datetime64(times, ref_year=1993)

    return times_np64


def find_orbits_in_day(*, times_np64, year, month, day, verbose=False):
    import numpy as np

    orbits = []

    target_day_begin_np64 = np.datetime64(f"{year:04d}-{month:02d}-{day:02d}T00:00")
    target_day_end_np64 = target_day_begin_np64 + np.timedelta64(24, "h")

    orbits = np.where(
        np.all(
            [(times_np64 > target_day_begin_np64), (times_np64 < target_day_end_np64)],
            axis=0,
        )
    )[0]
    orbits = np.insert(orbits, 0, orbits[0] - 1)
    orbits = np.append(orbits, orbits[-1] + 1)

    assert len(orbits < 19)  # This should always be true

    if verbose:
        print(times_np64[orbits])
    orbits = orbits + 1  # fix indexing offset in fortran file
    if verbose:
        print(orbits)
    return orbits


if __name__ == "__main__":

    times = read_amsr2_orbit_times()
    year = 2012
    month = 7
    day = 11
    orbits = find_orbits_in_day(times_np64=times, year=year, month=month, day=day)
    print(orbits)
    print
