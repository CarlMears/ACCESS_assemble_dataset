import datetime

import numpy as np
from numpy.typing import NDArray


def calendar_dates_from_datetime64(dt) -> dict[str, NDArray[np.uint32]]:
    """
    Convert array of datetime64 to a calendar array of year, month, day, hour,
    minute, seconds, microsecond with these quantites indexed on the last axis.

    Parameters
    ----------
    dt : datetime64 array (...)
        numpy.ndarray of datetimes of arbitrary shape

    Returns
    -------
    cal : dictionary of arrays
            'Y' = year, 'M' = month, 'D' =day, 'h' = hour,
            'm' = minute, 's' = second, 'us' = microsecond
            keys are the same of the datetime64 unit names except that
            I left out 'ms' and 'ns'
    """
    # decompose calendar floors
    Y, M, D, h, m, s = [dt.astype(f"M8[{x}]") for x in "YMDhms"]
    return {
        "Y": (Y + 1970).astype("u4"),  # Gregorian Year
        "M": ((M - Y) + 1).astype("u4"),  # month
        "D": ((D - M) + 1).astype("u4"),  # dat
        "h": ((dt - D).astype("m8[h]")).astype("u4"),  # hour
        "m": ((dt - h).astype("m8[m]")).astype("u4"),  # minute
        "s": ((dt - m).astype("m8[s]")).astype("u4"),  # second
        "us": ((dt - s).astype("m8[us]")).astype("u4"),  # microsecond
    }


def convert_to_sec_in_day(
    ob_time, date: datetime.date, ref_year: int = 2000
) -> NDArray[np.float32]:
    """Converts at time, in seconds since Jan 1, ref_year to
    seconds in day.  Reference year defaults to 2000.  Not
    sure how leap seconds affect this, so DO NOT USE for
    geolocation"""
    date_jan1_2000 = np.datetime64(f"{ref_year:04d}-01-01T00:00:00")
    start_of_day = np.datetime64(f"{date:%Y-%m-%d}T00:00:00")

    bad = ob_time < 100.0
    ob_time = ob_time.astype(np.int64)
    dt_obs = ob_time.astype("timedelta64[s]")
    obs_datetime = dt_obs + date_jan1_2000
    obtime_in_day = (obs_datetime - start_of_day).astype(np.float32)
    obtime_in_day[bad] = np.nan

    return obtime_in_day


def convert_to_np_datetime64(ob_time, ref_year: int = 2000) -> NDArray[np.timedelta64]:
    """converts at time, in seconds since Jan 1, ref_year to
    seconds in day.  Reference year defaults to 2000, the value
    used for the AMSR2 file.  I am not sure what this does about leap seconds,
    so DO NOT USE for geolocation until we check it out"""
    date_jan1_2000 = np.datetime64(f"{ref_year:04d}-01-01T00:00:00")
    # bad = ob_time < 0.00001  # ad data is often stored as zero

    ob_time = ob_time.astype(np.int64)
    dt_obs = ob_time.astype("timedelta64[s]")
    obs_datetime64 = dt_obs + date_jan1_2000

    return obs_datetime64
