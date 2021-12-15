import datetime
import os
import numpy as np
import xarray as xr
import pyproj as proj
from resampling_utils.AMSR2_Antenna_Gain import *
from pathlib import Path
import multiprocessing
import signal


NUM_LATS = 721
NUM_LONS = 1440


def init_worker():
    """
    This is a function which will allow the user to stop all workers
    with a Ctrl-C event.  Otherwise the code gets caught in an odd loop in
    which it does not terminate.  This appears to be a Python bug:
    https://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def read_imerg_half_hourly(
    *,
    minutes_of_day: int,
    date: datetime.date,
    target_path: Path,
) -> None:
    """
    Reading latitude, longitude, and precipitation from downloaded IMERG files.
    """
    minute_string = str(minutes_of_day).zfill(4)
    date_string = date.strftime("%Y%m%d")
    dir = target_path

    files = os.listdir(dir)

    if minute_string == "1440":
        date += datetime.timedelta(days=1)
        date_string = date.strftime("%Y%m%d")
        minute_string = "0000"
        print(minute_string, date_string)

        file = [i for i in files if (f".{minute_string}." in i) and (date_string in i)]
        filename = Path(dir / file[0])

    else:
        print(minute_string, date_string)
        file = [i for i in files if (f".{minute_string}." in i) and (date_string in i)]
        filename = Path(dir / file[0])
    print(filename)

    hr1 = xr.open_dataset(filename, group="Grid")
    rain = hr1["precipitationCal"].values
    lat = hr1["lat"].values
    lon = hr1["lon"].values
    hr1.close()

    return lat, lon, rain


def initialize_wgs84():
    """
    Initializing WGS84 projection that we wish to project onto.
    Using this function to initialize this once saves us run time.
    """
    crs_wgs = proj.crs.CRS("epsg:4326")  # assuming you're using WGS84 geographic

    return crs_wgs


def get_distance(latitude0, longitude0, ilat_submask, ilon_submask, crs_wgs):
    """
    Getting a projection for distance between a central grid point
    and an array of IMERG lats/lons. This is more accurate than
    assuming a completely spherical Earth.
    """

    # set up a custom projection, centered at latitude0,longitude0
    cust = proj.crs.CRS(
        f"+proj=aeqd +lat_0={latitude0} +lon_0={longitude0} +datum=WGS84 +units=m"
    )

    # set up a complete grid of lats/lons to be transformed.  the ilat_submask and ilon submask are 1-d arrays that define the grid.
    latv, lonv = np.meshgrid(ilat_submask, ilon_submask, indexing="ij")

    # apply the local projection to obtain the mesh points in meters from the center specified by  latitude0,longitude0
    transformer = proj.Transformer.from_crs(crs_wgs, cust, always_xy=True)
    xv, yv = transformer.transform(lonv, latv)

    return xv, yv


def resample_to_quarter(map_rain, lat_rain, lon_rain, mask, window=0.5):
    """
    Inputs a map of rain rates at a given time along with
    latitude and longitude, outputs resampled data on a 0.25x0.25 grid.
    """
    resampled_map = np.full((NUM_LATS, NUM_LONS), np.nan)
    crs_wgs = initialize_wgs84()

    for i in range(0, 1440):
        for j in range(0, 721):
            if mask[j, i] == False:
                continue
            lon_quarter = (i * 0.25) - 180.0
            lat_quarter = (j * 0.25) - 90.0

            lat_ok = np.where(
                (lat_rain >= (lat_quarter - window))
                & (lat_rain < (lat_quarter + window))
            )

            # Handling cases of longitude wrapping
            if (lon_quarter - window) < -180.0:
                wrap_lon = lon_quarter - window
                wrap_lon += 360.0
                lon_ok_reg = np.where(
                    (lon_rain >= (lon_quarter - window))
                    & (lon_rain < (lon_quarter + window))
                )
                lon_ok_wrap = np.where(lon_rain > wrap_lon)
                lon_ok = np.concatenate([lon_ok_reg[0], lon_ok_wrap[0]])
            elif (lon_quarter + window) > 180.0:
                wrap_lon = lon_quarter + window
                wrap_lon -= 360.0
                lon_ok_reg = np.where(
                    (lon_rain >= (lon_quarter - window))
                    & (lon_rain < (lon_quarter + window))
                )
                lon_ok_wrap = np.where(lon_rain < wrap_lon)
                lon_ok = np.concatenate([lon_ok_reg[0], lon_ok_wrap[0]])
            else:
                lon_ok = np.where(
                    (lon_rain >= (lon_quarter - window))
                    & (lon_rain < (lon_quarter + window))
                )

            x, y = np.meshgrid(lon_ok, lat_ok)
            rain = map_rain[0, x, y]

            # For percentage NaN checks - I am debating whether or not we should perform these checks at all
            # Currently removing windows around 0.25 degree grid box where over 50% of IMERG values
            # are NaN
            where_nan = np.where(np.isnan(rain))
            sz = rain.shape
            n_elements_rain = sz[0] * sz[1]
            n_elements_nan = len(where_nan[0])

            if (
                n_elements_nan / n_elements_rain
            ) > 0.5:  # if over 50% of window is NaN, do not include
                resampled_map[j, i] = np.nan
                continue

            if np.all(
                np.isnan(rain)
            ):  # if no IMERG data in window, set 0.25 deg pixel to NaN
                resampled_map[j, i] = np.nan
            elif (
                np.all(rain) == 0.0
            ):  # no need to do a weighted average calculation if all rain rates in window are zero
                resampled_map[j, i] = 0.0
            else:
                xv, yv = get_distance(
                    lat_quarter,
                    lon_quarter,
                    lat_rain[lat_ok],
                    lon_rain[lon_ok],
                    crs_wgs,
                )

                dist_km = ((xv / 1000.0) ** 2.0 + (yv / 1000.0) ** 2.0) ** 0.5
                good_rain = np.where(
                    rain > -0.01
                )  # a check based on Thomas' IMERG work. Removing any negative RR values

                gains = target_gain(dist_km, diameter_in_km=30.0)
                weighted_rain = np.nansum(
                    rain[good_rain] * gains[good_rain]
                ) / np.nansum(gains[good_rain])
                resampled_map[j, i] = weighted_rain

    return resampled_map


def resample_imerg_day(times, time_intervals, date, target_path=""):
    total_hour = np.full((NUM_LATS, NUM_LONS, 24), np.nan)

    for idx in range(0, 24):
        sat_time = times[:, :, idx]

        hour_beg = np.all(  # from beginning of hour to 15 minutes past the hour
            [
                (sat_time >= time_intervals[idx]),
                (sat_time < (time_intervals[idx + 1] - 2400)),
            ],
            axis=(0),
        )

        hour_mid = np.all(  # from 15 minutes past the hour to 45 minutes past the hour
            [
                (sat_time >= (time_intervals[idx + 1] - 2400)),
                (sat_time < (time_intervals[idx + 1] - 1200)),
            ],
            axis=(0),
        )

        hour_end = np.all(  # from 45 minutes past the hour to the following hour
            [
                (sat_time >= (time_intervals[idx + 1] - 1200)),
                (sat_time < time_intervals[idx + 1]),
            ],
            axis=(0),
        )

        # Opening all files needed to perform resampling for given hour (need 3 IMERG files)
        minutes_of_day_beg = int(time_intervals[idx] / 60)
        minutes_of_day_mid = int((time_intervals[idx + 1] - 1800) / 60)
        minutes_of_day_end = int(time_intervals[idx + 1] / 60)

        lat_beg, lon_beg, rain_beg = read_imerg_half_hourly(
            minutes_of_day=minutes_of_day_beg, date=date, target_path=target_path
        )
        lat_mid, lon_mid, rain_mid = read_imerg_half_hourly(
            minutes_of_day=minutes_of_day_mid, date=date, target_path=target_path
        )
        lat_end, lon_end, rain_end = read_imerg_half_hourly(
            minutes_of_day=minutes_of_day_end, date=date, target_path=target_path
        )

        for_parallel = [
            [rain_beg, lat_beg, lon_beg, hour_beg],
            [rain_mid, lat_mid, lon_mid, hour_mid],
            [rain_end, lat_end, lon_end, hour_end],
        ]

        # Processing 3 IMERG files for one hour in parallel
        p = multiprocessing.Pool(
            5, init_worker
        )  # 5 workers at the moment. Can be changed

        try:
            print("Starting jobs")
            res = p.starmap(resample_to_quarter, for_parallel)
            print("Waiting for results")
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            p.terminate()
        else:
            print("Normal termination")
            p.close()
        p.join()  # this waits for all worker processes to terminate.

        total_hour[:, :, idx] = np.nanmean(
            (res), axis=0
        )  # this should be ok because there is, in theory, no overlap between the three resampled IMERG maps (only data on satellite orbit tracks)

    return total_hour
