import concurrent.futures
import datetime
import signal
from pathlib import Path
from threading import Lock

import numpy as np
import pyproj as proj
import xarray as xr

from resampling_utils.AMSR2_Antenna_Gain import target_gain
from resampling_utils.resample_imerg_polar import ResampleIMERG
NUM_LATS = 721
NUM_LONS = 1440
NUM_LATS_EASE2 = 720
NUM_LONS_EASE2 = 720
NUM_HOURS = 24
hdf5_access = Lock()


def init_worker() -> None:
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
):
    """
    Reading latitude, longitude, and precipitation from downloaded IMERG files.
    """
    minute_string = f"{minutes_of_day:04}"
    date_string = date.strftime("%Y%m%d")

    if minute_string == "1440":
        date += datetime.timedelta(days=1)
        date_string = date.strftime("%Y%m%d")
        minute_string = "0000"
        print(minute_string, date_string)

        filename = Path(next(target_path.glob(f"*.{date_string}*.{minute_string}*")))

    else:
        print(minute_string, date_string)
        filename = Path(next(target_path.glob(f"*.{date_string}*.{minute_string}*")))

    last_modified = filename.stat().st_mtime
    # multiple threads do not play nice opening HDF5 files
    with hdf5_access:
        hr1 = xr.open_dataset(filename, group="Grid")
        rain = hr1["precipitationCal"].values
        lat = hr1["lat"].values
        lon = hr1["lon"].values
        hr1.close()

    return lat, lon, rain, last_modified


def initialize_wgs84():
    """
    Initializing WGS84 projection that we wish to project onto.
    Using this function to initialize this once saves us run time.
    """
    g = proj.Geod(ellps="WGS84")  # assuming you're using WGS84 geographic

    return g


def get_distance(latitude0, longitude0, ilat_submask, ilon_submask, g):
    """
    Calculating distance between a center point and surrounding lats/lons
    using WGS84 geographic.
    """

    lon0_array = np.repeat(longitude0, np.size(ilon_submask))
    lat0_array = np.repeat(latitude0, np.size(ilat_submask))

    latv, lonv = np.meshgrid(ilat_submask, ilon_submask, indexing="ij")
    lat0v, lon0v = np.meshgrid(lat0_array, lon0_array, indexing="ij")
    fw_az, bk_az, dist_grid = g.inv(lon0v, lat0v, lonv, latv)

    return dist_grid


def resample_to_quarter(
    map_rain, lat_rain, lon_rain, mask, footprint_diameter_km, window=0.5
):
    """
    Inputs a map of rain rates at a given time along with
    latitude and longitude, outputs resampled data on a 0.25x0.25 grid.
    """
    resampled_map = np.full((NUM_LATS, NUM_LONS), np.nan)
    g = initialize_wgs84()

    for i in range(NUM_LONS):
        for j in range(NUM_LATS):
            if not mask[j, i]:
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

            # For percentage NaN checks - I am debating whether or not we should
            # perform these checks at all. Currently removing windows around
            # 0.25 degree grid box where over 50% of IMERG values are NaN.
            where_nan = np.where(np.isnan(rain))
            n_elements_rain = rain.size
            n_elements_nan = len(where_nan[0])

            if (n_elements_nan / n_elements_rain) > 0.5:
                # if over 50% of window is NaN, do not include
                resampled_map[j, i] = np.nan
                continue

            if ~np.any(rain != 0.0):
                # no need to do a weighted average calculation if all rain rates
                # in window are zero
                resampled_map[j, i] = 0.0
            else:
                dist = get_distance(
                    lat_quarter, lon_quarter, lat_rain[lat_ok], lon_rain[lon_ok], g
                )

                dist_km = dist / 1000.0
                # a check based on Thomas' IMERG work. Removing any negative RR values
                good_rain = np.where(rain > -0.01)

                gains = target_gain(dist_km, diameter_in_km=footprint_diameter_km)
                weighted_rain = np.average(rain[good_rain], weights=gains[good_rain])
                resampled_map[j, i] = weighted_rain

    return resampled_map

def resample_to_ease(map_rain, mask, resampler):
    """
    Inputs a map of rain rates at a given time along with
    latitude and longitude, outputs resampled data on an EASE2 grid.
    """
    resampled_map = np.full((NUM_LATS_EASE2, NUM_LONS_EASE2), np.nan)

    var_in = np.squeeze(map_rain)
    var_in = np.transpose(var_in)
    var_in = np.roll(var_in, shift=1800, axis=1)

    var_out = resampler.resample_fortran(var_in)
    # var_out = np.transpose(var_out) # will need to be removed when Fortran is fiex
    resampled_map = var_out
    resampled_map[~mask] = np.nan


    return resampled_map

def resample_hour(
    hour, times, time_intervals, date, footprint_diameter_km, region, resampler, target_path
):
    sat_time = times[:, :, hour]

    hour_beg = np.all(  # from beginning of hour to 30 minutes past the hour
        [
            (sat_time >= time_intervals[hour]),
            (sat_time < (time_intervals[hour + 1] - 1800)),
        ],
        axis=(0),
    )

    hour_end = np.all(  # from 30 minutes past the hour to the following hour
        [
            (sat_time >= (time_intervals[hour + 1] - 1800)),
            (sat_time < time_intervals[hour + 1]),
        ],
        axis=(0),
    )

    # Opening all files needed to perform resampling for given hour (need 2 IMERG files)
    minutes_of_day_beg = int(time_intervals[hour] / 60)
    minutes_of_day_end = int((time_intervals[hour + 1] / 60) - 30)

    lat_beg, lon_beg, rain_beg, last_modified_beg = read_imerg_half_hourly(
        minutes_of_day=minutes_of_day_beg, date=date, target_path=target_path
    )

    lat_end, lon_end, rain_end, last_modified_beg = read_imerg_half_hourly(
        minutes_of_day=minutes_of_day_end, date=date, target_path=target_path
    )


    if region == 'global':
        res = np.full((2, NUM_LATS, NUM_LONS), np.nan)

        mp = resample_to_quarter(
            map_rain=rain_beg,
            lat_rain=lat_beg,
            lon_rain=lon_beg,
            mask=hour_beg,
            footprint_diameter_km=footprint_diameter_km,
            window=0.5,
        )
        res[0, :, :] = mp
        mp = resample_to_quarter(
            map_rain=rain_end,
            lat_rain=lat_end,
            lon_rain=lon_end,
            mask=hour_end,
            footprint_diameter_km=footprint_diameter_km,
            window=0.5,
        )
        res[1, :, :] = mp

        hour_map = np.nanmean((res), axis=0)
    elif region in ['north', 'south']:
        res = np.full((2, NUM_LATS_EASE2, NUM_LONS_EASE2), np.nan)
        mp = resample_to_ease(
            map_rain=rain_beg,
            mask=hour_beg,
            resampler=resampler,
        )

        res[0,:,:] = mp

        mp = resample_to_ease(
            map_rain=rain_end,
            mask=hour_end,
            resampler=resampler,
        )

        res[1,:,:] = mp

        hour_map = np.nanmean((res), axis=0)

    return (hour_map, hour, last_modified_beg)


def resample_imerg_day(
    times, 
    time_intervals,
    date,
    footprint_diameter_km,
    region, 
    target_path=Path("."),resampler=False
):
    if region == 'global':
        total_hour = np.full((NUM_LATS, NUM_LONS, NUM_HOURS), np.nan)
        resampler = False
    elif region in ['north', 'south']:
        total_hour = np.full((NUM_LATS_EASE2, NUM_LONS_EASE2, NUM_HOURS), np.nan)
        # Only initialize the resampler if we need to
        if not isinstance(resampler, ResampleIMERG):
            raise ValueError("Need to initialize resampler for EASE2 resampling")
    else:
        print(f"{region} not recognized")

    # Using process pool for resampling global

    if region == 'global':
        with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
            results = {
                executor.submit(
                    resample_hour,
                    hour,
                    times,
                    time_intervals,
                    date,
                    footprint_diameter_km,
                    target_path,
                ): hour
                for hour in range(0, NUM_HOURS)
            }
            for future in concurrent.futures.as_completed(results):
                try:
                    (map, idx, modtime) = future.result()
                    total_hour[:, :, idx] = map
                except KeyboardInterrupt:
                    return
                except Exception as e:
                    print(f"Error in run: {e}")
    elif region in ['north', 'south']:  # no parallel processing if EASE2 resampling
        for hour in range(0, NUM_HOURS):
            map, idx, modtime = resample_hour(hour,times,time_intervals,date,footprint_diameter_km, region, resampler, target_path)
            total_hour[:, :, idx] = map
    return total_hour, modtime
