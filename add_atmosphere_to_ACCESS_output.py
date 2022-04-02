"""Run the RTM and add the atmospheric terms to the resampled TB file.

ERA5 data is downloaded if missing.
"""

import argparse
import os
from datetime import date, timedelta
from pathlib import Path
from typing import Dict, Sequence

import numpy as np
from access_atmosphere import rtm
from access_atmosphere.download import Era5Downloader
from access_atmosphere.era5 import Era5DailyData, read_era5_data
from netCDF4 import Dataset
from numpy.typing import NDArray

from access_io.access_output import get_access_output_filename

# TODO: once Python 3.9 is the minimum supported version, remove the above
# imports from typing and switch to: "from collections.abc import Sequence" and
# use dict[] instead of Dict[]


# Reference frequencies (in GHz) to use
REF_FREQ: NDArray[np.float32] = np.array(
    [6.8, 10.7, 18.7, 23.8, 37.0, 89.0], np.float32
)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
REF_EIA: NDArray[np.float32] = np.array(
    [53.0, 53.0, 53.0, 53.0, 53.0, 53.0], np.float32
)


class DailyRtm:
    """RTM results for the entire day."""

    def __init__(
        self,
        lat: NDArray[np.float32],
        lon: NDArray[np.float32],
        time: NDArray[np.int32],
    ) -> None:
        """Allocate space for the results based on the size of the inputs."""
        self.lat = lat
        self.lon = lon
        self.time = time

        num_lat = len(self.lat)
        num_lon = len(self.lon)
        num_time = len(self.time)
        num_freq = len(REF_FREQ)

        # The shapes for the arrays all match the output netCDF variables
        shape_3d = (num_lat, num_lon, num_time)
        shape_4d = (num_lat, num_lon, num_time, num_freq)
        self.col_water_vapor: NDArray[np.float32] = np.full(
            shape_3d, np.nan, np.float32
        )
        self.col_cloud_liquid: NDArray[np.float32] = np.full(
            shape_3d, np.nan, np.float32
        )
        self.transmissivity: NDArray[np.float32] = np.full(shape_4d, np.nan, np.float32)
        self.tb_up: NDArray[np.float32] = np.full(shape_4d, np.nan, np.float32)
        self.tb_down: NDArray[np.float32] = np.full(shape_4d, np.nan, np.float32)

    def compute_atmosphere(
        self,
        era5_data: Era5DailyData,
        times: Sequence[int],
        valid_data: NDArray[np.bool_],
    ) -> None:
        """Compute the atmospheric RTM for the input ERA5 data.

        The `valid_data` is a boolean ndarray with dimensions (lats, lons, times).

        Note that this mutates `self`.
        """
        # Extract the ERA5 data corresponding to the valid data mask. Note that
        # the ERA5 profile data is dimensioned as (time, lats, lons, levels),
        # but the valid_mask has dimensions (lats, lons, time). So the
        # valid_mask is reshaped to match the ERA5 data.
        valid_in = np.moveaxis(valid_data[:, :, times], 2, 0)
        atmo_results = rtm.compute(
            era5_data.levels,
            era5_data.temperature[valid_in],
            era5_data.height[valid_in],
            era5_data.specific_humidity[valid_in],
            era5_data.liquid_content[valid_in],
            era5_data.surface_temperature[valid_in],
            era5_data.surface_height[valid_in],
            era5_data.surface_dewpoint[valid_in],
            era5_data.surface_pressure[valid_in],
            REF_EIA,
            REF_FREQ,
        )

        # Store the vectorized RTM outputs into the corresponding valid values
        # in the full output arrays
        valid_out = np.full_like(valid_data, False)
        valid_out[:, :, times] = valid_data[:, :, times]

        self.col_water_vapor[valid_out] = era5_data.columnar_water_vapor[valid_in]
        self.col_cloud_liquid[valid_out] = era5_data.columnar_cloud_liquid[valid_in]
        self.transmissivity[valid_out] = atmo_results.tran
        self.tb_up[valid_out] = atmo_results.tb_up
        self.tb_down[valid_out] = atmo_results.tb_down


class DailyAccessData:
    """One day of resampled TB data for a satellite."""

    def __init__(self, current_day: date, satellite: str, dataroot: Path) -> None:
        """Initialize based on the date, satellite, and root path."""
        self.tb_path = get_access_output_filename(current_day, satellite, dataroot)
        with Dataset(self.tb_path, "r") as f:
            # Each of these are 1d coordinate variables
            lat = f["latitude"][:]
            lon = f["longitude"][:]
            hours = f["hours"][:]
            # Dimensioned as (lat, lon, hours)
            time = f["time"][...]

        # The mask for the time data will be used as the mask for valid data. In
        # other words, only non-masked time values will be used to compute the
        # RTM.
        self.valid_data = ~np.ma.getmaskarray(time)
        self.time: NDArray[np.int64] = np.ma.getdata(time)

        self.reference_day = current_day

        # These are coordinate arrays, so none of them should have masked data
        if any(np.ma.is_masked(a) for a in (lat, lon, hours)):
            raise Exception("Unexpected masked data")
        self.lat: NDArray[np.float32] = np.ma.getdata(lat)
        self.lon: NDArray[np.float32] = np.ma.getdata(lon)
        self.hours: NDArray[np.int32] = np.ma.getdata(hours)

    @property
    def num_points(self) -> int:
        """Return the number of valid points."""
        return np.count_nonzero(self.valid_data)

    def append_results(self, data: DailyRtm) -> None:
        """Update the daily TB file in-place with the RTM results for the day."""
        with Dataset(self.tb_path, "a") as f:
            _ensure_rtm_vars(f)
            f.variables["columnar_water_vapor"][...] = data.col_water_vapor
            f.variables["columnar_cloud_liquid"][...] = data.col_cloud_liquid
            f.variables["transmissivity"][...] = data.transmissivity
            f.variables["upwelling_tb"][...] = data.tb_up
            f.variables["downwelling_tb"][...] = data.tb_down

    def interpolate_era5(
        self, hour: int, data_prev: Era5DailyData, data_next: Era5DailyData
    ) -> Era5DailyData:
        """Time-interpolate the ERA5 data to the measurement times."""
        # Check compatibility and inputs
        if hour not in range(24):
            raise Exception("Hour must be between 0 and 23")
        if not all(
            np.array_equal(getattr(data_prev, name), getattr(data_next, name))
            for name in ["levels", "lats", "lons"]
        ):
            raise Exception("Mismatched ERA5 inputs")
        if data_prev.time.shape != (1,) or data_next.time.shape != (1,):
            raise Exception("ERA5 data should only contain one hour each")

        valid_data = self.valid_data[:, :, hour]

        # Convert measurement times to the ERA5 epoch: from seconds since
        # 1900-01-01 to hours since 1900-01-01.
        time = self.time[:, :, hour].astype(np.float32) / (60 * 60)
        time[~valid_data] = np.nan

        # Interpolate all the content. It's actually a little awkward to use the
        # numpy or scipy interpolators for this, it's easier to just do it manually.
        fractional_time = (time - data_prev.time) / (data_next.time - data_prev.time)
        fractional_time[~valid_data] = 0

        # The "fractional time" should vary from 0 to 1, but if any of the time
        # data are outside the ERA5 min/max values for this hour bin, it will go
        # outside that range. If this happens, it's probably a problem in the
        # file data to fix. At any rate, check for those out-of-bounds
        # conditions and warn about it, but then clamp the data so it stays
        # between 0 and 1 so that interpolation will work.
        too_low = fractional_time < 0
        too_high = fractional_time > 1
        if too_low.any() or too_high.any():
            out_of_range = too_low.sum() + too_high.sum()
            print(
                f"WARNING: clamping {out_of_range} fractional times "
                "to between 0 and 1"
            )
            np.clip(fractional_time, 0, 1, out=fractional_time)

        fractional_time_3d = fractional_time[..., np.newaxis]

        def lerp(
            t: NDArray[np.float32], x0: NDArray[np.float32], x1: NDArray[np.float32]
        ) -> NDArray[np.float32]:
            """Linear interpolation."""
            return (1 - t) * x0 + t * x1

        interpolated_3d_data = {
            variable: lerp(
                fractional_time_3d,
                getattr(data_prev, variable),
                getattr(data_next, variable),
            )
            for variable in (
                "temperature",
                "specific_humidity",
                "height",
                "liquid_content",
            )
        }
        interpolated_2d_data = {
            variable: lerp(
                fractional_time,
                getattr(data_prev, variable),
                getattr(data_next, variable),
            )
            for variable in (
                "surface_pressure",
                "surface_temperature",
                "surface_dewpoint",
                "surface_height",
                "columnar_water_vapor",
                "columnar_cloud_liquid",
            )
        }

        for data in interpolated_3d_data.values():
            data[0, ~valid_data] = np.nan
        for data in interpolated_2d_data.values():
            data[0, ~valid_data] = np.nan

        return Era5DailyData(
            data_prev.levels,
            data_prev.lats,
            data_prev.lons,
            time.mean(),
            interpolated_3d_data["temperature"],
            interpolated_3d_data["specific_humidity"],
            interpolated_3d_data["height"],
            interpolated_3d_data["liquid_content"],
            interpolated_2d_data["surface_pressure"],
            interpolated_2d_data["surface_temperature"],
            interpolated_2d_data["surface_dewpoint"],
            interpolated_2d_data["surface_height"],
            interpolated_2d_data["columnar_water_vapor"],
            interpolated_2d_data["columnar_cloud_liquid"],
        )


def _ensure_rtm_vars(f: Dataset) -> None:
    """Ensure RTM output variables are defined in the output file."""
    if "frequency" not in f.dimensions:
        f.createDimension("frequency", len(REF_FREQ))
        v = f.createVariable("frequency", np.float32, ("frequency",))
        v[:] = REF_FREQ
        v.setncatts(
            {
                "standard_name": "sensor_band_central_radiation_frequency",
                "long_name": "frequency",
                "units": "GHz",
            }
        )

    for varname, std_name, long_name in [
        (
            "columnar_water_vapor",
            "atmosphere_mass_content_of_water_vapor",
            "columnar water vapor",
        ),
        (
            "columnar_cloud_liquid",
            "atmosphere_mass_content_of_cloud_liquid_water",
            "columnar liquid cloud content",
        ),
    ]:
        if varname not in f.variables:
            v = f.createVariable(
                varname,
                np.float32,
                ("latitude", "longitude", "hours"),
                zlib=True,
            )
            v.setncatts(
                {
                    "standard_name": std_name,
                    "long_name": long_name,
                    "units": "kg m-2",
                }
            )

    for varname, long_name, units in [
        ("transmissivity", "atmospheric transmissivity", None),
        ("upwelling_tb", "upwelling brightness temperature", "kelvin"),
        ("downwelling_tb", "downwelling brightness temperature", "kelvin"),
    ]:
        if varname not in f.variables:
            v = f.createVariable(
                varname,
                np.float32,
                ("latitude", "longitude", "hours", "frequency"),
                zlib=True,
            )
            v.setncattr("long_name", long_name)
            if units is not None:
                v.setncattr("units", units)
            v.setncattr(
                "reference_incidence_angle",
                f"{REF_EIA[0]} degrees",
            )


def append_atmosphere_to_daily_ACCESS(
    current_day: date,
    satellite: str,
    dataroot: Path,
    downloader: Era5Downloader,
    verbose: bool = False,
) -> None:
    """Append the atmospheric terms to daily ACCESS resampled-TB dataset.

    The file is read in order to determine the valid points where the RTM should
    be called. The resulting values are written by appending the new variables
    to the file.

    ERA5 data, both on the surface and as profiles, is required. Since the ERA5
    data is time-interpolated to the ACCESS data, two days of ERA5 data are
    required for each ACCESS data file.

    If the ACCESS resampled TB file doesn't exist, a `FileNotFound` error will
    be raised.

    The configured `downloader` is used to download the ERA5 data, if needed.
    """
    if verbose:
        print(f"Reading data for {satellite} on {current_day} in {dataroot}")
    daily_data = DailyAccessData(current_day, satellite, dataroot)
    if verbose:
        valid_points = daily_data.num_points
        total_points = np.prod(daily_data.valid_data.shape)
        print(
            f"Number of valid points: {valid_points} "
            f"({valid_points / total_points:0.2%})"
        )

    next_day = current_day + timedelta(days=1)
    downloader.download_day(current_day, verbose)
    downloader.download_day(next_day, verbose)

    rtm_data = DailyRtm(daily_data.lat, daily_data.lon, daily_data.hours)

    era5_path = downloader.out_dir
    era5_surface_current_day = era5_path / f"era5_surface_{current_day.isoformat()}.nc"
    era5_levels_current_day = era5_path / f"era5_levels_{current_day.isoformat()}.nc"
    era5_surface_next_day = era5_path / f"era5_surface_{next_day.isoformat()}.nc"
    era5_levels_next_day = era5_path / f"era5_levels_{next_day.isoformat()}.nc"

    # To reduce the memory usage, process the data one hour at a time. However,
    # the ERA5 data is interpolated in time to the satellite measurement times,
    # so two hours from the ERA5 data are read for each output hour. In order to
    # avoid re-reading ERA5 data for the same hour, they're cached in a dict.
    NUM_HOURS = 24
    era5_cache: Dict[int, Era5DailyData] = {}
    for hour in range(NUM_HOURS):
        if verbose:
            print(f"Processing hour {hour+1}/{NUM_HOURS}")

        # Read the ERA5 data for the two hours that bracket the measurement
        # data, noting that for the last hour, data from the next day is
        # required.
        #
        # The era5_cache dict is used to look up the previous-hour data, which
        # will always succeed except for the first hour of this loop when the
        # cache is empty. For the next-hour data, the cache is not checked since
        # the next-hour data will never be present in the cache. The values are
        # popped from the dict since a cached entry is only reused once and this
        # reduces the memory usage by not storing data from multiple hours ago.
        try:
            era5_data_prev = era5_cache.pop(hour)
        except KeyError:
            era5_data_prev = read_era5_data(
                era5_surface_current_day,
                era5_levels_current_day,
                [hour],
                verbose,
            )
        if hour + 1 < NUM_HOURS:
            era5_data_next = read_era5_data(
                era5_surface_current_day,
                era5_levels_current_day,
                [hour + 1],
                verbose,
            )
        else:
            era5_data_next = read_era5_data(
                era5_surface_next_day,
                era5_levels_next_day,
                [0],
                verbose,
            )
        era5_cache[hour + 1] = era5_data_next

        if verbose:
            print("Interpolating ERA5 data to measurement times")
        era5_data = daily_data.interpolate_era5(hour, era5_data_prev, era5_data_next)

        if verbose:
            print("Computing RTM")
        rtm_data.compute_atmosphere(era5_data, (hour,), daily_data.valid_data)

    if verbose:
        print(f"Appending results to: {daily_data.tb_path}")
    daily_data.append_results(rtm_data)


if __name__ == "__main__":
    cds_help = (
        "For downloading ERA5 data from CDS, the UID and API key "
        "must be set as arguments or in the 'CDS_UID' and 'CDS_API_KEY` "
        "environment variables"
    )
    parser = argparse.ArgumentParser(
        description=(
            "Compute and append atmospheric RTM terms to ACCESS output file. "
            "ERA5 data is downloaded if required."
        ),
        epilog=cds_help,
    )
    parser.add_argument(
        "access_root", type=Path, help="Root directory to ACCESS project"
    )
    parser.add_argument(
        "date", type=date.fromisoformat, help="Day to process, as YYYY-MM-DD"
    )
    parser.add_argument("sensor", choices=["amsr2"], help="Microwave sensor to use")
    parser.add_argument(
        "--user",
        metavar="UID",
        help="CDS UID (overrides CDS_UID environment variable)",
    )
    parser.add_argument(
        "--key",
        metavar="API_KEY",
        help="CDS API key (overrides CDS_API_KEY environment variable)",
    )
    args = parser.parse_args()

    if args.user is not None:
        cds_uid = args.user
    else:
        try:
            cds_uid = os.environ["CDS_UID"]
        except KeyError:
            parser.error(cds_help)

    if args.key is not None:
        cds_api_key = args.key
    else:
        try:
            cds_api_key = os.environ["CDS_API_KEY"]
        except KeyError:
            parser.error(cds_help)

    access_root: Path = args.access_root
    era5_dir = access_root / "era5"
    downloader = Era5Downloader(cds_uid, cds_api_key, era5_dir)

    append_atmosphere_to_daily_ACCESS(
        args.date, args.sensor, access_root, downloader, verbose=True
    )
