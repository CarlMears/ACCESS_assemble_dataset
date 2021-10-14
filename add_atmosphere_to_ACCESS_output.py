"""Run the RTM and add the atmospheric terms to the resampled TB file.

This depends on the atmospheric RTM package, access-atmosphere:
http://gitlab.remss.com/access/atmospheric-rtm
"""

import os
from datetime import date
from pathlib import Path
from typing import Sequence, NamedTuple, cast

import numpy as np
from netCDF4 import Dataset

from access_io.access_output import get_access_output_filename
from access_atmosphere.download import Era5Downloader
from access_atmosphere.era5 import Era5DailyData, read_era5_data
from access_atmosphere import rtm

# from resampled_tbs.read_resampled_orbit import read_resampled_tbs

# Reference frequencies (in GHz) to use
REF_FREQ = np.array([6.8, 10.7, 18.7, 23.8, 37.0, 89.0], np.float32)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
REF_EIA = np.array([53.0, 53.0, 53.0, 53.0, 53.0, 53.0], np.float32)


class HourlyRtm(NamedTuple):
    """RTM results for some number of hour inputs."""

    # Note, the type annotations here are sub-optimal until numpy has better
    # typing support. They are all ndarrays.

    # Latitudes dimensioned as (lat, ).
    lat: Sequence[np.float32]

    # Longitudes dimensioned as (lon, ).
    lon: Sequence[np.float32]

    # Hours since midnight, dimensioned as (time, ).
    hours: Sequence[int]

    # Total column water vapor in kg/m^2, dimensioned as (lat, lon, time).
    columnar_water_vapor: Sequence[np.float32]

    # Total column cloud liquid water in kg/m^2, dimensioned as (lat, lon,
    # time).
    columnar_cloud_liquid: Sequence[np.float32]

    # Atmospheric transmissivity, unitless, dimensioned as (lat, lon, time, freq).
    transmissivity: Sequence[np.float32]

    # Atmospheric upwelling brightness temperature in K, dimensioned as
    # (lat, lon, time, freq).
    tb_up: Sequence[np.float32]

    # Atmospheric downwelling brightness temperature in K, dimensioned as
    # (lat, lon, time, freq).
    tb_down: Sequence[np.float32]


class DailyAccessData:
    """One day of resampled TB data for a satellite."""

    def __init__(self, current_day: date, satellite: str, dataroot: Path) -> None:
        """Initialize based on the date, satellite, and root path."""
        self.tb_path = get_access_output_filename(current_day, satellite, dataroot)
        with Dataset(self.tb_path, "r") as f:
            self.lat = f["latitude"][:]
            self.lon = f["longitude"][:]
            self.time = f["hours"][:]
            # Dimensioned as (lat, lon, time, channel)
            tb = f["brightness_temperature"][...]

        # This boolean array is True whereever there is valid data. The channel
        # dimension is collapsed so if there is any valid data in that axis the
        # valid mask is set. Thus the resulting shape is (lat, lon, time).
        self.valid_data = np.any(~np.ma.getmaskarray(tb), axis=3)

    @property
    def num_points(self) -> int:
        """Return the number of valid points."""
        return np.count_nonzero(self.valid_data)

    def compute_atmosphere(
        self, times: Sequence[int], era5_data: Era5DailyData
    ) -> HourlyRtm:
        """Compute the atmospheric RTM for the input ERA5 data."""
        # Extract the ERA5 data corresponding to the valid data mask. Note that
        # the ERA5 profile data is dimensioned as (time, lats, lons, levels),
        # but the valid_mask has dimensions (lats, lons, time). So the
        # valid_mask is reshaped to match the ERA5 data.
        valid_at_times = np.moveaxis(self.valid_data[:, :, times], 2, 0)
        atmo_results = rtm.compute(
            era5_data.levels,
            era5_data.temperature[valid_at_times],
            era5_data.height[valid_at_times],
            era5_data.relative_humidity[valid_at_times],
            era5_data.liquid_content[valid_at_times],
            era5_data.surface_temperature[valid_at_times],
            era5_data.surface_height[valid_at_times],
            era5_data.surface_relative_humidity[valid_at_times],
            era5_data.surface_pressure[valid_at_times],
            REF_EIA,
            REF_FREQ,
        )

        # Reshape the vectorized RTM outputs to full arrays, filling in missing
        # values with NaNs
        num_lat = len(self.lat)
        num_lon = len(self.lon)
        num_time = len(times)
        num_freq = len(REF_FREQ)
        valid_out = self.valid_data[:, :, times]
        col_water_vapor = np.full((num_lat, num_lon, num_time), np.nan)
        col_cloud_liquid = np.full((num_lat, num_lon, num_time), np.nan)
        transmissivity = np.full((num_lat, num_lon, num_time, num_freq), np.nan)
        tb_up = np.full((num_lat, num_lon, num_time, num_freq), np.nan)
        tb_down = np.full((num_lat, num_lon, num_time, num_freq), np.nan)

        col_water_vapor[valid_out] = era5_data.columnar_water_vapor[valid_at_times]
        col_cloud_liquid[valid_out] = era5_data.columnar_cloud_liquid[valid_at_times]
        transmissivity[valid_out] = atmo_results.tran
        tb_up[valid_out] = atmo_results.tb_up
        tb_down[valid_out] = atmo_results.tb_down
        return HourlyRtm(
            self.lat,
            self.lon,
            times,
            cast(Sequence[np.float32], col_water_vapor),
            cast(Sequence[np.float32], col_cloud_liquid),
            cast(Sequence[np.float32], transmissivity),
            cast(Sequence[np.float32], tb_up),
            cast(Sequence[np.float32], tb_down),
        )

    def append_results(self, hours: Sequence[int], data: HourlyRtm) -> None:
        """Update the daily TB file in-place with the RTM results."""
        with Dataset(self.tb_path, "a") as f:
            _ensure_rtm_vars(f)
            s_3d = np.s_[:, :, hours]
            s_4d = np.s_[:, :, hours, :]
            f.variables["columnar_water_vapor"][s_3d] = data.columnar_water_vapor
            f.variables["columnar_cloud_liquid"][s_3d] = data.columnar_cloud_liquid
            f.variables["transmissivity"][s_4d] = data.transmissivity
            f.variables["upwelling_tb"][s_4d] = data.tb_up
            f.variables["downwelling_tb"][s_4d] = data.tb_down


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

    ERA5 data, both on the surface and as profiles, is required.

    If the resampled TB file doesn't exist, a `FileNotFound` error will be
    raised.

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

    downloader.download_day(current_day, verbose)
    era5_path = downloader.out_dir

    # Process only one hour at a time. This assumes there are 24 hours in the
    # file...but that may not be correct.
    for hour in range(24):
        if verbose:
            print(f"Processing hour {hour+1}/24")

        era5_data = read_era5_data(
            era5_path / f"era5_surface_{current_day.isoformat()}.nc",
            era5_path / f"era5_levels_{current_day.isoformat()}.nc",
            (hour,),
            verbose,
        )

        if verbose:
            print("Computing RTM")
        atmosphere_results = daily_data.compute_atmosphere((hour,), era5_data)

        if verbose:
            print(f"Appending results to: {daily_data.tb_path}")
        daily_data.append_results((hour,), atmosphere_results)


if __name__ == "__main__":
    # TODO: use argparse

    if os.name == "nt":
        ACCESS_ROOT = Path("L:/access/")
    elif os.name == "posix":
        # ACCESS_ROOT = Path("/mnt/ops1p-ren/l/access")
        ACCESS_ROOT = Path("/mnt/wunda/access")

    try:
        cds_uid = os.environ["CDS_UID"]
        cds_api_key = os.environ["CDS_API_KEY"]
    except KeyError:
        print("CDS_UID and CDS_API_KEY environment variables need to be set")

    era5_dir = ACCESS_ROOT / "era5"
    downloader = Era5Downloader(cds_uid, cds_api_key, era5_dir)

    append_atmosphere_to_daily_ACCESS(
        date(2012, 7, 11), "amsr2", ACCESS_ROOT, downloader, verbose=True
    )
