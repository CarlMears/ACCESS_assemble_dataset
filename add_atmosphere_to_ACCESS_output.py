"""Run the RTM and add the atmospheric terms to the resampled TB file.

This depends on the atmospheric RTM package, access-atmosphere:
http://gitlab.remss.com/access/atmospheric-rtm
"""

import os
from datetime import date
from pathlib import Path
from typing import Any, Sequence

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


class DailyRtm:
    """RTM results for the entire day."""

    def __init__(
        self,
        lat: Sequence[np.float32],
        lon: Sequence[np.float32],
        time: Sequence[np.float32],
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
        self.col_water_vapor = np.full(shape_3d, np.nan, np.float32)
        self.col_cloud_liquid = np.full(shape_3d, np.nan, np.float32)
        self.transmissivity = np.full(shape_4d, np.nan, np.float32)
        self.tb_up = np.full(shape_4d, np.nan, np.float32)
        self.tb_down = np.full(shape_4d, np.nan, np.float32)

    def compute_atmosphere(
        self, era5_data: Era5DailyData, times: Sequence[int], valid_data: Any
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

    def append_results(self, data: DailyRtm) -> None:
        """Update the daily TB file in-place with the RTM results for the day."""
        with Dataset(self.tb_path, "a") as f:
            _ensure_rtm_vars(f)
            f.variables["columnar_water_vapor"][...] = data.col_water_vapor
            f.variables["columnar_cloud_liquid"][...] = data.col_cloud_liquid
            f.variables["transmissivity"][...] = data.transmissivity
            f.variables["upwelling_tb"][...] = data.tb_up
            f.variables["downwelling_tb"][...] = data.tb_down


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

    rtm_data = DailyRtm(daily_data.lat, daily_data.lon, daily_data.time)

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
        rtm_data.compute_atmosphere(era5_data, (hour,), daily_data.valid_data)

    if verbose:
        print(f"Appending results to: {daily_data.tb_path}")
    daily_data.append_results(rtm_data)


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
