"""Run the RTM and add the atmospheric terms to the resampled TB file.

This depends on the atmospheric RTM package, access-atmosphere:
http://gitlab.remss.com/access/atmospheric-rtm
"""

import os
from datetime import date
from pathlib import Path
from typing import Sequence, NamedTuple

import numpy as np
from numpy.typing import ArrayLike
from netCDF4 import Dataset

from access_io.access_output import get_access_output_filename, write_daily_tb_netcdf
from access_atmosphere.download import Era5Downloader
from access_atmosphere.era5 import Era5DailyData, read_era5_data

# from resampled_tbs.read_resampled_orbit import read_resampled_tbs

# Reference frequencies (in GHz) to use
REF_FREQ = np.array([6.8, 10.7, 18.7, 23.8, 37.0, 89.0], np.float32)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
REF_EIA = np.array([53.0, 53.0, 53.0, 53.0, 53.0, 53.0], np.float32)


# NUM_LATS = 721
# NUM_LONS = 1440
# NUM_HOURS = 24
# NUM_CHANNELS = 14  # all possible AMSR2 channels
# AVAILABLE_CHANNELS = [
#     "time",
#     "6V",
#     "6H",
#     "7V",
#     "7H",
#     "11V",
#     "11H",
#     "19V",
#     "19H",
#     "24V",
#     "24H",
#     "37V",
#     "37H",
#     "89V",
#     "89H",
# ]


class RtmResults(NamedTuple):
    """Results after running the atmospheric RTM."""

    # Note, the type annotations here are sub-optimal until numpy has better
    # typing support. They are all 1d or 2d numpy arrays.

    # Total column water vapor in kg/m^2, dimensioned as (num_points, ).
    columnar_water_vapor: Sequence[np.float32]

    # Total column cloud liquid water in kg/m^2, dimensioned as (num_points, )
    columnar_cloud_liquid: Sequence[np.float32]

    # Atmospheric transmissivity, unitless, dimensioned as (num_points, freq).
    transmissivity: Sequence[np.float32]

    # Atmospheric upwelling brightness temperature in K, dimensioned as
    # (num_points, freq).
    tb_up: Sequence[np.float32]

    # Atmospheric downwelling brightness temperature in K, dimensioned as
    # (num_points, freq).
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

        # This boolean array is True whereever there is valid data
        self.valid_data = ~np.ma.getmaskarray(tb)

    @property
    def num_points(self) -> int:
        """Return the number of valid points."""
        return np.count_nonzero(self.valid_data)

    def compute_atmosphere(self, era5_data: Era5DailyData) -> RtmResults:
        raise NotImplementedError

    def append_results(
        self, hours: Union[SupportsIndex, slice], atmosphere_results: RtmResults
    ) -> None:
        with Dataset(self.tb_path, "a") as f:
            _ensure_rtm_vars(f)
            # TODO: write variables
            ...
        raise NotImplementedError


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
    era5_data = read_era5_data(
        era5_path / f"era5_surface_{current_day.isoformat()}.nc",
        era5_path / f"era5_levels_{current_day.isoformat()}.nc",
        verbose,
    )

    atmosphere_results = daily_data.compute_atmosphere(era5_data)
    daily_data.append_results(atmosphere_results)


if __name__ == "__main__":
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
