import datetime
from typing import Sequence, Union, Optional
from pathlib import Path


def common_global_attributes_access(date : datetime.datetime, 
                                    version: str = 'v00r00') -> dict:

    day_boundary = datetime.datetime.combine(date, datetime.time())
    start_date = day_boundary - datetime.timedelta(minutes=30)
    end_date = day_boundary + datetime.timedelta(minutes=1410.0)

    attrs_glob = dict(
        Conventions="CF-1.8",
        standard_name_vocabulary="CF Standard Name Table (v78, 21 September 2021)",
        product_version=version,
        date_issued="2021-10-01",
        keywords_vocabulary=(
            "NASA Global Change Master Directory (GCMD) "
            "Earth Science Keywords, Version 6.0"
        ),
        cdm_data_type="Grid",
        program=(
            "NASA ACCESS-0031 > Machine Learning Datasets "
            "for the Earth's Natural Microwave Emission"
        ),
        date_created=datetime.datetime.now().isoformat(),
        creator_name="Carl Mears",
        creator_url="http://www.remss.com/ ",
        creator_email="mears@remss.com",
        institution="Remote Sensing Systems",
        processing_level="NASA Level 4",
        references="None",
        geospatial_lat_min=-90.0,
        geospatial_lat_max=90.0,
        geospatial_lon_min=0.0,
        geospatial_lon_max=359.9999,
        geospatial_lat_units="degrees_north",
        geospatial_lon_units="degrees_east",
        spatial_resolution="30 km X 30 km",
        time_coverage_start=start_date.isoformat(),
        time_coverage_end=end_date.isoformat(),
        time_coverage_duration="P24H",
        license="This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.",
        contributor_name="Frank Wentz, Carl Mears",
        contributor_role=(
            "Principal investigator and originator of "
            "input/source or antenna temperature data, "
            "Principle investigator of ACCESS project which resamples "
            "RSS native brightness temperatures and ancillary to a "
            "fixed Earth grid"
        ),
    )

    return attrs_glob


def resamp_tb_attributes_access(
                        satellite: str, 
                        version="v00r00"
                                ):

    attrs_glob_tb = (
        dict(
            summary=(
                "Remote Sensing Systems (RSS) Resampled brightness temperature; "
                "intercalibrated and homogenized brightness temperature "
                "polar-orbiting resampled to a regular Earth grid"
            ),
            keywords=(
                "EARTH SCIENCE > SPECTRAL/ENGINEERING > MICROWAVE > BRIGHTNESS TEMPERATURE"
            ),
            title=f"Resampled {satellite} brightness temperatures",
            version=version,
            history=(
                datetime.datetime.now().isoformat()
                + " Created Resampled Brightness Temperature from RSS AMSR2 L1A data"
            ),
            platform="GCOM-W1, JAXA",
            sensor="AMSR2 > Advanced Microwave Scanning Radiometer 2",
        ),
    )
    return attrs_glob_tb


def skt_era5_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_rr = dict(
        summary=(
            f"Skin Temperature corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=(
            "EARTH SCIENCE > LAND SURFACE > SURFACE THERMAL PROPERTIES > SKIN TEMPERATURE, "
            "EARTH SCIENCE > OCEANS > OCEAN TEMPERATURE > SEA SURFACE TEMPERATURE > SEA SURFACE SKIN TEMPERATURE"
        ),
        title=f"ERA5 skin temperature, resampled to circular gaussian footprints 0.25 degree Earth grid",
        history=f"{datetime.datetime.now()} Created from ERA5 data downloaded from the Copernicus Climate Data Store",
    )
    return attrs_glob_rr


def rr_imerg_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_rr = dict(
        summary=(
            f"Rain Rate corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=("EARTH SCIENCE > ATMOSPHERE > PRECIPITATION > PRECIPITATION RATE"),
        title=f"IMERG Rain Rate, resampled to circular gaussian footprints 0.25 degree Earth grid",
        history=f"{datetime.datetime.now()} Created from on-line IMERG data",
    )
    return attrs_glob_rr


def tcwv_era5_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_tcwv = dict(
        summary=(
            f"Total Column Water Vapor corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=(
            "EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WATER VAPOR > WATER VAPOR INDICATORS > TOTAL PRECIPITABLE WATER"
        ),
        title=f"Total Column Water Vapor from ERA5 on a 0.25 degree Earth grid, time interpolated to satellite overpass time",
        history=f"{datetime.datetime.now()} Created from ERA5 data downloaded from the Copernicus Climate Data Store",
    )
    return attrs_glob_tcwv


def tclw_era5_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_tclw = dict(
        summary=(
            f"Total Column Cloud Water corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=(
            "EARTH SCIENCE > ATMOSPHERE > CLOUDS > CLOUD MICROPHYSICS > CLOUD LIQUID WATER"
        ),
        title=f"Total Column Cloud Water from ERA5 on a 0.25 degree Earth grid, time interpolated to satellite overpass time",
        history=f"{datetime.datetime.now()} Created from ERA5 data downloaded from the Copernicus Climate Data Store",
    )
    return attrs_glob_tclw


def u10n_era5_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_u10n = dict(
        summary=(
            f"Eastward wind corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=(
            "EARTH SCIENCE > OCEANS > OCEAN WINDS > SURFACE WINDS",
            "EARTH SCIENCE > ATMOSHERE > ATMOSPHERIC WINDS > SURFACE WINDS",
        ),
        title=f"10m neutral zonal winds from ERA5 on a 0.25 degree Earth grid, time interpolated to satellite overpass time",
        history=f"{datetime.datetime.now()} Created from ERA5 data downloaded from the Copernicus Climate Data Store",
    )
    return attrs_glob_u10n


def v10n_era5_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_v10n = dict(
        summary=(
            f"Northward wind corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=(
            "EARTH SCIENCE > OCEANS > OCEAN WINDS > SURFACE WINDS",
            "EARTH SCIENCE > ATMOSHERE > ATMOSPHERIC WINDS > SURFACE WINDS",
        ),
        title=f"10m neutral meridional winds from ERA5 on a 0.25 degree Earth grid, time interpolated to satellite overpass time",
        history=f"{datetime.datetime.now()} Created from ERA5 data downloaded from the Copernicus Climate Data Store",
    )
    return attrs_glob_v10n

def atm_pars_era5_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_atm = dict(
        summary=(
            f"Atmopspheric radiation parameters corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=(
            "EARTH SCIENCE > ATMOSHERE > ATMOSPHERIC RADIATION > SPECTRAL IRRADIANCE",
            "EARTH SCIENCE > ATMOSHERE > ATMOSPHERIC RADIATION > TRANSMITTANCE",
        ),
        title=f"Atmospheric radiation parameters computed from ERA5 profile on a 0.25 degree Earth grid, time interpolated to satellite overpass time",
        history=f"{datetime.datetime.now()} Created from ERA5 data downloaded from the Copernicus Climate Data Store using the RSS radiative transfer model",
    )
    return attrs_glob_atm

def tbup_era5_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_tbup = dict(
        summary=(
            f"Upwelling radiance (brightness temperature) corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=(
            "EARTH SCIENCE > ATMOSHERE > ATMOSPHERIC RADIATION > SPECTRAL IRRADIANCE",
        ),
        title=f"Upwelling Radiance computed from ERA5 profile on a 0.25 degree Earth grid, time interpolated to satellite overpass time",
        history=f"{datetime.datetime.now()} Created from ERA5 data downloaded from the Copernicus Climate Data Store using the RSS radiative transfer model",
    )
    return attrs_glob_tbup

def tbdown_era5_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_tbdown= dict(
        summary=(
            f"Downwelling radiance (brightness temperature) corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=(
            "EARTH SCIENCE > ATMOSHERE > ATMOSPHERIC RADIATION > SPECTRAL IRRADIANCE",
        ),
        title=f"Donwelling Radiance computed from ERA5 profile on a 0.25 degree Earth grid, time interpolated to satellite overpass time",
        history=f"{datetime.datetime.now()} Created from ERA5 data downloaded from the Copernicus Climate Data Store using the RSS radiative transfer model",
    )
    return attrs_glob_tbdown

def trans_era5_attributes_access(satellite:str, version="v00r00"):

    attrs_glob_trans= dict(
        summary=(
            f"Atmospheric Transmissivty corresponding to Remote Sensing Systems (RSS) Resampled brightness temperature for {satellite}; "
        ),
        keywords=(
            "EARTH SCIENCE > ATMOSHERE > ATMOSPHERIC RADIATION > TRANSMITTANCE",
        ),
        title=f"Donwelling Radiance computed from ERA5 profile on a 0.25 degree Earth grid, time interpolated to satellite overpass time",
        history=f"{datetime.datetime.now()} Created from ERA5 data downloaded from the Copernicus Climate Data Store using the RSS radiative transfer model",
    )
    return attrs_glob_tbdown
