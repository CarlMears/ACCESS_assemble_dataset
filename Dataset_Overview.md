# Microwave Brightness Temperatures Resampled to a regular Earth grid.

This document describes the format, data content and metadata for RSS Earth-grid microwave radiance datasets.  

## Overview and purpose
Microwave radiances measured by conically-scanning satellite radiometers are available from Remote Sensing Systems as L1B files with measurements arranged in the native "swath" format.  The measurement geometry is complex and it is not trivial to collocate these data between different satellites, or with other types of Earth data, such as model output.  To facilitate collocation of these microwave radiance with other sources of Earth data, we have accurately resampled the L1B radiances onto circular footprints on Earth-referenced latitude/longitude or polar grids.

## Satellites Available to Date
* AMSR2 on GCOM-W
    - 0.25 x 0.25 Latitude/Longitude Grid, 30 km Circular Footprints
    -  0.25 x 0.25 Latitude/Longitude Grid, 60 km Circular Footprints (under development)

## Dataset Structure
Data for each satellite/grid type are arranged into daily CF-compliant NetCDF4 files. The "base" or "resamp_tb" file contains resampled radiances, expressed as brightness temperatures (Tb's) for a number of microwave channels (observation frequency and polarization).  The brightess temperature variable for the global latitude/longitude grid has the following dimensions:

* Latitude (721)
* Longitude (1440)
* Hour (0..23)
* Channel

Each Hour=hh/Channel=ch slice is a map of the satellite observations for channel=ch than fall into the 1-hour time period starting at hh:00.  Any given slice contains a lot of missing data.  The effect of the missing data on file size is reduced by using the built-in compression in the NetCDF4.  The advantage of this data structure is that no measurements are discarded, and the time-ordered structure makes it easy find collocated data from other data sources.  The channels available vary between different sensor types and the size of the target footprint and are described by the "channels" dimension in the NetCDF4 file.

A more precise time-of-observation is contained in the "time" variable.  The variation of observation time for the different native footprints that are combined to construct the resampled footprint is on the order of 10's of seconds.  The reported time is the time of the actaul satellite observation with the maximum weight.

Each base file also contains land-fraction for each resampled footprint, constructed by averaging a MODIS-derived land/water mask over the resampled footprint.

The target footprints are circular gaussians with the FWHM diameter denoted by the "footprint size", e.g. a "30 km footprint" means a gaussian with FWHW of 30 km.

## Channels Available
* AMSR2, 30 km - 11,19,24,37,89 GHz (V and H polarization)
* AMSR2, 60 km - 6.9,7.2,11,19,24,37,89 GHz (V and H polarization)

    
## Collocated Ancillary Data
For each daily "base" file, we provide a number of ancillary data types on the same grid.  The exact nature of the spatial extent of these ancillary data depends on the data source.  Here we briefly describe the ancillary data currently available for each sensor type.  More detailed descriptions are available below.

| data type | description |
|-----------|-------------|
|IMERGE rain rate (mm/hr)| Resampled onto target footprint|
|ERA5 Skin Temperature (K)|ERA5 value regridded by ECMWF to footprint location|
|ERA5 Total Column Water Vapor (kg/m^2)|ERA5 value regridded by ECMWF to footprint location|
|ERA5 Total Column Cloud Water (kg/m^2)|ERA5 value regridded by ECMWF to footprint location|
|ERA5 Zonal Wind (m/s)|ERA5 value regridded by ECMWF to footprint location|
|ERA5 Meridional Wind (m/s)|ERA5 value regridded by ECMWF to footprint location|
|Atmospheric Parameters for each channel|Calculated using ERA5 profiles as input to the RSS RTM|
|- Transmissivity|Calculated using ERA5 profiles as input to the RSS RTM|
|- Upwelling Radiance (K)|Calculated using ERA5 profiles as input to the RSS RTM|
|- Downwelling Radiance (K)|Calculated using ERA5 profiles as input to the RSS RTM|

## Date range available
2012-07-02 to 2021-12-31

## Directory Structure and file name conventions
Files are locating in a folder for eac each inside a year - month - day structure e.g.
***/2013/05/04***

| data type | description |
|-----------|-------------|
|Resampled Satellite Brightness Temperatures|amsr2_resamp_tbs_2012_10_04.nc| 
|IMERGE rain rate (mm/hr)| amsr2_rain_rate_imerge_2012_10_04.nc|
|ERA5 Skin Temperature (K)|amsr2_skt_era5_2012_10_04.nc|
|ERA5 Total Column Water Vapor (kg/m^2)|amsr2_tcwv_era5_2012_10_04.nc|
|ERA5 Total Column Cloud Water (kg/m^2)|amsr2_tclw_era5_2012_10_04.nc|
|ERA5 Zonal Wind (m/s)|amsr2_u10n_era5_2012_10_04.nc|
|ERA5 Meridional Wind (m/s)|amsr2_v10n_era5_2012_10_04.nc|
|Atmospheric Parameters for each channel (Upwelling Radiance, Downwelling Radiance, Tranmissivity)|amsr2_atm_par_era5_2012_10_04.nc|
|- Transmissivity|Calculated using ERA5 profiles as input to the RSS RTM|
|- Upwelling Radiance (K)|Calculated using ERA5 profiles as input to the RSS RTM|
|- Downwelling Radiance (K)|Calculated using ERA5 profiles as input to the RSS RTM|

## Documentation
The methods to construct this data collection are described in two papers in preparation
* "A computationally-efficient for resampling microwave radianc-es from conical scanners to a regular Earth grid"
* "A multi-variate data collection of AMSR2 Satellite radiances and collocated ancillary variables on a regular Earth grid"









