# Assemble ACCESS Dataset

This set of scripts combines all the constituent components to finish assembling
the ACCESS dataset.

All processing has been moved to linux, though most of these scripts will also run on windows with the exception of add_wind_emiss_ACCESS.py

## Instructions

The first step is to assemble a daily brightness temperature file.  The other steps,
which append ancillary data to the daily file, can be performed in any order.  Because
as writing to netCDF files is done using file locking, the data append jobs can be run 
simultaneously.

These are the required dependencies to run the scripts:

```bash
# From PyPI
pip install numpy netCDF4 matplotlib cdsapi xarray Cartopy
# From RSS repos
pip install access-atmosphere --index-url http://gitlab.remss.com/api/v4/projects/68/packages/pypi/simple --trusted-host gitlab.remss.com
pip install git+http://gitlab.remss.com/Carl/plotting.git
pip install git+http://gitlab.remss.com/Carl/rss_lock.git
```

### Assemble measured brightness temperatures

The `make_daily_ACCESS_files.py` script assembles the daily files from orbit files that contain the 
resampled circular footprints. Each file contains 24 hour-long slices of data arranged in lat/lon maps
on a 0.25 degree by 0.25 degree grid.  This results in a lot of missing data - the effect of this on file 
size if minimized by writing using slightly lossy built-in netCDF4 compression.

**Positional Arguments:**
- `access_root`: Path to location of daily ACCESS files
- `temp_root`: Path to location to store temporary files
- `start_date`: first day to process, in YYYY-MM-DD format
- `end_date`: last day to process, in YYYY-MM-DD format
- `sensor`: name of the sensor - currently only 'amsr2' is supported
- `target_size`: diameter of resampled footprints in km
- `version`: version string to place in metadata

**Optional Arguments:**
- `--overwrite`: if set, process and write file even if a file for this day already exists

**Example Command:**
```
on linux

python make_daily_ACCESS_files.py /mnt/ops1p-ren/l/access/amsr2_out_test /mnt/ops1p-ren/l/access/_temp 2012-07-02 2012-07-31 amsr2 30 v01r00--overwrite 

```

### 2D variables from ERA5
The `add_ERA5_2D_vars_to_ACCESS_output.py` script adds selected 2D variables to ERA5 to the daily file.

Currently implemented variables are:
- skt: skin temperatures
- tcwv: total column water vapor
- tclw: total column liquid water
- u10n: 10m neutral stability zonal wind
- v10n: 10m neutral stability meridional wind

If not available locally, the data are automatically downloaded from ECMWF using
an active account on the [Climate Data Store](https://cds.climate.copernicus.eu/). 
The CDS UID and API key are required inputs to the script, either as environment 
variables or as arguments.

**Positional Arguments:**
- `access_root`: Path to location of daily ACCESS files
- `output_root`: Path to location to write files
- `temp_root`: Path to location to store temporary files
- `start_date`: first day to process, in YYYY-MM-DD format
- `end_date`: last day to process, in YYYY-MM-DD format
- `sensor`: name of the sensor - currently only 'amsr2' is supported
- `target_size`: diameter of resampled footprints in km
- `version`: version string to place in metadata
- `target_size`: diameter of resampled footprints in km
- `vars`: list of variables to include

**Optional Arguments:**
- `--overwrite`: if set, process and write file even if a file for this day already exists

**Example Command:**
```

on linux:

python add_ERA5_surface_temperature_to_ACCESS_output.py /mnt/ops1p-ren/l/access/amsr2_out_test /mnt/ops1p-ren/l/access/_temp 2012-07-02 2012-07-31 amsr2 30 v10r00 "skt tclw u10n v10n tcwv"

```
### Ocean Emissivity
The `add_wind_emiss_ACCESS.py` script adds ocean emissivity layer to the daily netCDF4 file.  The input wind vectors need to already exist in the access_root.  (Currently, these are winds from ERA5.) The calculation is done using Meissner-Wentz surface model, implemented by a python wrapper about lightly modified versions of geomod10b.f and geomod10c.f

**Positional Arguments:**
- access_root: Path to location of daily ACCESS files
- output_root: Path to location to write new files
- wind_source: source to use for wind data, currently must be era5
- temperautre_source: source to use for wind data, currently must be era5
- start_date: first day to process, in YYYY-MM-DD format
- end_date: last day to process, in YYYY-MM-DD format
- sensor: name of the sensor - currently only 'amsr2' is supported
- target size: size of the resampled footprint in km
- version: version of the dataset, e.g. v01r00

**Example Command:**
```
linux:

python add_wind_emiss_to_ACCESS.py /mnt/ops1p-ren/l/access/amsr2_out_test /mnt/ops1p-ren/l//access/amsr2_out_test era5 era5 2012-07-02 2012-07-31 amsr2 30 v01r00 --verbose
```


### Land fraction
The `add_land_fraction_to_ACCESS_output.py` script adds a land fraction layer to the daily netCDF4 file.  The preferred dataset,
"modis" is mostly from the MODIS land/water dataset, with a few areas filled in using NSIDC data.

**Positional Arguments:**
- access_root: Path to location of daily ACCESS files
- temp_root: Path to location to store temporary files
- start_date: first day to process, in YYYY-MM-DD format
- end_date: last day to process, in YYYY-MM-DD format
- sensor: name of the sensor - currently only 'amsr2' is supported
- target size: size of the resampled footprint in km
- version: version of the dataset, e.g. r01v00
- lf_version: chooses the version of the land fraction dataset to use.  Currently preferred option is "modis"

**Optional Arguments:**
- --overwrite: if set, process and write layer even if it already exists in this file.
- --verbose: if set, print more verbose informational messages

**Example Command:**
```
python add_land_fraction_to_ACCESS_output.py L:\access\amsr2_out_test L:\access\_temp 2012-07-02 2012-07-31 amsr2 30 v01r00 modis --verbose --overwrite
```

### Rain Rate
The `add_imerg_rain_Rate_to_ACCESS_output.py` script adds a rain rate to the dataset, dereived by resmapling the IMERG 3 hourly rain rate product.  The required data is downloaded if needed.

**Positional Arguments:**
- access_root: Path to location of daily ACCESS files
- temp_root: Path to location to store temporary files
- start_date: first day to process, in YYYY-MM-DD format
- end_date: last day to process, in YYYY-MM-DD format
- sensor: name of the sensor - currently only 'amsr2' is supported
- target size: size of the resampled footprint in km
- version (not implemented yet): version of the dataset, e.g. r01v00
**Optional Arguments:**
- --overwrite: if set, process and write layer even if it already exists in this file.

**Example Command:**
```
python add_imerg_rain_rate_fraction_to_ACCESS_output.py L:\access\amsr2_out_test L:\access\_temp 2012-07-02 2012-07-31 amsr2 30 --overwrite
```

### Atmosphere

The `add_atmosphere_to_ACCESS_output.py` script computes the atmospheric terms
(transmissivity, upwelling brightness temperature, downwelling brightness
temperature) required for valid resampled TB.

The atmospheric data is from ERA5 and it needs to be already downloaded

Once the ERA5 data is present, the script runs the atmospheric RTM for all valid input
data and then appends the results to the output ACCESS file. The RTM is
implemented as a Python extension written in Fortran.

The ERA5 downloading and RTM calculations are contained in a dependent package,
[`access-atmosphere`](http://gitlab.remss.com/access/atmospheric-rtm). It can be
installed using:

```
pip install access-atmosphere --index-url http://gitlab.remss.com/api/v4/projects/68/packages/pypi/simple --trusted-host gitlab.remss.com
```

(The `--trusted-host` permits pip to use HTTP instead of HTTPS and is required until we get a TLS certificate for `gitlab.remss.com`...)

Required dependencies such as `numpy` and `netCDF4` should be automatically installed when installing `access-atmosphere`.

As an example, for updating the AMSR2 file on 2012-07-11:

```bash
python add_atmosphere_to_ACCESS_output.py L:/access 2012-07-11 amsr2 --user $CDS_UID --key $CDS_API_KEY
```

Instead of command-line arguments, the CDS credentials can instead be given as environment variables:

```bash
export CDS_UID=XXX
export CDS_API_KEY=XXX
python add_atmosphere_to_ACCESS_output.py L:/access 2012-07-11 amsr2
```
 ### Attributes for the netCDF4 output files
 The constant attributes are now defined in .json files in
 ~~~
 M:\job_access\python\dataset_assembly\access_io\attr_define_json
 ~~~
 Some attributes still need to be constructed programatically, so are not in the .json

 Here is an example json file.
 ~~~json
 {
    "global":{
        "summary": "Modis Land Mask resampled to gaussian footprints",
        "keywords": "EARTH SCIENCE > LAND SURFACE > LAND USE/LAND COVER > LAND/OCEAN/ICE MASK",
        "title": "Modis Land Mask resampled to circulargaussian footprintsgaussian footprints 0.25 degree Earth grid",
        "history": "Created from on-line MODIS data from the CMR EarthData API"},
    "var":{
        "standard_name": "land fraction",
        "valid_min": "0.0",
        "valid_max": "1.0",
        "units": "dimensionless"
        }
}
~~~
The idea is that attributes that are added to the global attributes for the file are in the "global" part, and attributes for the specific variable are in the "var" part.  

### .json file naming conventions

The files should be named with the following conventions:

~~~python
filename = f"attr.{project}.{satellite}.{variable}.json"
~~~
In this case, "project" is set to "access".
"variable" can be set to a variable, a coordinate, or "common" for common attributes shared by all files.
These files are located in 
~~~
assemble-dataset/access_io/attr_define_json/
~~~
