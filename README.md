# Assemble ACCESS Dataset

This set of scripts combines all the constituent components to finish assembling
the ACCESS dataset.

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

**Optional Arguments:**
- `--overwrite`: if set, process and write file even if a file for this day already exists
- `--plot_map`: if set, plot an example map as a debugging feature
- `--verbose`: if set, print more verbose informational messages

**Example Command:**
```
python make_daily_ACCESS_files.py L:\access\amsr2_out_test L:\access\_temp 2012-07-02 2012-07-31 amsr2 --overwrite --verbose
```

### Surface temperature
The `add_ERA5_surface_temperature_to_ACCESS_output.py` script adds surface skin temperature from 
ERA5 to the daily file.
If not available locally, the data are automatically downloaded from ECMWF using
an active account on the [Climate Data Store](https://cds.climate.copernicus.eu/). 
The CDS UID and API key are required inputs to the script, either as environment 
variables or as arguments.

**Positional Arguments:**
- `access_root`: Path to location of daily ACCESS files
- `temp_root`: Path to location to store temporary files
- `start_date`: first day to process, in YYYY-MM-DD format
- `end_date`: last day to process, in YYYY-MM-DD format
- `sensor`: name of the sensor - currently only 'amsr2' is supported

**Optional Arguments:**
- `--verbose`: if set, print more verbose informational messages

**Example Command:**
```
python add_ERA5_surface_temperature_to_ACCESS_output.py L:\access\amsr2_out_test L:\access\_temp 2012-07-02 2012-07-31 amsr2 --verbose
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
- version: chooses the version of the land fraction dataset to use.  Currently preferred option is "modis"

**Optional Arguments:**
- --overwrite: if set, process and write layer even if it already exists in this file.
- --verbose: if set, print more verbose informational messages

**Example Command:**
```
python add_land_fraction_to_ACCESS_output.py L:\access\amsr2_out_test L:\access\_temp 2012-07-02 2012-07-31 amsr2 modis --verbose --overwrite
```

### Atmosphere

The `add_atmosphere_to_ACCESS_output.py` script computes the atmospheric terms
(transmissivity, upwelling brightness temperature, downwelling brightness
temperature) required for valid resampled TB.

The atmospheric data is from ERA5 and it downloads the necessary surface and
profile data if needed. (Note that this will usually take quite some time.) In
order to download data, an active account on the [Climate Data
Store](https://cds.climate.copernicus.eu/) needs to be used. The CDS UID
and API key are required inputs to the script, either as environment variables
or as arguments.

Once the ERA5 data is present, it runs the atmospheric RTM for all valid input
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
