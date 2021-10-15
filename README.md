# Assemble ACCESS Dataset

This set of scripts combines all the constituent components to finish assembling
the ACCESS dataset.

## Instructions

### Surface temperature

### Land fraction

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

To install the other required packages:
```
pip install numpy netCDF4
```

As an example, for updating the AMSR2 file on 2012-07-11:

```bash
python add_atmosphere_to_ACCESS_output.py L:/access 2012-07-11 amsr2 --user $CDS_UID --key $CDS_API_KEY
```

The CDS credentials can also be given as environment variables:

```bash
export CDS_UID=XXX
export CDS_API_KEY=XXX
python add_atmosphere_to_ACCESS_output.py L:/access 2012-07-11 amsr2
```