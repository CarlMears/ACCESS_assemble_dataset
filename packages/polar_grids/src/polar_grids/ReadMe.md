# Polar Grid Module

This Python module provides functions for interpolating values from a polar map onto specific latitude and longitude coordinates. The module uses the os and pathlib libraries to determine the appropriate file paths based on the operating system.

## Installation

To use this module, make sure you have Python installed on your system. You can then import the module using the following statement:

```
python

import polar_grid
```

## Usage

The module provides the following functions:

- polar_stereo_interp_SP(polar_map, lats, lons

This function interpolates values given in a south pole polar stereographic polar map onto the specified latitude and longitude coordinates. 
### Parameters
- polar_map: input map data
- lats: array of latitude to be interpolated to
- lons: array of longitude to be interpolarter too

### Return Value
array containing the interpolated values.

-------------------
- polar_stereo_interp(polar_map, lats, lons)


Similar to polar_stereo_interp_SP, this function interpolates values from a polar map onto latitude and longitude coordinates. However, it uses a different set of parameters and calculations specific to a different projection for the North Pole. The polar_map, lats, and lons parameters have the same meaning as in polar_stereo_interp_SP.

~~~
polarstereo_fwd(lats, lons, r_e=6378.2730, e=0.081816153, std_parallel=70.0, lon_y=-45.0)
~~~

This function performs the forward transformation for the polar stereographic projection. It converts latitude and longitude coordinates to Cartesian coordinates on a polar stereographic map. The lats and lons parameters are arrays containing the latitude and longitude coordinates, respectively. Optional parameters r_e, e, std_parallel, and lon_y define the projection parameters and have default values.
polarstereo_fwd_SP(lats, lons, r_e=6378.2730, e=0.081816153, std_parallel=70.0, lon_y=180.0)

Similar to polarstereo_fwd, this function performs the forward transformation for a different polar stereographic projection. It converts latitude and longitude coordinates to Cartesian coordinates on a polar stereographic map. The parameters have the same meaning as in polarstereo_fwd.
File Paths

The module uses the os library to determine the operating system and sets the appropriate file paths for the polar grids. If the operating system is "nt" (Windows), the module assumes the polar grids are located in "M:/job_access/polar_grids/". If the operating system is "posix" (Unix-like systems), the module assumes the polar grids are located in "/mnt/m/job_access/polar_grids/". If the operating system is neither "nt" nor "posix", a ValueError is raised.

Make sure to adjust the file paths in the module based on your system configuration.
Dependencies

The module requires the numpy library to perform calculations and array operations.

To install numpy, use the following command:

pip install numpy

License

This module is provided under the MIT License. Feel free to use and modify it according to your needs.