import sys

sys.path.append(
    "m:/job_access/python/resample_wts/AMSR2"
)  # contains the target footprint pattern code

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pyproj as proj
import warnings

# from bin_ndarray import bin_ndarray
# from gauss_kruger import local_geogauss

from AMSR2_Antenna_Gain import target_gain

# using this routine ensures that the target footprint for the land fraction is the same
# as is used for the resampling target pattern.  It is not quite a gaussian, but close.

warnings.filterwarnings("ignore")

if __name__ == "__main__":

    warnings.filterwarnings("ignore")
    restart = False

    footprint_diameter_km = 30.0
    weight_map = np.full((721, 1440, 100), np.nan, dtype=np.float32)
    ilon_map = np.full((721, 1440, 100), -9999, dtype=np.int16)
    ilat_map = np.full((721, 1440, 100), -9999, dtype=np.int16)

    for latitude0 in np.arange(-88.0, 88.1, 0.25):
        lat_index = int(np.round((latitude0 + 90.0) / 0.25))
        latitude0_radians = np.deg2rad(latitude0)
        print(latitude0)
        for longitude0 in np.arange(0.0, 360.0, 0.25):
            lon_index = int(np.round(longitude0 / 0.25))
            longitude0_radians = np.deg2rad(longitude0)

            # figure out indices of the ERA5 submap needed

            half_numy = 10
            # eventually, this will be a function of the footprint size
            half_numx = int(np.ceil(half_numy / np.cos(np.deg2rad(latitude0))))

            numy = 1 + 2 * half_numy
            numx = 1 + 2 * half_numx

            sub_map_weight = np.zeros((numy, numx), dtype=np.float32)
            sub_map_ilat = np.zeros((numy, numx), dtype=np.int32)
            sub_map_ilon = np.zeros((numy, numx), dtype=np.int32)

            ilat_array = lat_index + np.arange(-half_numy, half_numy + 1)
            ilat_edge_array = -0.5 + lat_index + np.arange(-half_numy, half_numy + 2)
            ilon_array = lon_index + np.arange(-half_numx, half_numx + 1)
            ilon_edge_array = -0.5 + lon_index + np.arange(-half_numx, half_numx + 2)

            # wrap the longitudes
            ilon_array[ilon_array < 0] = ilon_array[ilon_array < 0] + 1440
            ilon_array[ilon_array > 1439] = ilon_array[ilon_array > 1439] - 1440

            ilon_edge_array[ilon_edge_array < 0] = (
                ilon_edge_array[ilon_edge_array < 0] + 1440
            )
            ilon_edge_array[ilon_edge_array > 1439] = (
                ilon_edge_array[ilon_edge_array > 1439] - 1440
            )

            for i in np.arange(numx):
                sub_map_ilat[:, i] = ilat_array
            for i in np.arange(numy):
                sub_map_ilon[i, :] = ilon_array

            sub_map_lat = 0.25 * (sub_map_ilat - 360)
            sub_map_lon = 0.25 * (sub_map_ilon)

            sub_map_edge_ilat = np.zeros((numy + 1, numx + 1), dtype=np.float32)
            sub_map_edge_ilon = np.zeros((numy + 1, numx + 1), dtype=np.float32)

            for i in np.arange(numx + 1):
                sub_map_edge_ilat[:, i] = ilat_edge_array
            for i in np.arange(numy + 1):
                sub_map_edge_ilon[i, :] = ilon_edge_array

            sub_map_edge_lat = 0.25 * (sub_map_edge_ilat - 360)
            sub_map_edge_lon = 0.25 * (sub_map_edge_ilon)

            # could use code here to wrap over the pole but for now, only do -88 to 88

            # The data is in lat/lon coordinates.  To integrate over the footprint, we need it to be
            # in x,y (kilometer) coordinates.  Fortunately, the pyproj projection package does just this.
            # The warnings from this package are suppressed.  Something about the way things are called
            # here are deprecated (including internally to the package).  This needs to be cleaned up at
            # some point, but for now it works.
            #
            # construct a local projection, centered at the center of the target gaussian.

            crs_wgs = proj.Proj(
                init="epsg:4326"
            )  # assuming you're using WGS84 geographic
            cust = proj.Proj(
                f"+proj=aeqd +lat_0={latitude0} +lon_0={longitude0} +datum=WGS84 +units=m"
            )

            # apply the local projection to obtain the mesh points in meters.
            xv, yv = proj.transform(crs_wgs, cust, sub_map_lon, sub_map_lat)

            # convert to km
            yv = yv / 1000.0
            xv = xv / 1000.0

            # calculate the footprint - don't care about normalization -- normalization is
            # done "by hand" later.

            dist_from_center = np.sqrt(np.square(xv) + np.square(yv))
            # this is the sample routine used to calculate the target footprints
            footprint_amplitude = target_gain(
                dist_from_center, diameter_in_km=footprint_diameter_km
            )

            indices_y, indices_x = np.unravel_index(
                np.argsort(footprint_amplitude, axis=None), footprint_amplitude.shape
            )
            indices_y = indices_y[::-1][0:100]
            indices_x = indices_x[::-1][0:100]

            ilat_out = sub_map_ilat[indices_y, indices_x]
            ilon_out = sub_map_ilon[indices_y, indices_x]
            wt_out = footprint_amplitude[indices_y, indices_x]
            wt_out = wt_out / np.sum(wt_out)

            weight_map[lat_index, lon_index, :] = wt_out
            ilon_map[lat_index, lon_index, :] = ilon_out.astype(np.int16)
            ilat_map[lat_index, lon_index, :] = ilat_out.astype(np.int16)

        lon_array = np.arange(0.0, 360.0, 0.25)
        lat_array = np.arange(-90.0, 90.1, 0.25)
        sample_array = np.arange(0, 100)
        land_mask_DS = xr.Dataset(
            data_vars=dict(
                weights=(["Latitude", "Longitude", "Sample"], weight_map),
                lat_index=(["Latitude", "Longitude", "Sample"], ilat_map),
                lon_index=(["Latitude", "Longitude", "Sample"], ilon_map),
            ),
            coords={
                "Latitude": lat_array,
                "Longitude": lon_array,
                "Sample": sample_array,
            },
        )

        nc_file = f"L:/access/era5/resample_weights/resamp_wts_1440_721_rect_{int(footprint_diameter_km)}km.nc"
        encoding = {"weights": {"zlib": True, "complevel": 4}}
        land_mask_DS.to_netcdf(nc_file, encoding=encoding)
        print(f"Finished Lat = {latitude0}")
        print(f"Wrote: {nc_file}")

    print()
