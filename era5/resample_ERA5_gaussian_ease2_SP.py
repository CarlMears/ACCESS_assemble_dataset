from AMSR2Antenna.AMSR2_Antenna_Gain import target_gain

import numpy as np
import xarray as xr
import pyproj as proj
import warnings

from polar_grids import NSIDC_ease2_grids

# using this routine ensures that the target footprint for the land fraction is the same
# as is used for the resampling target pattern.  It is not quite a gaussian, but close.

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    restart = False

    footprint_diameter_km = 30.0
    weight_map = np.full((720, 720, 100), np.nan, dtype=np.float32)
    ilon_map = np.full((720, 720, 100), -9999, dtype=np.int16)
    ilat_map = np.full((720, 720, 100), -9999, dtype=np.int16)

    ease2_grid_25km_north = NSIDC_ease2_grids(pole="north", resolution="25km")
    ease2_grid_25km_south = NSIDC_ease2_grids(pole="south", resolution="25km")

    for iy in np.arange(0, 720):
        for ix in np.arange(0, 720):
            latitude0 = ease2_grid_25km_south.latitude[iy, ix]
            longitude0 = ease2_grid_25km_south.longitude[iy, ix]
            if longitude0 < 0.0:
                longitude0 = longitude0 + 360
            if not np.isfinite(latitude0 + longitude0):
                continue
            if np.abs(latitude0) > 88.0:
                continue

            lat_index = int(np.round((latitude0 + 90.0) / 0.25))
            lon_index = int(np.round(longitude0 / 0.25))

            # figure out indices of the ERA5 submap needed

            half_numy = 10  # eventually, this will be a function of the footprint size
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

            # could use code here to wrap over the pole but for now,
            # only do -88 to 88

            # The data is in lat/lon coordinates.  To integrate
            # over the footprint, we need it to be
            # in x,y (kilometer) coordinates.
            # Fortunately, the pyproj projection
            # package does just this.
            # The warnings from this package are
            # suppressed.  Something about the way things are called
            # here are deprecated (including internally
            # to the package).  This needs to be cleaned up at
            # some point, but for now it works.
            #
            # construct a local projection, centered at the center of the target
            # gaussian.

            crs_wgs = proj.Proj(
                init="epsg:4326"
            )  # assuming you're using WGS84 geographic
            cust = proj.Proj(
                (
                    f"+proj=aeqd +lat_0={latitude0} +lon_0={longitude0} "
                    f"+datum=WGS84 +units=m"
                )
            )

            # apply the local projection to obtain the mesh points in meters.
            xv, yv = proj.transform(crs_wgs, cust, sub_map_lon, sub_map_lat)

            # convert to km
            yv = yv / 1000.0
            xv = xv / 1000.0

            # calculate the footprint - don't care about normalization --
            # normalization is
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

            weight_map[iy, ix, :] = wt_out
            ilon_map[iy, ix, :] = ilon_out.astype(np.int16)
            ilat_map[iy, ix, :] = ilat_out.astype(np.int16)

        sample_array = np.arange(0, 100)
        ix_array = np.arange(0, 720)
        iy_array = np.arange(0, 720)

        resample_wts_DS = xr.Dataset(
            data_vars=dict(
                weights=(["Y", "X", "Sample"], weight_map),
                latitude=(["Y", "X"], ease2_grid_25km_north.latitude),
                longitude=(["Y", "X"], ease2_grid_25km_north.longitude),
                lat_index=(["Latitude", "Longitude", "Sample"], ilat_map),
                lon_index=(["Latitude", "Longitude", "Sample"], ilon_map),
            ),
            coords={"Y": iy_array, "X": ix_array, "Sample": sample_array},
        )

        nc_file = (
            f"L:/access/era5/resample_weights/resamp_wts_ease2_25km_SP"
            f"_{int(footprint_diameter_km)}km.nc"
        )
        encoding = {"weights": {"zlib": True, "complevel": 4}}
        resample_wts_DS.to_netcdf(nc_file, encoding=encoding)
        print(f"Finished Row = {iy}")
        print(f"Wrote: {nc_file}")

    print()
