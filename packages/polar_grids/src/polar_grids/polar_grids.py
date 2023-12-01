import os
from pathlib import Path

if os.name == "nt":
    polar_grid_root = Path("M:/job_access/polar_grids/")
elif os.name == "posix":
    polar_grid_root = Path("/mnt/m/job_access/polar_grids/")
else:
    raise ValueError(f"os.name = {os.name} not valid for polar_grid.py")


def polar_stereo_interp_SP(polar_map, lats, lons):
    """This routine interpolates values given in polar map onto the locations
    given by lats and lons"""

    import numpy as np

    x, y = polarstereo_fwd_SP(lats, lons)
    xscale = x * 0.04
    yscale = y * 0.04

    x_lower = np.floor(xscale - 0.5).astype(np.int32)
    y_lower = np.floor(yscale - 0.5).astype(np.int32)

    x_upper = x_lower + 1
    y_upper = y_lower + 1

    wt_x_upper = (xscale - 0.5) - x_lower
    wt_y_upper = (yscale - 0.5) - y_lower

    interp_ok = np.all(
        [
            ((y_lower + 158) >= 0),
            ((y_upper + 158) <= 331),
            ((x_lower + 158) >= 0),
            ((x_upper + 158) <= 315),
        ],
        axis=(0),
    )
    interp_value = np.zeros(lats.shape)

    wt_upper_right = wt_y_upper[interp_ok] * wt_x_upper[interp_ok]
    wt_upper_left = wt_y_upper[interp_ok] * (1.0 - wt_x_upper)[interp_ok]
    wt_lower_right = (1.0 - wt_y_upper[interp_ok]) * wt_x_upper[interp_ok]
    wt_lower_left = (1.0 - wt_y_upper[interp_ok]) * (1.0 - wt_x_upper[interp_ok])

    upper_right = polar_map[y_upper[interp_ok] + 158, x_upper[interp_ok] + 158]
    upper_left = polar_map[y_upper[interp_ok] + 158, x_lower[interp_ok] + 158]
    lower_right = polar_map[y_lower[interp_ok] + 158, x_upper[interp_ok] + 158]
    lower_left = polar_map[y_lower[interp_ok] + 158, x_lower[interp_ok] + 158]

    total = np.zeros((upper_right.size))
    total_wt = np.zeros((upper_right.size))

    ok_finite = np.isfinite(upper_right)
    total[ok_finite] = (
        total[ok_finite] + upper_right[ok_finite] * wt_upper_right[ok_finite]
    )
    total_wt[ok_finite] = total_wt[ok_finite] + wt_upper_right[ok_finite]

    ok_finite = np.isfinite(upper_left)
    total[ok_finite] = (
        total[ok_finite] + upper_left[ok_finite] * wt_upper_left[ok_finite]
    )
    total_wt[ok_finite] = total_wt[ok_finite] + wt_upper_left[ok_finite]

    ok_finite = np.isfinite(lower_right)
    total[ok_finite] = (
        total[ok_finite] + lower_right[ok_finite] * wt_lower_right[ok_finite]
    )
    total_wt[ok_finite] = total_wt[ok_finite] + wt_lower_right[ok_finite]

    ok_finite = np.isfinite(lower_left)
    total[ok_finite] = (
        total[ok_finite] + lower_left[ok_finite] * wt_lower_left[ok_finite]
    )
    total_wt[ok_finite] = total_wt[ok_finite] + wt_lower_left[ok_finite]

    interp_value[interp_ok] = total / total_wt

    return interp_value


def polar_stereo_interp(polar_map, lats, lons):
    """This routine interpolates values given in polar map onto the locations
    given by lats and lons"""

    import numpy as np

    x, y = polarstereo_fwd(lats, lons)
    xscale = x * 0.04
    yscale = y * 0.04

    x_lower = np.floor(xscale - 0.5).astype(np.int32)
    y_lower = np.floor(yscale - 0.5).astype(np.int32)

    x_upper = x_lower + 1
    y_upper = y_lower + 1

    wt_x_upper = (xscale - 0.5) - x_lower
    wt_y_upper = (yscale - 0.5) - y_lower

    interp_ok = np.all(
        [
            ((y_lower + 234) >= 0),
            ((y_upper + 234) <= 447),
            ((x_lower + 154) >= 0),
            ((x_upper + 154) <= 303),
        ],
        axis=(0),
    )
    interp_value = np.zeros(lats.shape)

    wt_upper_right = wt_y_upper[interp_ok] * wt_x_upper[interp_ok]
    wt_upper_left = wt_y_upper[interp_ok] * (1.0 - wt_x_upper)[interp_ok]
    wt_lower_right = (1.0 - wt_y_upper[interp_ok]) * wt_x_upper[interp_ok]
    wt_lower_left = (1.0 - wt_y_upper[interp_ok]) * (1.0 - wt_x_upper[interp_ok])

    upper_right = polar_map[y_upper[interp_ok] + 214, x_upper[interp_ok] + 154]
    upper_left = polar_map[y_upper[interp_ok] + 214, x_lower[interp_ok] + 154]
    lower_right = polar_map[y_lower[interp_ok] + 214, x_upper[interp_ok] + 154]
    lower_left = polar_map[y_lower[interp_ok] + 214, x_lower[interp_ok] + 154]

    total = np.zeros((upper_right.size))
    total_wt = np.zeros((upper_right.size))

    ok_finite = np.isfinite(upper_right)
    total[ok_finite] = (
        total[ok_finite] + upper_right[ok_finite] * wt_upper_right[ok_finite]
    )
    total_wt[ok_finite] = total_wt[ok_finite] + wt_upper_right[ok_finite]

    ok_finite = np.isfinite(upper_left)
    total[ok_finite] = (
        total[ok_finite] + upper_left[ok_finite] * wt_upper_left[ok_finite]
    )
    total_wt[ok_finite] = total_wt[ok_finite] + wt_upper_left[ok_finite]

    ok_finite = np.isfinite(lower_right)
    total[ok_finite] = (
        total[ok_finite] + lower_right[ok_finite] * wt_lower_right[ok_finite]
    )
    total_wt[ok_finite] = total_wt[ok_finite] + wt_lower_right[ok_finite]

    ok_finite = np.isfinite(lower_left)
    total[ok_finite] = (
        total[ok_finite] + lower_left[ok_finite] * wt_lower_left[ok_finite]
    )
    total_wt[ok_finite] = total_wt[ok_finite] + wt_lower_left[ok_finite]

    interp_value[interp_ok] = total / total_wt

    return interp_value


def polarstereo_fwd(
    lats, lons, r_e=6378.2730, e=0.081816153, std_parallel=70.0, lon_y=-45.0
):
    #
    #   The National Snow and Ice Data Center (NSIDC) and Scientific Committee
    #   on Antarctic Research (SCAR) use a version of the polar stereographic
    #   projection.
    #
    #   Equations from: Map Projections - A Working manual - by J.P. Snyder. 1987
    #   http://kartoweb.itc.nl/geometrics/Publications/Map!20Projections!20-!20A!20Working!20manual!20-!20by!20J.P.!20Snyder.pdf
    #   P 160-161
    #   See the section on Polar Stereographic, with a south polar aspect and
    #   known phi_c not at the pole.
    #
    import numpy as np

    with np.errstate(invalid="ignore"):
        # constants for NSIDC projection
        xphi_c = std_parallel  # standard parallel, latitude of true scale in degrees
        xlambda_0 = lon_y  # meridian along positive Y axis

        # internal variables

        # convert to radians
        phi = np.pi / 180.0 * lats
        phi_c = np.pi / 180.0 * xphi_c * np.ones_like(phi)
        phi_c[phi < 0.0] = -phi_c[phi < 0]

        rlambda = np.pi / 180.0 * lons
        lambda_0 = np.pi / 180.0 * xlambda_0 * np.ones_like(phi)

        #   pm=-1 !plus or minus, north lat. or south
        phi[phi_c < 0.0] = -phi[phi_c < 0.0]
        rlambda[phi_c < 0.0] = -rlambda[phi_c < 0.0]
        lambda_0[phi_c < 0.0] = -lambda_0[phi_c < 0.0]
        phi_c[phi_c < 0.0] = -phi_c[phi_c < 0.0]
        # Equation numbers are from Snyder 1987 "Map Projections a Working Manual"
        # A copy of the paper is in
        # X:\DOCUMENT\Journal_Papers\Snyder_1987_projections_USGS.pdf

        t = np.tan(np.pi / 4.0 - phi / 2.0) / np.power(
            ((1 - e * np.sin(phi)) / (1 + e * np.sin(phi))), (e / 2)
        )  # eq 15-9, pg 161
        t_c = np.tan(np.pi / 4 - phi_c / 2.0) / np.power(
            ((1 - e * np.sin(phi_c)) / (1 + e * np.sin(phi_c))), (e / 2)
        )  # eq 15-9, pg 161

        # m = np.cos(phi) / np.sqrt(
        #     1 - e * e * np.square(np.sin(phi))
        # )  # eq 14-15, pg 160
        m_c = np.cos(phi_c) / np.sqrt(
            1 - e * e * np.square(np.sin(phi_c))
        )  # eq 14-15, pg 160

        rho = r_e * m_c * t / t_c

        x = rho * np.sin(rlambda - lambda_0)
        y = -rho * np.cos(rlambda - lambda_0)
        y[lats < 0.0] = -y[lats < 0.0]

        # k=rho/(a*m)

        return x, y


def polarstereo_fwd_SP(
    lats, lons, r_e=6378.2730, e=0.081816153, std_parallel=70.0, lon_y=180.0
):
    #
    #   The National Snow and Ice Data Center (NSIDC) and Scientific Committee
    #   on Antarctic Research (SCAR) use a version of the polar stereographic
    #   projection.
    #
    #   Equations from: Map Projections - A Working manual - by J.P. Snyder. 1987
    #   http://kartoweb.itc.nl/geometrics/Publications/Map!20Projections!20-!20A!20Working!20manual!20-!20by!20J.P.!20Snyder.pdf
    #   P 160-161
    #   See the section on Polar Stereographic, with a south polar aspect and
    #   known phi_c not at the pole.
    #
    import numpy as np

    with np.errstate(invalid="ignore"):
        # constants for NSIDC projection
        xphi_c = std_parallel  # standard parallel, latitude of true scale in degrees
        xlambda_0 = lon_y  # meridian along positive Y axis

        # internal variables

        # convert to radians
        phi = np.pi / 180.0 * lats
        phi_c = np.pi / 180.0 * xphi_c * np.ones_like(phi)
        phi_c[phi < 0.0] = -phi_c[phi < 0]

        rlambda = np.pi / 180.0 * lons
        lambda_0 = np.pi / 180.0 * xlambda_0 * np.ones_like(phi)

        #   pm=-1 !plus or minus, north lat. or south
        phi[phi_c < 0.0] = -phi[phi_c < 0.0]
        rlambda[phi_c < 0.0] = -rlambda[phi_c < 0.0]
        lambda_0[phi_c < 0.0] = -lambda_0[phi_c < 0.0]
        phi_c[phi_c < 0.0] = -phi_c[phi_c < 0.0]
        # Equation numbers are from Snyder 1987 "Map Projections a Working Manual"
        # A copy of the paper is in
        # X:\DOCUMENT\Journal_Papers\Snyder_1987_projections_USGS.pdf

        t = np.tan(np.pi / 4.0 - phi / 2.0) / np.power(
            ((1 - e * np.sin(phi)) / (1 + e * np.sin(phi))), (e / 2)
        )  # eq 15-9, pg 161
        t_c = np.tan(np.pi / 4 - phi_c / 2.0) / np.power(
            ((1 - e * np.sin(phi_c)) / (1 + e * np.sin(phi_c))), (e / 2)
        )  # eq 15-9, pg 161

        # m = np.cos(phi) / np.sqrt(
        #     1 - e * e * np.square(np.sin(phi))
        # )  # eq 14-15, pg 160
        m_c = np.cos(phi_c) / np.sqrt(
            1 - e * e * np.square(np.sin(phi_c))
        )  # eq 14-15, pg 160

        rho = r_e * m_c * t / t_c

        x = rho * np.sin(rlambda - lambda_0)
        y = rho * np.cos(rlambda - lambda_0)
        y[lats < 0.0] = -y[lats < 0.0]

        return x, y


def polarstereo_inv(x, y, r_e=6378.2730, e=0.081816153, std_parallel=70.0, lon_y=-45.0):
    #
    #   The National Snow and Ice Data Center (NSIDC) and Scientific Committee
    #   on Antarctic Research (SCAR) use a version of the polar stereographic
    #   projection.

    # Equation numbers are from Snyder 1987 "Map Projections a Working Manual"
    # A copy of the paper is in
    # X:\DOCUMENT\Journal_Papers\Snyder_1987_projections_USGS.pdf
    # import numpy as np

    with np.errstate(invalid="ignore"):
        sgn = 1.0
        e2 = e * e
        e4 = e2 * e2
        e6 = e2 * e4
        e8 = e2 * e6

        sl = std_parallel * np.pi / 180.0

        rho = np.sqrt(x * x + y * y)

        m_c = np.cos(sl) / np.sqrt(
            1.0 - e2 * (np.power(np.sin(sl), 2))
        )  # eq 14-15, pg 160
        t_c = np.tan((np.pi / 4.0) - (sl / (2.0))) / np.power(
            ((1.0 - e * np.sin(sl)) / (1.0 + e * np.sin(sl))), (e / 2.0)
        )  # eq 15-9, pg 161
        if np.abs(std_parallel - 90.0) < 1.0e-5:
            t = (
                rho
                * np.sqrt(
                    (np.power((1.0 + e), (1.0 + e)) * np.power((1.0 - e), (1.0 - e)))
                )
                / (2.0 * r_e)
            )  # eq 21-39, pg 162
        else:
            t = rho * t_c / (r_e * m_c)  # eq 21-40 pg 162

        chi = (np.pi / 2.0) - 2.0 * np.arctan(t)  # eq 7-13, pg 162
        alat = chi + (
            (e2 / 2.0) + (5.0 * e4 / 24.0) + (e6 / 12.0) + (13 * e8 / 360)
        ) * np.sin(
            2.0 * chi
        )  # eq 3-5, pg 162
        alat = alat + (
            (7.0 * e4 / 48.0) + (29.0 * e6 / 240.0) + (811.0 * e8 / 11520)
        ) * np.sin(
            4.0 * chi
        )  # this is an approximation to the true inverse
        alat = alat + ((7.0 * e6 / 120.0) + (81.0 * e8 / 1120)) * np.sin(6.0 * chi)
        # it is accurate to more or less machine precision for
        # locations near the pole.
        alat = alat + (4279.0 * e8 / 161280.0) * np.sin(8.0 * chi)

        alat = sgn * alat
        along = np.arctan2(sgn * x, -sgn * y)
        # along = np.arctan(x/-y)
        along = sgn * along

        lat = alat * 180.0 / np.pi
        lon = along * 180.0 / np.pi + lon_y

        if lat.size > 1:
            lat[rho < 0.1] = 90.0
            lon[rho < 0.1] = 0.0
            lon[lon > 180.0] = lon[lon > 180.0] - 360.0
            lon[lon < -180.0] = lon[lon < -180.0] + 360.0
        else:
            if rho < 0.1:
                lat = 90.0
                lon = 0.0
            if lon > 180.0:
                lon = lon - 360.0
            if lon < -180.0:
                lon = lon + 360.0

        return lon, lat


def ease2(lat, lon):
    import numpy as np

    #     from WGS84
    a = 6378.13700  # Equatorial radius in meters
    q_theta_90 = 1.995531087
    ecc = 0.0818191908426
    #      integer(4) isign
    #      real(8) rho, radicand
    #      real(8) q_theta, sin_lat, e_sin_lat
    isign = np.ones_like(lat)
    isign[lat > 0.0] = -1.0

    sin_lat = np.sin(np.pi / 180.0 * lat)
    e_sin_lat = ecc * sin_lat

    # q_theta = (1 - ecc**2) *
    # (sin_lat / (1 - e_sin_lat**2) - 1 / (2 * ecc) *
    # log((1 - e_sin_lat)/(1 + e_sin_lat)))

    q_theta = (1 - ecc * ecc) * (
        sin_lat / (1 - e_sin_lat * e_sin_lat)
        - 1 / (2 * ecc) * np.log((1 - e_sin_lat) / (1 + e_sin_lat))
    )

    radicand = q_theta_90 + isign * q_theta
    rho = a * np.sqrt(radicand)
    rho[np.abs(radicand) < 1.0e-8] = 0.0

    x = rho * np.sin(lon * np.pi / 180.0)
    y = rho * np.cos(lon * np.pi / 180.0)

    return x, y


class NSIDC_ease2_grids:
    def __init__(self, pole="north", resolution="25km"):
        lats, lons, crs = _load_NSIDC_ease2_grids(pole=pole, resolution=resolution)

        self.latitude = lats
        self.longitude = lons
        self.crs = crs
        self.initialized = True

        return


def _load_NSIDC_ease2_grids(pole="north", resolution="25km"):
    import xarray as xr

    if pole == "north":
        if resolution == "25km":
            filename = polar_grid_root / "NSIDC0772_LatLon_EASE2_N25km_v1.0.nc"
        elif resolution == "12.5km":
            filename = polar_grid_root / "NSIDC0772_LatLon_EASE2_N12.5km_v1.0.nc"
        else:
            raise ValueError(
                f"resolution = {resolution} "
                "not valid in polar_grids.load_NSIDC_ease2_grids"
            )
    elif pole == "south":
        if resolution == "25km":
            filename = polar_grid_root / "NSIDC0772_LatLon_EASE2_S25km_v1.0.nc"
        elif resolution == "12.5km":
            filename = polar_grid_root / "NSIDC0772_LatLon_EASE2_S12.5km_v1.0.nc"
        else:
            raise ValueError(
                f"resolution = {resolution} "
                "not valid in polar_grids.load_NSIDC_ease2_grids"
            )
    else:
        raise ValueError(
            f"pole = {pole} not valid in polar_grids.load_NSIDC_ease2_grids"
        )

    ds = xr.open_dataset(filename)
    lats = ds["latitude"].values
    lons = ds["longitude"].values
    crs = ds["crs"].attrs
    return lats, lons, crs


def ease1_fwd(lat, lon, r=6371.22800, lambda0=0.0, south_pole=False):
    import numpy as np

    #     from WGS84
    r = 6371.22800  # Equatorial radius in kilometers

    if south_pole:
        x = (
            2.0
            * r
            * np.sin((lon - lambda0) * np.pi / 180.0)
            * np.sin(np.pi / 4 - 0.5 * lat * np.pi / 180.0)
        )  # eq 24-8, pg 186
        y = (
            2.0
            * r
            * np.cos((lon - lambda0) * np.pi / 180.0)
            * np.sin(np.pi / 4 - 0.5 * lat * np.pi / 180.0)
        )  # eq 24-9, pg 186
    else:
        x = (
            2.0
            * r
            * np.sin((lon - lambda0) * np.pi / 180.0)
            * np.sin(np.pi / 4 - 0.5 * lat * np.pi / 180.0)
        )  # eq 24-3, pg 186
        y = (
            -2.0
            * r
            * np.cos((lon - lambda0) * np.pi / 180.0)
            * np.sin(np.pi / 4 - 0.5 * lat * np.pi / 180.0)
        )  # eq 24-4, pg 186

    return x, y


def ease1_inv(x, y, r=6371.22800, lambda0=0.0, south_pole=False):
    import numpy as np

    rho = np.sqrt(x * x + y * y)

    c = 2.0 * np.arcsin(rho / (2.0 * r))

    if south_pole:
        phi = np.arcsin(-np.cos(c))  # eq 20-14, pg 186, assuming phi_1 = 90 degrees
        lambda_out = lambda0 * np.pi / 180.0 + np.arctan2(x, y)
    else:  # north pole
        phi = np.arcsin(np.cos(c))  # eq 20-14, pg 186, assuming phi_1 = 90 degrees
        lambda_out = lambda0 * np.pi / 180.0 + np.arctan2(x, -y)
    #    if south_pole:
    #        phi[rho < 0.00001] = -np.pi/2.0
    #        lambda_out[rho < 0.00001] = lambda0
    #    else:  # north pole
    #        phi[rho < 0.00001] = np.pi/2.0
    #        lambda_out[rho < 0.00001] = lambda0

    lat = phi * 180.0 / np.pi
    lon = lambda_out * 180.0 / np.pi

    return lon, lat


def modis_tile_inv(
    vertical_tile,
    horizontal_tile,
    line,
    sample,
    r=6371.22800,
    lambda0=0.0,
    south_pole=False,
):
    # Computes the latitude and longitude of grid points in the MODIS 1km
    # sea-ice tiles from NSIDC

    # The code was tested against results from the on-line calculator located at
    # https://landweb.modaps.eosdis.nasa.gov/cgi-bin/developer/tilemap.cgi
    # the scale factor (1.0027010) was found here:
    # https://nsidc.org/data/MYD29P1D/versions/6 , Table 3

    # currently only works for the North pole projection....

    x = (sample - 475) * 1.0027010
    y = (475 - line) * 1.0027010

    x = x + (horizontal_tile - 9) * 951.0 * 1.0027010
    y = y + (9 - vertical_tile) * 951.0 * 1.0027010

    lon, lat = ease1_inv(x, y, south_pole=False)

    return lon, lat


# import numpy as np

if __name__ == "__main__":
    import numpy as np

    test_polar_stereo_interp = True
    test_polar_projection = True

    if test_polar_stereo_interp:
        xlocs = np.arange(-154, 150, dtype=np.float64)
        xlocs = 12.5 + 25.0 * xlocs
        ylocs = np.arange(-214, 234, dtype=np.float64)
        ylocs = 12.5 + 25.0 * ylocs
        xgrid, ygrid = np.meshgrid(xlocs, ylocs)

        lon_grid, lat_grid = polarstereo_inv(xgrid, ygrid)

        test_lons = np.array([-45.1, 135.2, 45.2, 0.12, 123.4, 120.0])
        test_lats = np.array([80.0, 80.0, 80.0, 70.0, 67.2, 89.5])

        test_lon_result = polar_stereo_interp(lon_grid, test_lats, test_lons)
        test_lat_result = polar_stereo_interp(lat_grid, test_lats, test_lons)
        np.set_printoptions(precision=8)
        print("------- lons")
        print(["{:10.5f}".format(x) for x in test_lons])
        print(["{:10.5f}".format(x) for x in test_lon_result])
        print(["{:10.5f}".format(x) for x in test_lon_result - test_lons])
        print("------- lats")
        print(["{:10.5f}".format(x) for x in test_lats])
        print(["{:10.3f}".format(x) for x in test_lat_result])
        print(["{:10.5f}".format(x) for x in test_lat_result - test_lats])

        print()

    if test_polar_projection:
        """
            h_tiles = np.array([9,10,12,4],dtype='Int32')
            v_tiles = np.array([9,9,11,8],dtype='Int32')

            lines = np.array([0,23,123,45],dtype='Int32')
            samples= np.array([0,600,566,3],dtype='Int32')

            test_lons = np.array([-135.0,14.063257,38.550684,-164.653289])
            test_lats = np.array([83.939869,85.360129,60.849892,39.939949])

            lons,lats = modis_tile_inv(h_tiles,v_tiles,lines,samples)

            print('Lons  ',("{:10.3f}".format(lons)))
            print('Test  ',("{:10.3f}".format(test_lons)))
            print('Lats  ',("{:10.3f}".format(lats)))
            print('Test  ',("{:10.3f}".format(test_lats)))

            print('-------------------------------------------')

            test=np.array([-890000.0, -629000.0,  79.9641229, -99.7495626,
            1720000.0, -629000.0,  73.2101233,  24.9126514,
                          -890000.0, -3410000.0, 58.2706251, -59.6277136,
            20000.0, -10000.0,     89.7592932, -18.2336765])

            test = np.reshape(test,((4,4)))
            lats  = test[:,2]
            lons  = test[:,3]

            np.set_printoptions(precision = 8)
            x_ease1,y_ease1 = ease1_fwd(lats,lons,lambda0 = -45.0)
            print('X,EASE1 ', x_ease1)
            print('Y,EASE1 ', y_ease1)
            print('RHO,EASE1 ',np.sqrt(x_ease1*x_ease1 + y_ease1*y_ease1))

            lons2,lats2 = ease1_inv(x_ease1,y_ease1,lambda0 = -45.0)

            print('Lats inv',lats2)
            print('Lats orig',lats)
            print('Lons Inv',lons2)
            print('Lons orig',lons)
            print('----------------------------------------------------------------')
            x_stereo,y_stereo = polarstereo_fwd(lats,lons)

            print('Lats',lats)
            print('Lons',lons)
            print('X,polar_stereo ', x_stereo)
            print('Y,polar_stereo ', y_stereo)
            print('RHO,Pol_St ',np.sqrt(x_stereo*x_stereo + y_stereo*y_stereo))

            lons2,lats2 = polarstereo_inv(x_stereo,y_stereo)

            print('Lats inv',lats2)
            print('Lats orig',lats)
            print('Lons Inv',lons2)
            print('Lons orig',lons)


            lons = lons + 45.
            x_ease1,y_ease1 = ease1_fwd(lats,lons)
            x_ease1 = x_ease1/1000.0
            y_ease1 = y_ease1/1000.0

            x_ease2,y_ease2 = ease2(lats,lons)
            x_ease2 = x_ease2/1000.0
            y_ease2 = y_ease2/1000.0
            np.set_printoptions(precision=3)
            print('X,polar_stereo ', x_stereo)
            print('X,EASE         ', x_ease1)
            print('X,EASE2        ', x_ease2)
            print()
            print('Y,polar_stereo ', y_stereo)
            print('Y,EASE         ', y_ease1)
            print('Y,EASE2        ', y_ease2)
            print()

        #

            #print('check the code using Snyders example.
            #Should get x=-1540033.6; y=-560526.4')
            #phi=-75; lambda=150;
            #x,y=polarstereo_fwd(phi,lambda,6378388.0,0.0819919,-71,-100)
            #print(x,y)
            ##
            ##!!!!!!!!!!!
            ##!check with AntDEM
            ##!!!!!!!!!!!
            ##!http://nsidc.org/data/docs/daac/nsidc0304_0305_glas_dems.gd.html
            ##!Center Point of Corner Grid Cell
            ##!x	y	Latitude	Longitude
            ##! test=[-2812000.0  2299500.0   -57.3452815 -50.7255753
            ##    ! 2863500.0   2299500.0   -57.0043684 51.2332036
            ##    ! -2812000.0  -2384000.0  -56.8847122 -130.2911169
            ##    ! 2863500.0   -2384000.0  -56.5495152  129.7789915];
            ##! [x,y]=polarstereo_fwd(test(:,3),test(:,4),6378137.0,
            ##! axes2ecc(6378137.0, 6356752.3),-70,0);
            ##! figure,hold on,plot(test(:,1),test(:,2),'.'),plot(x,y,'r+')
            ##! [test(:,1) test(:,1)-x],[test(:,2) test(:,2)-y]
            ##!error is less than half a meter (probably just round-off error).
            ##
            ##!!!!!!!!!!!
            ##!check with Greenland
            ##!!!!!!!!!!!
            ##!projected from the WGS 84 Ellipsoid, with 70 N as the latitude
            ##!of true scale and a rotation of 45.
            ##! test=[-890000.0 -629000.0 79.9641229 -99.7495626
            ##!center point of cell
            ##    ! 1720000.0 -629000.0 73.2101233 24.9126514
            ##    ! -890000.0 -3410000.0 58.2706251 -59.6277136
            ##    ! 1720000.0 -3410000.0 55.7592932 -18.2336765];
            ##! [x,y]=polarstereo_fwd(test(:,3),test(:,4),
            ##! 6378273,0.081816153,70,-45);
            ##!slightly off
            ##![x2,y2]=polarstereo_fwd(test(:,3),test(:,4),
            # #! 6378137.0,0.08181919,70,-45);
            ##!correct
            ##! figure,hold on,plot(test(:,1),test(:,2),'.'),plot(x,y,'r+'),
            ##! plot(x2,y2,'gx')
            ##! [test(:,1) test(:,1)-x test(:,1)-x2],[test(:,2)
            ##! test(:,2)-y test(:,2)-y2]
            ##!error is less than half a meter (probably just round-off error).
            ##!}
        """
        # north pole corner test

        test_lats = np.array([30.98, 39.43, 31.37, 56.35, 34.35, 43.28, 33.92, 55.50])
        test_lons = np.array(
            [168.35, 135.00, 102.34, 45.00, 350.03, 315.00, 279.26, 225.0]
        )
        test_x = np.array([-3850.0, 0.0, 3750.0, 3750.0, 3750.0, 0.0, -3850.0, -3850.0])
        test_y = np.array([5850.0, 5850.0, 5850.0, 0.0, -5350.0, -5350.0, -5350.0, 0.0])

        x, y = polarstereo_fwd(
            test_lats, test_lons, r_e=6378.2730, e=0.081816153, std_parallel=70.0
        )
        print("North Pole")
        for i, lat in enumerate(test_lats):
            print(
                (
                    "{:10.3f}".format(test_lats[i]),
                    "{:10.3f}".format(test_lons[i]),
                    "{:10.3f}".format(x[i]),
                    "{:10.3f}".format(y[i]),
                    "{:10.3f}".format(x[i] - test_x[i]),
                    "{:10.3f}".format(y[i] - test_y[i]),
                )
            )

        # south pole corner test
        test_lons = np.array([317.76, 0.0, 42.24, 90.0, 135.0, 180.0, 225.0, 270.0])
        test_lats = np.array(
            [-39.23, -51.32, -39.23, -54.66, -41.45, -54.66, -41.45, -54.66]
        )
        test_x = np.array([-3950.0, 0.0, 3950.0, 3950.0, 3950.0, 0.0, -3950.0, -3950.0])
        test_y = np.array([4350.0, 4350.0, 4350.0, 0.0, -3950.0, -3950.0, -3950.0, 0.0])

        x, y = polarstereo_fwd_SP(
            test_lats,
            test_lons,
            r_e=6378.2730,
            e=0.081816153,
            std_parallel=70.0,
            lon_y=180.0,
        )
        print("South Pole")
        for i, lat in enumerate(test_lats):
            print(
                (
                    "{:10.3f}".format(test_lats[i]),
                    "{:10.3f}".format(test_lons[i]),
                    "{:10.3f}".format(x[i]),
                    "{:10.3f}".format(y[i]),
                    "{:10.3f}".format(x[i] - test_x[i]),
                    "{:10.3f}".format(y[i] - test_y[i]),
                )
            )
