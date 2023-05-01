import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt

from rss_plotting.plot_2d_array import plot_2d_array

fwhm_to_alpha = 1.0 / (2 * np.sqrt(2 * np.log(2)))


def find_sigma_from_FWHM(fwhm):
    return fwhm * fwhm_to_alpha


def integral_of_2D_gaussian(fwhm, x1, x2, y1, y2):

    sigma = find_sigma_from_FWHM(fwhm)
    sigma_sqrt2 = find_sigma_from_FWHM(fwhm) * np.sqrt(2)

    # print(np.exp(-(x2*x2)/(2.0*sigma*sigma)))

    x_part = 0.5 * (erf(x2 / sigma_sqrt2) - erf(x1 / sigma_sqrt2))
    y_part = 0.5 * (erf(y2 / sigma_sqrt2) - erf(y1 / sigma_sqrt2))

    return x_part * y_part


def find_weights_1440_721(*, lat, lon, fwhm, compact=True):
    """finds the weights for a 1/4-degree edge-centered grid for a arbitrary latitude/longitude location
    Input:
         lat: latitude
         lon: longitude
         fwhm: FWHM of the target gaussian
     Return:
         721x1440 array of weights
         if compact is true, this is in a form of a dictionary of 3 arrays
             rows
             cols
             wts
         if compact is False
             the full 721x1440 array
    """
    grid_delta = 0.25
    grid_offset_lat = -90.0
    grid_offset_lon = 0.0
    ilat_ll = np.floor((lat - grid_offset_lat) / grid_delta)
    ilon_ll = np.floor((lon - grid_offset_lon) / grid_delta)


if __name__ == "__main__":

    fwhm = 30.0
    lat = 60.0

    wts_all = np.zeros((15, 15, 721), dtype=np.float32)

    for ilat, lat in enumerate(np.arange(-90.0, 90.1, 0.25)):

        delta_y = 0.25 * 111
        delta_x = 0.25 * 111 * np.cos(lat * np.pi / 180.0)

        wts = np.zeros((15, 15), dtype=np.float32)

        for ix in np.arange(-7, 8):
            for iy in np.arange(-7, 8):
                x1 = (ix - 0.5) * delta_x
                x2 = x1 + delta_x
                y1 = (iy - 0.5) * delta_y
                y2 = y1 + delta_y
                wts[iy + 7, ix + 7] = integral_of_2D_gaussian(fwhm, x1, x2, y1, y2)

        print(ilat, lat, np.sum(wts))
        wts = wts / np.sum(wts)
        wts_all[:, :, ilat] = wts

        if ilat % 40 == 0:
            xvals = 0.25 * np.arange(-7, 8)
            yvals = 0.25 * np.arange(-7, 8)
            fig, ax = plot_2d_array(
                wts, xvals, yvals, zrange=(0.0, 1.2 * np.max(wts)), cmap="viridis"
            )
            plt.show()

    print
