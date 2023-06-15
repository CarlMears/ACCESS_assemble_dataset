import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
from pathlib import Path
import concurrent.futures

"""
To use the python-wrapped fortran version of the resampler,
the f2py-based resamp_using_wts must be installed

To do this:
    clone the project from http://gitlab.remss.com/access/
    resample_using_weights_fortran_python
    (or cd to http://gitlab.remss.com/access/resample_using_weights_fortran_python)
    issue the command

    python -m build

    Ideally, this will build the project with no errors.  There are a LOT of warnings.

    then

    python -m pip install dist/samp-0.0.1-cp310-cp310-linux_x86_64.whl
    (or whatever wheel was generated for your flavor of python)

"""


def resample_iy_row(xindex, yindex, weights, var, iy):
    if iy % 10 == 0:
        print(iy)
    numx = xindex.shape[0]
    num_time_steps = var.shape[0]

    var_resamp = np.full((num_time_steps, numx), np.nan, dtype=np.float32)
    for ix in np.arange(numx):
        wts = weights[ix, :]
        x_indices = xindex[ix, :]
        y_indices = yindex[ix, :]
        if np.min(x_indices) < 0:
            continue
        if np.min(y_indices) < 0:
            continue
        for itime in range(num_time_steps):
            var_sub = var[itime, :][y_indices, x_indices]
            temp = np.sum(var_sub * wts)
            var_resamp[itime, ix] = temp
    return var_resamp, iy


class ResampleERA5:
    def __init__(self, *, target_size, region):
        if os.name == "nt":
            resample_wt_path = Path("L:/access/era5/resample_weights")
        elif os.name == "posix":
            resample_wt_path = Path("/mnt/ops1p-ren/l/access/era5/resample_weights")
        if region in ["north", "south"]:
            grid_type = "ease2"
            if region == "north":
                nc_file = (
                    resample_wt_path
                    / f"resamp_wts_ease2_25km_NP_{int(target_size)}km.nc"
                )
            else:
                nc_file = (
                    resample_wt_path
                    / f"resamp_wts_ease2_25km_SP_{int(target_size)}km.nc"
                )
        elif region == "global":
            grid_type = "equirectangular"
            nc_file = (
                resample_wt_path / f"resamp_wts_1440)721)rect_{int(target_size)}km.nc"
            )

        print(f"Initializing resampler for {target_size} targets for {region}")
        print(f"Loading weights and indices from {nc_file}")

        ds = xr.open_dataset(nc_file)
        self.weights = ds["weights"].values
        self.yindex = ds["lat_index"].values
        self.xindex = ds["lon_index"].values

        # before passing to the FORTRAN routine make fortran versions of the arrays
        # 1.  Swap axes so it is in column major order with weight dimension first
        # 2.  Ensure that it stored in this order using "asfortranarray"
        # 3.  Add one to the indices to agree with fortran default

        self.weightsF = np.asfortranarray(np.swapaxes(ds["weights"].values, 0, 2))
        self.yindexF = np.asfortranarray(np.swapaxes(ds["lat_index"].values, 0, 2) + 1)

        self.xindexF = np.asfortranarray(np.swapaxes(ds["lon_index"].values, 0, 2) + 1)
        # self.xindexF[self.xindexF > 1440] = self.xindexF[self.xindexF > 1440] - 1440

        try:
            self.lats = ds["latitude"].values
        except KeyError:
            self.lats = ds["Latitude"].values
        try:
            self.lats = ds["longitude"].values
        except KeyError:
            self.lats = ds["Longitude"].values
        self.target_size = target_size
        self.grid_type = grid_type
        self.region = region
        self.numy = self.weights.shape[0]
        self.numx = self.weights.shape[1]

    def resample(self, var):
        """This python only version is much, much slower than the fortran version,
        even though it uses multiple processes"""

        max_workers = 10

        num_time_steps, _, _ = var.shape

        var_resamp = np.full(
            (num_time_steps, self.numy, self.numx), np.nan, dtype=np.float32
        )

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers
        ) as executor:
            results = {
                executor.submit(
                    resample_iy_row,
                    self.xindex[iy, :, :],
                    self.yindex[iy, :, :],
                    self.weights[iy, :, :],
                    var,
                    iy,
                )
                for iy in range(0, self.numy)
                # for iy in range(0, 60)
            }
        for future in concurrent.futures.as_completed(results):
            try:
                (row, iy) = future.result()
                var_resamp[:, iy, :] = row
            except KeyboardInterrupt:
                return
            except Exception as e:
                print(f"Error in run: {e}")
        return var_resamp

    def resample_fortran(self, var, verbose=False):
        # To import the following, the module needs to be installed using build/f2py
        #
        # see instruction in the README.md in
        #
        # gitlab.remss.com/access/resample_using_weights_fortran_python

        from resamp import resamp_using_wts

        if verbose:
            iverbose = 1
        else:
            iverbose = 0

        var_shape = var.shape
        if len(var_shape) == 3:
            num_time_steps, _, _ = var.shape
            var_resamp = np.full(
                (num_time_steps, self.numy, self.numx), np.nan, dtype=np.float32
            )

            for itime in range(num_time_steps):
                var_resamp_step = resamp_using_wts(
                    self.xindexF,
                    self.yindexF,
                    self.weightsF,
                    np.transpose(var[itime, :, :]),
                    iverbose,
                )
                var_resamp[itime, :, :] = var_resamp_step
        else:
            var_resamp = resamp_using_wts(
                self.xindexF, self.yindexF, self.weightsF, np.transpose(var), iverbose
            )

        if verbose:
            print("Finished Resampling in Fortran")
        return var_resamp


if __name__ == "__main__":
    from rss_plotting.plot_2d_array import plot_2d_array

    resampler = ResampleERA5(target_size=70, region="north")

    plot_2d_array(
        resampler.xindexF[0, :, :],
        xvals=np.arange(0, 720),
        yvals=np.arange(0, 720),
        zrange=[-1440, 1440],
        plt_colorbar=True,
    )
    plot_2d_array(
        resampler.yindexF[0, :, :],
        xvals=np.arange(0, 720),
        yvals=np.arange(0, 720),
        zrange=[-720, 720],
        plt_colorbar=True,
    )

    plt.show()

    var_in = np.zeros((721, 1440), dtype=np.float32)

    for ilat in range(0, 721):
        var_in[ilat, :] = -90.0 + ilat * 0.25

    var_out = np.zeros((720, 720), dtype=np.float32)
    var_out = resampler.resample_fortran(var_in, verbose=False)

    print(f"Max = {np.nanmax(var_out)}")
    print(f"Min = {np.nanmin(var_out)}")

    plot_2d_array(
        var_out, xvals=np.arange(0, 720), yvals=np.arange(0, 720), zrange=[-90.0, 90.0]
    )
    plt.show()

    print
