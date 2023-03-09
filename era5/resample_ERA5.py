import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
from pathlib import Path
import concurrent.futures

def init_worker() -> None:
    """
    This is a function which will allow the user to stop all workers
    with a Ctrl-C event.  Otherwise the code gets caught in an odd loop in
    which it does not terminate.  This appears to be a Python bug:
    https://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def resample_iy_row(xindex,yindex,weights,var,iy):

    if iy%10 == 0:
        print(iy)
    numx = xindex.shape[0]
    num_time_steps = var.shape[0]

    var_resamp = np.full((num_time_steps,numx),np.nan,dtype=np.float32)
    for ix in np.arange(numx):
            wts = weights[ix,:]
            x_indices = xindex[ix,:]
            y_indices = yindex[ix,:]
            if np.min(x_indices) < 0:
                continue
            if np.min(y_indices) < 0:
                continue
            for itime in range(num_time_steps):
                var_sub = var[itime,:][y_indices,x_indices]
                temp = np.sum(var_sub*wts)
                var_resamp[itime,ix] = temp
    return var_resamp,iy

class ResampleERA5:

    def __init__(self,*,target_size,region):

        if os.name == "nt":
            resample_wt_path = Path("L:/access/era5/resample_weights")
        elif os.name == "posix":
            resample_wt_path = Path("/mnt/ops1p-ren/l/access/era5/resample_weights")
        if region in ['north','south']:
            grid_type = 'ease2'
            if region == 'north':
                nc_file = resample_wt_path / f'resamp_wts_ease2_25km_NP_{int(target_size)}km.nc'
            else:
                nc_file = resample_wt_path / f'resamp_wts_ease2_25km_NP_{int(target_size)}km.nc'
        elif region == "global":
            grid_type = 'equirectangular'
            nc_file = resample_wt_path / f'resamp_wts_1440)721)rect_{int(target_size)}km.nc'

        ds = xr.open_dataset(nc_file)
        self.weights = ds['weights'].values
        self.yindex = ds['lat_index'].values
        self.xindex = ds['lon_index'].values

        self.weightsF = np.swapaxes(ds['weights'].values,0,2)
        self.yindexF = np.swapaxes(ds['lat_index'].values,0,2)
        self.xindexF = np.swapaxes(ds['lon_index'].values,0,2)

        try:
            self.lats = ds['latitude'].values
        except KeyError:
            self.lats = ds['Latitude'].values
        try:
            self.lats = ds['longitude'].values
        except KeyError:
            self.lats = ds['Longitude'].values
        self.target_size = target_size
        self.grid_type = grid_type
        self.region = region
        self.numy = self.weights.shape[0]
        self.numx = self.weights.shape[1]
        


    def resample(self,var):

        max_workers = 10

        num_time_steps,num_var_y,num_var_x = var.shape

        var_resamp = np.full((num_time_steps,self.numy,self.numx),np.nan,dtype=np.float32)

        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = {
            executor.submit(
                resample_iy_row,
                self.xindex[iy,:,:],
                self.yindex[iy,:,:],
                self.weights[iy,:,:],
                var,
                iy
            )
            for iy in range(0, self.numy)
            #for iy in range(0, 60)
        }
        for future in concurrent.futures.as_completed(results):
            try:
                (row,iy) = future.result()
                var_resamp[:,iy,:] = row
            except KeyboardInterrupt:
                return
            except Exception as e:
                print(f"Error in run: {e}")
        return var_resamp

    def resample_fortran(self,var):

        # To import the following, the module needs to be installed using build/f2py
        #
        # see instruction in the README.md in
        #
        # gitlab.remss.com/access/resample_using_weights_fortran_python
    
        from resamp import resamp_using_wts

        num_time_steps,_,_ = var.shape
        var_resamp = np.full((num_time_steps,self.numy,self.numx),np.nan,dtype=np.float32)

        for itime in range(num_time_steps):
            var_resamp_step = resamp_using_wts(np.asfortranarray(self.xindexF),
                                               np.asfortranarray(self.yindexF),
                                               np.asfortranarray(self.weightsF),
                                               np.asfortranarray(np.transpose(var[itime,:,:])))
            var_resamp[itime,:,:] = var_resamp_step

        return var_resamp

        

                

        