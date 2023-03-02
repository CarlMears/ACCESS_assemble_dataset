# script to make a set of images for a specific data,hour and location
# images are 960 x360 to conform to AGU poster slideshow parameters


import datetime
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
import os

from rss_plotting.global_map.regional_map import plot_regional_subset_map as plt_region


day_to_do = datetime.date(2017, 7, 12)
hour_of_day = 7
footprint_size = 30
channels = [2, 3, 4, 5]
freqs = [10.65, 18.7, 23.8, 36.5]
pols = ["V", "H"]

lat_center = 45.0
lon_center = 280.0
lat_range = 5.0
lon_range = 12.5

path_to_ACCESS_data = Path("L:/access/amsr2_out_test3")

path_to_this_day = (
    path_to_ACCESS_data
    / f"Y{day_to_do.year}"
    / f"M{day_to_do.month:02d}"
    / f"D{day_to_do.day:02d}"
)

atm_par_file = (
    path_to_this_day
    / f"amsr2_atm_par_era5_{day_to_do.year}_{day_to_do.month:02d}"
      f"_{day_to_do.day:02d}.{footprint_size:03d}km.nc"
)

print(atm_par_file)

ds = xr.open_dataset(atm_par_file)

vars = ["transmissivity", "upwelling_tb", "downwelling_tb"]
varnames = ["Transmissivity", "Tb_up", "Tb_down"]
varmin = [0.5, 0.0, 0.0]
varmax = [1.0, 150.0, 150.0]
varunits = ["", "K", "K"]
freqs = [10.65, 18.7, 23.8, 36.5]

for ivar, var in enumerate(vars):
    var_maps = ds[var]
    for ifreq, freq in enumerate(freqs):
        glb_map = var_maps[:, :, hour_of_day, ifreq + 2].values
        glb_map[glb_map < -50.0] = np.nan

        if len(varunits[ivar]) > 0:
            title_str = f"Calculated AMSR2 (from ERA5 profiles) "\
                        f"{varnames[ivar]} {freq:.2f} ({varunits[ivar]})"
        else:
            title_str = (
                f"Calculated AMSR2 (from ERA5 profiles) {varnames[ivar]} {freq:.2f}"
            )

        fig = plt.figure(figsize=(9.6, 3.6))
        fig, ax = plt_region(
            glb_map,
            panel_label_loc=[0.03, 0.90],
            fig_in=fig,
            vmin=varmin[ivar],
            vmax=varmax[ivar],
            cmap="viridis",
            title=None,
            central_longitude=lon_center,
            central_latitude=lat_center,
            longitude_size=2.0 * lon_range,
            latitude_size=2.0 * lat_range,
            units=title_str,
            return_map=False,
            panel_label=None,
            plt_colorbar=True,
        )
        date_str = f"{day_to_do.year:04d}-{day_to_do.month:02d}-{day_to_do.day:02d}"
        png_path = Path(f"M:/job_access/plots/image_sets/{date_str}")
        os.makedirs(png_path, exist_ok=True)
        png_file = png_path / f"atm_par_{var}_{freq:.2f}_{date_str}.png"
        fig.savefig(png_file)
        plt.close(fig)


vars = ["tclw_era5", "tcwv_era5", "rainfall_rate", "u10n_era5", "land_frac_modis"]
varmin = [0.0, 0.0, 0.0, -15.0, 0.0]
varmax = [1.0, 60.0, 10.0, 15.0, 1.0]
var_units = [
    "$\mathregular{kg/m^{2}}$",
    "$\mathregular{kg/m^{2}}$",
    "mm/hour",
    "m/s",
    "",
]
var_names = [
    "Total Column Cloud Water from ERA5",
    "Total Column Water Vapor from ERA5",
    "Rainfall Rate from IMERG",
    "Zonal Wind from ERA5",
    "Land Fraction from MODIS Land Mask",
]
for ivar, var in enumerate(vars):
    file = (
        path_to_this_day
        / f"amsr2_{var}_{day_to_do.year}_{day_to_do.month:02d}_"
          f"{day_to_do.day:02d}.{footprint_size:03d}km.nc"
    )
    print(file)
    ds = xr.open_dataset(file)
    if var == "land_frac_modis":
        glb_var_map = ds["land_fraction"].values
    else:
        glb_var_map = ds[var].values
        glb_var_map = glb_var_map[:, :, hour_of_day]

    glb_var_map[glb_var_map < -50.0] = np.nan

    title_str = f"{var_names[ivar]}, {var_units[ivar]}"
    fig = plt.figure(figsize=(9.6, 3.6))
    fig, ax = plt_region(
        glb_var_map,
        panel_label_loc=[0.03, 0.90],
        fig_in=fig,
        vmin=varmin[ivar],
        vmax=varmax[ivar],
        cmap="viridis",
        title=None,
        central_longitude=lon_center,
        central_latitude=lat_center,
        longitude_size=2.0 * lon_range,
        latitude_size=2.0 * lat_range,
        units=title_str,
        return_map=False,
        panel_label=None,
        plt_colorbar=True,
    )

    date_str = f"{day_to_do.year:04d}-{day_to_do.month:02d}-{day_to_do.day:02d}"
    png_path = Path(f"M:/job_access/plots/image_sets/{date_str}")
    os.makedirs(png_path, exist_ok=True)
    png_file = png_path / f"{var}_{date_str}.png"
    fig.savefig(png_file)
    plt.close(fig)

print


resamp_tb_file = (
    path_to_this_day
    / f"amsr2_resamp_tbs_{day_to_do.year}_{day_to_do.month:02d}_"
      f"{day_to_do.day:02d}.{footprint_size:03d}km.nc"
)

print(resamp_tb_file)

ds = xr.open_dataset(resamp_tb_file)

tbs = ds["brightness_temperature"]
lats = ds["latitude"]
lons = ds["longitude"]
ds.keys()

for ifreq, freq in enumerate(freqs):
    for ipol, pol in enumerate(pols):
        glb_tb_map = tbs[:, :, hour_of_day, ifreq + 2, ipol].values
        glb_tb_map[glb_tb_map < 50.0] = np.nan

        title_str = f"AMSR2 Radiance {freq:.2f}{pol} (K)"
        fig = plt.figure(figsize=(9.6, 3.6))
        fig, ax = plt_region(
            glb_tb_map,
            panel_label_loc=[0.03, 0.90],
            fig_in=fig,
            vmin=100.0,
            vmax=300.0,
            cmap="viridis",
            title=None,
            central_longitude=lon_center,
            central_latitude=lat_center,
            longitude_size=2.0 * lon_range,
            latitude_size=2.0 * lat_range,
            units=title_str,
            return_map=False,
            panel_label=None,
            plt_colorbar=True,
        )
        date_str = f"{day_to_do.year:04d}-{day_to_do.month:02d}-{day_to_do.day:02d}"
        png_path = Path(f"M:/job_access/plots/image_sets/{date_str}")
        os.makedirs(png_path, exist_ok=True)
        png_file = png_path / f"resamp_tb_{freq:.2f}{pol}_{date_str}.png"
        fig.savefig(png_file)
        plt.close(fig)
