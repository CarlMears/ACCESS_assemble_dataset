# polar_grids module

## Functions
polarstereo_fwd(
    lats, lons, r_e=6378.2730, e=0.081816153, std_parallel=70.0, lon_y=-45.0
)
polarstereo_inv(
    x, y, r_e=6378.2730, e=0.081816153, std_parallel=70.0, lon_y=-45.0
)

polar_stereo_interp(polar_map, lats, lons)
polar_stereo_interp_SP(polar_map, lats, lons)


## Required Modules
os

