import numpy as np

land_fraction_path = '//ops1p-to-be-renamed/M/Module_Data/Land_Fraction/'

def read_land_fraction_1440_720(path = land_fraction_path):

    bin_file = land_fraction_path + 'land_fraction_1440x720.dat'

    with open(bin_file, mode='rb') as file:
        land_frac_raw = np.fromfile(file,dtype = np.float32)

    land_frac_raw.shape
    land_frac = np.reshape(land_frac_raw,(720,1440))
    return land_frac

def read_land_fraction_polar_stereographic(pole='north'):
    if pole == 'north':
        return read_land_fraction_polar_stereographic_NP()
    if pole == 'south':
        return read_land_fraction_polar_stereographic_SP()

    raise ValueError('arg "pole" should be "north" or "south"')
    
def read_land_fraction_polar_stereographic_NP():

    import numpy as np
    
    land_fraction_file = 'l:/sea_ice/land_fraction/nsidc_polar_stereographic_land_fraction_north_pole.dat'
    
    lf = np.fromfile(land_fraction_file,dtype = 'float64')
    lf = np.reshape(lf,(448,304))
    
    return lf


def read_land_fraction_polar_stereographic_SP():
    import numpy as np

    land_fraction_file = 'l:/sea_ice/land_fraction/nsidc_polar_stereographic_land_fraction_south_pole.dat'

    lf = np.fromfile(land_fraction_file, dtype='float64')
    lf = np.reshape(lf, (332,316))
    return lf