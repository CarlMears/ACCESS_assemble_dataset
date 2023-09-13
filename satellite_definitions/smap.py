import numpy as np

SAT_NAME = "SMAP"

# SMAP Frequencies
REF_FREQ = np.array([1.4], np.float32)

# Mapping to footprint sizes
REF_FREQ_mapping = np.array([0], np.int32)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
REF_EIA = np.array([40.22], np.float32)

CHANNEL_TO_FREQ_MAP = np.array([0, 0, 0, 0],np.int32)
CHANNEL_TO_POL_MAP =  np.array([0, 1, 2, 3],np.int32) # 0 = V, 1 = H, 2 = 3rd Stokes, 3 = 4th Stokes
