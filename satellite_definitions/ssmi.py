import numpy as np

SAT_NAME = "SSMI"

# SSMI Frequencies
REF_FREQ = np.array([19.35,22.235,37.0], np.float32)

# Mapping to footprint sizes
REF_FREQ_mapping = np.array([0, 1, 2], np.int32)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
REF_EIA = np.array([53.1, 53.1, 53.1, 53.1, 53.1, 53.1], np.float32)

CHANNEL_TO_FREQ_MAP = np.array([0, 0, 1, 1, 2, 2],np.int32)
CHANNEL_TO_POL_MAP =  np.array([0, 1, 0, 1, 0, 1],np.int32) # 0 = V, 1 = H
