import numpy as np

SAT_NAME = 'AMSR2'

#AMSR2 Frequencies
REF_FREQ = np.array([6.9, 7.3, 10.7, 18.7, 23.8, 37.0], np.float32)

#Mapping to footprint sizes
REF_FREQ_mapping = np.array([1,1,2,3,4,5],np.int32)

# Reference Earth incidence angle to use for each reference frequency (in degrees)
REF_EIA = np.array([53.0, 53.0, 53.0, 53.0, 53.0, 53.0], np.float32)

CHANNEL_TO_FREQ_MAP = np.array([0,0,1,1,2,2,3,3,4,4,5,5,6,6])
CHANNEL_TO_POL_MAP =  np.array([0,1,0,1,0,1,0,1,0,1,0,1,0,1])