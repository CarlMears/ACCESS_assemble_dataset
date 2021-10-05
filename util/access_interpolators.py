def time_interpolate_synoptic_maps_ACCESS(map_array,map_times,time_map):

    import numpy as np
    import datetime

    '''this routine interpolates an array of gridded data (map_array), assumed cover an entrire day, including both 0:00 and 24:00
     time slots.  The map_times array is a 1-d array of the times for each map in seconds since midnight'''

    sz = map_array.shape
    num_maps = sz[0]
    num_maps2 = len(map_times)

    if (num_maps != num_maps2):
        raise ValueError('# of maps in map array not equal to number of map times')

    time_step = np.rint(map_times[1] - map_times[0])

    num_steps_in_day = len(map_times)-1
    hr_interp = np.copy(time_map/60.0)           # convert to hours after midnight ZZ without side effects

    map_interp = np.zeros_like(map_array[0,:,:])*np.nan     # output array -- initialize to be all NAN
    
    for interval in range(0,num_steps_in_day):
        ok = np.all([(time_map >= map_times[interval]),(time_map < map_times[interval+1])],axis = (0))
        num_ok = np.sum(ok)
        if num_ok > 0:
            y = np.zeros((num_ok,2))
            y[:,0] = map_array[interval,:,:][ok]
            y[:,1] = map_array[interval+1,:,:][ok]

            wt_lower = (map_times[interval+1]-time_map[ok])/time_step

            assert(np.all(wt_lower >= 0.0)),"wt_lower < 0.0"  # check to make sure weights are in bounds
            assert(np.all(wt_lower <= 1.0)),"wt_lower > 1.0"

            wt_upper = 1.0 - wt_lower

            y_interp = wt_lower*y[:,0] + wt_upper*y[:,1]
            map_interp[ok] = y_interp


    return map_interp #,test_map