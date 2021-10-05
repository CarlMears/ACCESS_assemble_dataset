def get_access_output_filename(*,year,month,day,satellite,dataroot):

    return f'{dataroot}{satellite}_out/{satellite}_resamp_tbs_{year:04d}_{month:02d}_{day:02d}.nc'

def append_var_to_daily_tb_netcdf(*,year,month,day,satellite,var,
                                    var_name,
                                    standard_name=None,
                                    long_name=None,
                                    valid_min=None,
                                    valid_max=None,
                                    units=None,
                                    v_fill = -999.0,
                                    dataroot='L:/access/',
                                    overwrite=False):

    from netCDF4 import Dataset as netcdf_dataset

   
    filename = get_access_output_filename(year=year,month=month,day=day,satellite=satellite,dataroot=dataroot)

    root_grp = netcdf_dataset(filename, 'a', format='NETCDF4')
    try:
        v = root_grp.createVariable(var_name, 'f4', ('latitude', 'longitude','hours',),zlib=True)
    except:
        if overwrite:
            print(f'Warning -- variable {var_name} already exists - overwriting')
            v = root_grp.variables[var_name]
        else:
            raise RuntimeError(f'Error -- Variable {var_name} already exists and overwrite is False, so raising error')

    if standard_name is not None:
        v.standard_name = standard_name
    else:
        v.standard_name = var_name
    
    if long_name is not None:
        v.long_name = long_name

    v.FillValue = v_fill
    # the next few lines are because an class attribute can not start with "_", but that is what CF requires.
    try:
        v.renameAttribute('FillValue', '_FillValue')
    except:
        #print('can not rename FillValue to _FillValue because it exists')
        v.delncattr('_FillValue')
        v.renameAttribute('FillValue', '_FillValue')

    v.missing = v_fill

    if valid_min is not None:
        v.valid_min = valid_min
    if valid_max is not None:
        v.valid_max = valid_max
    if units is not None:
        v.units=units
    v.coordinates = 'latitude longitude hours'

    v[:,:,:] = var
    root_grp.close()

def append_lf_daily_tb_netcdf(*,year,month,day,satellite,land_fraction,dataroot='L:/access/'):


    from netCDF4 import Dataset as netcdf_dataset

    lf_fill = -999.0
    filename = get_access_output_filename(year=year,month=month,day=day,satellite=satellite,dataroot=dataroot)

    root_grp = netcdf_dataset(filename, 'a', format='NETCDF4')
    lf = root_grp.createVariable('land_fraction', 'f4', ('latitude', 'longitude',),zlib=True)

    lf.standard_name = 'land_fraction'
    lf.long_name = 'land_fraction'
    lf.FillValue = -999.0
    lf.renameAttribute('FillValue', '_FillValue')
    lf.missing = -999.0
    lf.valid_min = 0.0
    lf.valid_max = 1.0
    lf.coordinates = 'latitude longitude'

    lf[:,:] = land_fraction
    root_grp.close()

def write_daily_tb_netcdf(*,year,month,day,satellite,tb_array_by_hour,time_array_by_hour,dataroot='L:/access/',file_list=None):

    tb_fill = -999.0

    import numpy as np
    import os
    import datetime
    from netCDF4 import Dataset as netcdf_dataset
    

    #figure out the output file
    filename = get_access_output_filename(year=year,month=month,day=day,satellite=satellite,dataroot=dataroot)

    list_of_channels=['time','6V','6H','7V','7H','11V','11H','19V','19H',
                     '24V','24H','37V','37H','89V','89H']
    
    lats = np.arange(0,721)*0.25 - 90.0
    lons = np.arange(0,1440)*0.25 

    root_grp = netcdf_dataset(filename, 'w', format='NETCDF4')

    root_grp.Conventions = "CF-1.8"
    root_grp.standard_name_vocabulary = "CF Standard Name Table (v78, 21 September 2021)"
    root_grp.id = os.path.basename(filename)
    root_grp.title = f"Resampled {satellite} brightness temperatures"
    root_grp.product_version = "v00r00"
    root_grp.date_issued = "2021-10-01"
    root_grp.summary = "Remote Sensing Systems (RSS) Resampled brightness temperature; intercalibrated and homogenized brightness temperature polar-orbiting resampled to a regular Earth grid"
    root_grp.keywords = "EARTH SCIENCE > SPECTRAL/ENGINEERING > MICROWAVE > BRIGHTNESS TEMPERATURE"
    root_grp.keywords_vocabulary = "NASA Global Change Master Directory (GCMD) Earth Science Keywords, Version 6.0"
    root_grp.platform = "GCOM-W1, JAXA"
    root_grp.sensor = "AMSR2 > Advanced Microwave Scanning Radiometer 2"
    root_grp.cdm_data_type = "Grid"
    root_grp.program = "NASA ACCESS-0031 > Machine Learning Datasets for the Earth’s Natural Microwave Emission"
    if file_list is not None:
        source_string=''
        for source_file in file_list:
            source_string += f"{source_file}, "
        root_grp.source = source_string
    root_grp.date_created = datetime.datetime.now().isoformat()
    root_grp.creator_name = "Carl Mears"
    root_grp.creator_url = "http://www.remss.com/"
    root_grp.creator_email = "mears@remss.com"
    root_grp.institution = "Remote Sensing Systems"
    root_grp.processing_level = "NASA Level 4"
    root_grp.references = "None"
    root_grp.history = "Created Resampled Brightness Temperature from RSS AMSR2 L1A data"
    root_grp.geospatial_lat_min = -90.0
    root_grp.geospatial_lat_max = 90.0
    root_grp.geospatial_lon_min = 0.0
    root_grp.geospatial_lon_max = 359.9999
    root_grp.geospatial_lat_units = "degrees_north"
    root_grp.geospatial_lon_units = "degrees_east"
    root_grp.spatial_resolution = "30 km X 30 km"
    day_boundary = datetime.datetime(year,month,day,0,0,0,0)
    start_date = day_boundary - datetime.timedelta(minutes=30)
    end_date = day_boundary + datetime.timedelta(minutes=1410.0)
    root_grp.time_coverage_start = start_date.isoformat()
    root_grp.time_coverage_end = end_date.isoformat()
    root_grp.time_coverage_duration = "P24H"
    root_grp.license = "No restrictions on access or use"
    root_grp.contributor_name = "Frank Wentz, Carl Mears"
    root_grp.contributor_role = "Principal investigator and originator of input/source or antenna temperature data, Processor and author of entire driver routine which resamples RSS native brightness temperature to a fixed Earth grid";


    root_grp.createDimension('latitude', 721)
    root_grp.createDimension('longitude', 1440)
    root_grp.createDimension('hours', 24)
    root_grp.createDimension('channels', 14)
    

    latitude = root_grp.createVariable('latitude', 'f4', ('latitude',))
    longitude = root_grp.createVariable('longitude', 'f4', ('longitude',))
    hours = root_grp.createVariable('hours', 'i4', ('hours',))
    channels = root_grp.createVariable('channels', 'i4', ('channels',))

    time = root_grp.createVariable('second_since_midnight', 'i4', ('latitude', 'longitude','hours',),zlib=True)
    tbs = root_grp.createVariable('brightness_temperature', 'f4', ('latitude', 'longitude','hours','channels',),zlib=True,least_significant_digit=2)
    

    latitude.standard_name = "latitude"
    latitude.long_name = "latitude"
    latitude.units = "degrees_north"
    latitude.valid_range = (-90.0,90.0)
    latitude.FillValue = -999.0
    latitude.renameAttribute('FillValue', '_FillValue')

    longitude.standard_name = "longitude"
    longitude.long_name = "longitude"
    longitude.units = "degrees_east"
    longitude.valid_range = (0.0, 360.0)
    longitude.FillValue = -999.0
    longitude.renameAttribute('FillValue', '_FillValue')

    hours.standard_name = 'hours_since_midnight'
    hours.units = 'hours'
    hours.valid_min = 0
    hours.valid_max = 23


    time.standard_name = 'seconds_since_midnight'
    time.long_name = 'seconds_since_midnight'
    time.FillValue = -999999
    time.renameAttribute('FillValue', '_FillValue')
    time.missing = -999999
    time.valid_min = -1800
    time.valid_max = 84600  #not 86400 because last 1/2 hour is in the next day.
    time.coordinates = 'latitude longitude'


    tbs.standard_name = 'brightness_temperature'
    tbs.units = 'degrees kelvin'
    tbs.FillValue = tb_fill
    tbs.renameAttribute('FillValue', '_FillValue')
    tbs.missing = tb_fill
    tbs.valid_min = 50.0
    tbs.valid_max = 350.0
    tbs.long_name = f'resampled {satellite} brightness temperature on lat/lon grid'
    channel_names = ''
    for channel_name in list_of_channels[1:]:
        channel_names += f'{channel_name}, '
    channel_names = channel_names[0:-2]
    tbs.channel_names = channel_names
    tbs.coordinates = 'latitude longitude'

    latitude[:] = lats
    longitude[:] = lons
    hours[:] = np.arange(0,24)
    channels[:] = np.arange(1,15)

    time_to_put = np.nan_to_num(time_array_by_hour,nan=-999998.99,posinf=-999998.99,neginf=-999998.99)
    time_to_put = np.floor(time_to_put).astype(np.int32)
    time[:, :, :] =  time_to_put

    tbs_to_put = np.nan_to_num(tb_array_by_hour,nan=tb_fill,posinf=tb_fill,neginf=tb_fill).astype(np.float32)
    tbs[:,:,:,:] = tbs_to_put
    root_grp.close()