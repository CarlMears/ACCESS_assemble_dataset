access_root=/mnt/l/access/amsr2_out_daily_test
output_root=/mnt/l/access/amsr2_out_daily_test
temp_root=/mnt/a/data/_access_temp
rtm_data_root=/mnt/a/data/_access_temp
start_date=2012-07-26
end_date=2012-07-31
satellite=amsr2
target_size=70
region=global
land_mask_source=modis
era5_vars_to_include="-v skt tclw u10n v10n tcwv"
wind_source=era5
version=test_01

# python add_land_fraction_to_ACCESS_output.py \
#                       $output_root \
#                       $temp_root \
#                       $start_date \
#                       $end_date \
#                       $satellite \
#                       $target_size \
#                       $version \
#                       $land_mask_source 

python make_daily_ACCESS_files.py \
                    "--access_root" $output_root \
                    "--temp_root" $temp_root \
                    "--start_date" $start_date \
                    "--end_date" $end_date \
                    "--sensor" $satellite \
                    "--region" $region \
                    "--target_size" $target_size \
                    "--version" $version \
                    "--update"

# python add_atmosphere_to_ACCESS_output_no_compute.py \
#                     $access_root \
#                     $output_root \
#                     $rtm_data_root \
#                     $start_date \
#                     $end_date \
#                     $satellite \
#                     $target_size \
#                     $version \
#                     "--update"
                    
# python add_ERA5_2D_vars_ACCESS_output.py \
#                     $access_root \
#                     $output_root \
#                     $temp_root \
#                     $start_date \
#                     $end_date \
#                     $satellite \
#                     $target_size \
#                     $version \
#                     $era5_vars_to_include \
#                     "--update"
                   

# python add_wind_emiss_ACCESS.py \
#                     $access_root \
#                     $output_root \
#                     $wind_source \
#                     $wind_source \
#                     $start_date \
#                     $end_date \
#                     $satellite \
#                     $target_size \
#                     $version \
#                     "--update"

# python add_imerg_rain_rate_to_ACCESS_output.py \
#                     $access_root \
#                     $output_root \
#                     $temp_root \
#                     $start_date \
#                     $end_date \
#                     $satellite \
#                     $target_size \
#                     "--update"
                    
