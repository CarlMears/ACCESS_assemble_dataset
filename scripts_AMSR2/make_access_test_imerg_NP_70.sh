access_root=/mnt/l/access/amsr2_out_NP_70
output_root=/mnt/l/access/amsr2_out_NP_70
temp_root=/mnt/a/data/_access_temp
rtm_data_root=/mnt/a/data/_access_temp
start_date=2012-08-01
end_date=2021-12-31
satellite=amsr2
target_size=70
region=north
land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=test_01


python add_imerg_rain_rate_to_ACCESS_output.py \
                    $access_root \
                    $output_root \
                    $temp_root \
                    $start_date \
                    $end_date \
                    $satellite \
                    $target_size \
                    $region \
                    "--update"
                    
