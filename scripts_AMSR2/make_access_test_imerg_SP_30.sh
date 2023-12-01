access_root=/mnt/l/access/amsr2_out_SP_30
output_root=/mnt/l/access/amsr2_out_SP_30
temp_root=/mnt/a/data/_access_temp
rtm_data_root=/mnt/a/data/_access_temp
start_date=2012-08-01
end_date=2021-12-31
satellite=amsr2
target_size=30
region=south
land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=v01_r00


python add_imerg_rain_rate_to_ACCESS_output.py \
                    --access_root $access_root \
                    --output_root $output_root \
                    --temp_root $temp_root \
                    --start_date $start_date \
                    --end_date $end_date \
                    --sensor $satellite \
                    --footprint_diameter $target_size \
                    --region $region \
                    --overwrite
                    
