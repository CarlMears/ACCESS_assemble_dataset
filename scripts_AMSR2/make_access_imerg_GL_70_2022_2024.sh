access_root=/mnt/l/access/amsr2_out_GL_70
output_root=/mnt/l/access/amsr2_out_GL_70
temp_root=/mnt/flux-write/imerg
rtm_data_root=/mnt/a/data/_access_temp
start_date=2023-01-10
end_date=2024-07-31
satellite=amsr2
target_size=70
region=global
land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=v01r00

cd /mnt/m/job_access/python/dataset_assembly


python add_imerg_rain_rate_to_ACCESS_output.py \
                    --access_root $access_root \
                    --output_root $output_root \
                    --temp_root $temp_root \
                    --start_date $start_date \
                    --end_date $end_date \
                    --sensor $satellite \
                    --footprint_diameter $target_size \
                    --region $region \
                    --update
                    
