access_root=/mnt/ops1p-ren/l/access/ssmi_out_GL_70
output_root=/mnt/ops1p-ren/l/access/ssmi_out_GL_70
temp_root=/mnt/b/data/_access_temp
rmt_data_root=/mnt/a/data/_access_temp
start_date=2012-01-01
end_date=2012-12_31
satellite=ssmi
ksat=15
target_size=70
region=global
land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=test_01

python make_daily_ACCESS_files.py \
                        --access_root $output_root \
                        --temp_root $temp_root \
                        --start_date $start_date \
                        --end_date $end_date \
                        --sensor $satellite \
                        --ksat $ksat \
                        --target_size $target_size \
                        --version $version \
                        --region $region \
                        --update 
