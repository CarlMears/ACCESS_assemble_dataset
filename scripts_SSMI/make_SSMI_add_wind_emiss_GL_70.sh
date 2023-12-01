access_root=/mnt/l/access/ssmi_out_GL_70
output_root=/mnt/l/access/ssmi_out_GL_70
temp_root=/mnt/b/data/_access_temp
rtm_data_root=/mnt/a/data/_access_temp
start_date=$1-01-01
end_date=$1-12-31
satellite=ssmi
ksat=15
target_size=70
region=global
land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=test_01

cd /mnt/m/job_access/python/dataset_assembly

python add_wind_emiss_ACCESS.py \
                    --access_root $access_root \
                    --output_root $output_root \
                    --wind $wind_source \
                    --sst $wind_source \
                    --start_date $start_date \
                    --end_date $end_date \
                    --sensor $satellite \
                    --ksat $ksat \
                    --target_size $target_size \
                    --region $region \
                    --version $version \
                    "--overwrite"
                    
