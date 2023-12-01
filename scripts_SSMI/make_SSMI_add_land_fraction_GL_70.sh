access_root=/mnt/l/access/ssmi_out_GL_70
output_root=/mnt/l/access/ssmi_out_GL_70
temp_root=/mnt/b/data/_access_temp
rtm_data_root=/mnt/a/data/_access_temp
start_date=$1-01-01
end_date=$1-12-31
satellite=ssmi
ksat=$2
target_size=70
region=global
land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=test_01

cd /mnt/m/job_access/python/dataset_assembly

python add_land_fraction_to_ACCESS_output.py \
                       --output_root $output_root \
                       --temp_root $temp_root \
                       --start_date $start_date \
                       --end_date $end_date \
                       --sensor $satellite \
                       --ksat $ksat \
                       --target_size $target_size \
                       --version $version \
                       --region $region \
                       --lf_version $land_mask_source

# python make_daily_ACCESS_files.py \
#                     $output_root \
#                     $temp_root \
#                     $start_date \
#                     $end_date \
#                     $satellite \
#                     $target_size \
#                     $version \
#                     "--update"

# python add_atmosphere_to_ACCESS_output_no_compute.py \
#                     $access_root \
#                     $output_root \
#                     $rmt_data_root \
#                     $start_date \
#                     $end_date \
#                     $satellite \
#                     $target_size \
#                     $version \
#                     "--update" 

# python add_atmosphere_to_ACCESS_output_no_compute.py \
#                     --access_root $access_root \
#                     --output_root $output_root \
#                     --temp_root $rtm_data_root \
#                     --start_date $start_date \
#                     --end_date $end_date \
#                     --sensor $satellite \
#                     --target_size $target_size \
#                     --version $version \
#                     --region $region \
#                     --overwrite
                    
# python add_ERA5_2D_vars_ACCESS_output.py \
#                        --access_root $output_root \
#                        --output_root $output_root \
#                        --temp_root $temp_root \
#                        --start_date $start_date \
#                        --end_date $end_date \
#                        --sensor $satellite \
#                        --ksat $ksat \
#                        --target_size $target_size \
#                        --version $version \
#                        --region $region \
#                        --update \
#                        -v $era5_vars_to_include \

             

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
                    
