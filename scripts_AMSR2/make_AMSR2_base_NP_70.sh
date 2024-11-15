access_root=/mnt/l/access/AMSR2_out_NP_70
output_root=/mnt/l/access/AMSR2_out_NP_70
temp_root=/mnt/b/data/_access_temp
rmt_data_root=/mnt/a/data/_access_temp
tb_orbit_root=/mnt/l/access
start_date=$1-01-01
end_date=$1-12-31
satellite=amsr2
ksat=$2
target_size=70
region=north

land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=v01r00
echo $ksat
echo $start_date
echo $end_date
cd /mnt/m/job_access/python/dataset_assembly
python make_daily_ACCESS_files.py \
                        --access_root $output_root \
                        --temp_root $temp_root \
                        --tb_orbit_root $tb_orbit_root \
                        --start_date $start_date \
                        --end_date $end_date \
                        --sensor $satellite \
                        --ksat $ksat \
                        --target_size $target_size \
                        --version $version \
                        --region $region \
                        --update 
