access_root=/mnt/l/access/smap_out_GL_70
output_root=/mnt/l/access/smap_out_GL_70
temp_root=/mnt/b/data/_access_temp/smap_tb_orbits
tb_orbit_root=/mnt/b/data/_access_temp/smap_tb_orbits
rmt_data_root=/mnt/a/data/_access_temp
start_date=$1-01-01
end_date=$1-12-31
satellite=smap
target_size=70
region=global
look=1
land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=test_01

cd /mnt/m/job_access/python/dataset_assembly

python make_daily_ACCESS_files.py \
                        --access_root $output_root \
                        --temp_root $temp_root \
                        --tb_orbit_root $tb_orbit_root \
                        --start_date $start_date \
                        --end_date $end_date \
                        --sensor $satellite \
                        --target_size $target_size \
                        --look $look \
                        --version $version \
                        --region $region \
                        --update 
