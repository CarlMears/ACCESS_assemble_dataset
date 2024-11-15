
access_root=/mnt/l/access/amsr2_out_NP_30
output_root=/mnt/l/access/amsr2_out_NP_30
temp_root=/mnt/b/data/_access_temp
rmt_data_root=/mnt/a/data/_access_temp
tb_orbit_root=/mnt/l/access
start_date=2022-01-01
end_date=2024-07-31
satellite=amsr2
ksat=13
target_size=30
region=north

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
