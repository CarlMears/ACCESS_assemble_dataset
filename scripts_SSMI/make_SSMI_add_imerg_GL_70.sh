access_root=/mnt/l/access/ssmi_out_GL_70
output_root=/mnt/l/access/ssmi_out_GL_70
temp_root=/mnt/flux-write/imerg
start_date=$1-01-01
end_date=$1-12-31
satellite=ssmi
ksat=$2
target_size=70
region=global
version=v01r00

cd /mnt/m/job_access/python/dataset_assembly

python add_imerg_rain_rate_to_ACCESS_output.py \
                    --access_root $access_root \
                    --output_root $output_root \
                    --temp_root $temp_root \
                    --start_date $start_date \
                    --end_date $end_date \
                    --sensor $satellite \
                    --ksat $ksat \
                    --footprint_diameter $target_size \
                    --region $region \
                    --update
                    
