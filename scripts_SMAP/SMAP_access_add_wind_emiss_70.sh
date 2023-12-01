satellite=smap
access_root=/mnt/l/access/smap_out_GL_70
output_root=/mnt/l/access/smap_out_GL_70
temp_root=/mnt/a/data/_access_temp

start_date=2015-04-02
end_date=2015-04-04

target_size=70
region=global
version=v01r00
wind_source=era5
ksat=13

args=$(getopt --name "$0" --options s:e:l:v:r:h -- "$@")
eval set -- "$args"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s) start_date=$2; shift 2;;
        -e) end_date=$2; shift 2;;
        -l) look=$2; shift 2;;
        -v) version=$2; shift 2;;
        -r) region=$2; shift 2;;
        -h) echo "Usage: $0 -s <start_date> -e <end_date> -l <look> -v <version> -r <region>"
            exit 0;;
        --) shift; break;;
        *)
            echo "Invalid option: $1"
            exit 1
            ;;
    esac
done

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
                    --look $look \
                    "--overwrite"

                    
