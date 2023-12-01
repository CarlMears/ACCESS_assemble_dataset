


satellite=smap
target_size=70


temp_root=/mnt/b/data/_access_temp
rtm_data_root=/mnt/a/data/_access_temp

args=$(getopt --name "$0" --options s:e:l:v:r:h -- "$@")
eval set -- "$args"



start_date=2015-04-02
end_date=2015-04-30
look=0
version=test_01
region=global



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

land_mask_source=modis
era5_vars_to_include="skt tcwv tclw u10n v10n"

echo start_date: $start_date
echo end_date: $end_date
echo look: $look
echo version: $version
echo region: $region

region_short=GL
if [ "$region" = "north" ]; then
    region_short=NP
fi
if [ "$region" = "south" ]; then
    region_short=SP
fi


access_root="/mnt/l/access/"$satellite"_out_"$region_short"_"$target_size
output_root="/mnt/l/access/"$satellite"_out_"$region_short"_"$target_size

echo $access_root


cd /mnt/m/job_access/python/dataset_assembly
                    
python add_ERA5_2D_vars_ACCESS_output.py \
                       --access_root $output_root \
                       --output_root $output_root \
                       --temp_root $temp_root \
                       --start_date $start_date \
                       --end_date $end_date \
                       --sensor $satellite \
                       --target_size $target_size \
                       --look $look \
                       --version $version \
                       --region $region \
                       --overwrite \
                       -v $era5_vars_to_include \
