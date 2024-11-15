satellite=amsr2
target_size=70
region=north
if [ $region == "global" ]; then
    region_code=GL
fi
if [ $region == "north" ]; then
    region_code=NP
fi
if [ $region == "south" ]; then
    region_code=SP
fi

access_root="/mnt/l/access/"$satellite"_out_"$region_code"_"$target_size
output_root="/mnt/l/access/"$satellite"_out_"$region_code"_"$target_size

echo $access_root

temp_root=/mnt/b/data/_access_temp
rtm_data_root=/mnt/a/data/_access_temp
start_date=2022-01-01
end_date=2024-07-31

region=north
land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=v01r00

cd /mnt/m/job_access/python/dataset_assembly

python add_land_fraction_to_ACCESS_output.py \
                       --output_root $output_root \
                       --temp_root $temp_root \
                       --start_date $start_date \
                       --end_date $end_date \
                       --sensor $satellite \
                       --target_size $target_size \
                       --version $version \
                       --region $region \
                       --lf_version $land_mask_source

cd /mnt/m/job_access/python/dataset_assembly/scripts_AMSR2
