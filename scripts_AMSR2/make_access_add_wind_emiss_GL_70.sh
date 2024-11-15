satellite=amsr2
target_size=70
region=global
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
temp_root=/mnt/b/data/_access_temp
rtm_data_root=/mnt/a/data/_access_temp
start_date=2024-01-01
end_date=2024-07-31
satellite=amsr2
target_size=70
region=global
land_mask_source=modis
era5_vars_to_include="-v skt tcwv tclw u10n v10n"
wind_source=era5
version=v01r00

cd /mnt/m/job_access/python/dataset_assembly

python add_wind_emiss_ACCESS.py \
                    --access_root $access_root \
                    --output_root $output_root \
                    --wind $wind_source \
                    --sst $wind_source \
                    --start_date $start_date \
                    --end_date $end_date \
                    --sensor $satellite \
                    --target_size $target_size \
                    --region $region \
                    --version $version \
                    "--update"
                    
