access_root=/mnt/ops1p-ren/l/access/amsr2_out_test3
output_root=/mnt/ops1p-ren/l/access/amsr2_out_test3
temp_root=/mnt/ops1p-ren/l/access/_temp
rmt_data_root=/mnt/a/data/_access_temp
start_date=2018-01-02
end_date=2018-01-02
satellite=amsr2
target_size=30
land_mask_source=modis
era5_vars_to_include="-v skt tclw u10n v10n tcwv"
wind_source=era5
version=test_01

python make_daily_ACCESS_files.py \
                    $output_root \
                    $temp_root \
                    $start_date \
                    $end_date \
                    $satellite \
                    $target_size \
                    $version \
                    --overwrite

python add_atmosphere_to_ACCESS_output_no_compute.py \
                    $access_root \
                    $output_root \
                    $rmt_data_root \
                    $start_date \
                    $end_date \
                    $satellite \
                    $target_size \
                    $version \
                    --overwrite

python add_land_fraction_to_ACCESS_output.py \
                    $output_root \
                    $temp_root \
                    $start_date \
                    $end_date \
                    $satellite \
                    $target_size \
                    $version \
                    $land_mask_source \
                    --overwrite

python add_ERA5_2D_vars_ACCESS_output.py \
                    $access_root \
                    $output_root \
                    $temp_root \
                    $start_date \
                    $end_date \
                    $satellite \
                    $target_size \
                    $version \
                    $era5_vars_to_include \
                    --overwrite


python add_wind_emiss_ACCESS.py \
                    $access_root \
                    $output_root \
                    $wind_source \
                    $wind_source \
                    $start_date \
                    $end_date \
                    $satellite \
                    $target_size \
                    $version \
                    --overwrite

python add_imerg_rain_rate_to_ACCESS_output.py \
                    $access_root \
                    $output_root \
                    $temp_root \
                    $start_date \
                    $end_date \
                    $satellite \
                    $target_size \
                    --overwrite