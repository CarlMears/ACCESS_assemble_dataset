@echo off
set access_root=L:\access\amsr2_out_GL_30
set output_root=L:\access\amsr2_out_GL_30
set temp_root=\\flux\cm\era5\1hr
set start_date=%1-01-01
set end_date=%1-12-31
set satellite=amsr2
set target_size=30
set region=global
set land_mask_source=modis
set era5_vars_to_include=-v skt tcwv tclw u10n v10n
set wind_source=era5
set version=v01r00

cd /d M:\job_access\python\dataset_assembly

python add_ERA5_2D_vars_ACCESS_output.py ^
    --access_root %output_root% ^
    --output_root %output_root% ^
    --temp_root %temp_root% ^
    --start_date %start_date% ^
    --end_date %end_date% ^
    --sensor %satellite% ^
    --target_size %target_size% ^
    --version %version% ^
    --region %region% ^
    --update ^
    -v01r00 %era5_vars_to_include%



                    
