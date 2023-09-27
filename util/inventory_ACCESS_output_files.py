from pathlib import Path
import datetime


def look_for_missing_files(
    access_output_root: Path,
    start_date: datetime.date,
    end_date: datetime.date,
    list_of_files: list[str],
    footprint_size: int,
) -> dict:
    date = start_date

    while date <= end_date:
        num_good = 0
        num_missing = 0
        access_output_root_this_day = (
            access_output_root
            / f"Y{date.year:04d}"
            / f"M{date.month:02d}"
            / f"D{date.day:02d}"
        )

        for file_template in list_of_files:
            file_name = (
                f"{file_template}"
                f"{date.year:04d}_{date.month:02d}_"
                f"{date.day:02d}.{footprint_size:03d}km.nc"
            )
            path_to_file = access_output_root_this_day / file_name
            if path_to_file.is_file():
                num_good += 1
            else:
                if num_missing == 0:
                    missing_files = []

                missing_files.append(file_name)
                num_missing += 1

        print()


if __name__ == "__main__":
    start_date = datetime.date(2013, 1, 1)
    end_date = datetime.date(2020, 12, 31)
    access_output_root = Path("/mnt/l/access/amsr2_out_test3")
    list_of_files = [
        "amsr2_atm_par_era5_",
        "amsr2_land_frac_modis_",
        "amsr2_ocean_emiss_era5_",
        "amsr2_rainfall_rate_",
        "amsr2_resamp_tbs_",
        "amsr2_skt_era5_",
        "amsr2_tclw_era5_",
        "amsr2_tcwv_era5_",
        "amsr2_u10n_era5_",
        "amsr2_v10n_era5_",
    ]
    footprint_size = 30

    missing_files_dict = look_for_missing_files(
        access_output_root=access_output_root,
        start_date=start_date,
        end_date=end_date,
        list_of_files=list_of_files,
        footprint_size=footprint_size,
    )

    print()
