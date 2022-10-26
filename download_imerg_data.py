"""Add IMERG rain rates to an existing daily ACCESS data file."""

import datetime
from pathlib import Path
from imerg_request.imerg_requests import imerg_half_hourly_request


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Interpolate and append IMERG rainfall data to ACCESS output file. "
            "IMERG data is downloaded if required."
        )
    )
    parser.add_argument(
        "temp_root", type=Path, help="Root directory store temporary files"
    )
    parser.add_argument(
        "start_date",
        type=datetime.date.fromisoformat,
        help="First Day to process, as YYYY-MM-DD",
    )
    parser.add_argument(
        "end_date",
        type=datetime.date.fromisoformat,
        help="Last Day to process, as YYYY-MM-DD",
    )

    args = parser.parse_args()

    temp_root: Path = args.temp_root

    START_DAY = args.start_date
    END_DAY = args.end_date

    date = START_DAY
    while date <= END_DAY:
        imerg_half_hourly_request(date=date, target_path=temp_root / "imerg")
        date += datetime.timedelta(days=1)
