"""
    Code designed to download 30 minute IMERG data for a given day
    along with the IMERG data from 23:30 UTC on the previous day
    and 00:30 UTC on the following day (for interpolation purposes).

    An EarthData login is required.  Details on how make a login as well
    as how to set up a required .netrc file which works with Python requests is
    found at the following URL: https://disc.gsfc.nasa.gov/data-access#python-requests.

"""

import requests
import datetime
from pathlib import Path


def try_download(
    *,
    date: datetime.date,
    url: str,
    target_path: Path,
    start_time: str,
    end_time: str,
    time: str,
) -> Path:
    """
        Trying multiple IMERG file types to see if any
        have data for a given time.
    """

    year = date.strftime("%Y")
    month = date.strftime("%m")
    day = date.strftime("%d")

    if "HHL" in url:
        filename = f"3B-HHR-L.MS.MRG.3IMERG.{year}{month}{day}-S{start_time}-E{end_time}.{time}.V06B.HDF5"
    elif "HHE" in url:
        filename = f"3B-HHR-E.MS.MRG.3IMERG.{year}{month}{day}-S{start_time}-E{end_time}.{time}.V06B.HDF5"
    else:
        filename = f"3B-HHR.MS.MRG.3IMERG.{year}{month}{day}-S{start_time}-E{end_time}.{time}.V06B.HDF5"

    target = target_path / filename

    if target.exists():
        print(f"File: {target} already exists, skipping")
        return target
    else:
        print(f"Getting: {target}")
        print(url + filename)

        result = requests.get(url + filename)
        try:
            result.raise_for_status()
            f = open(target, "wb")
            f.write(result.content)
            f.close()
            print(f"contents of URL written to {target}")
            return target
        except:
            error = f"{str(result.status_code)}"
            print(f"requests.get() returned an error code {error}")
            return error


def imerg_half_hourly_request(*, date: datetime.date, target_path: Path) -> Path:

    first_time = datetime.datetime.combine(
        date - datetime.timedelta(days=1), datetime.time(23, 30)
    )
    start_times = [
        first_time + datetime.timedelta(minutes=x) for x in range(0, 1500, 30)
    ]
    end_times = [
        first_time + datetime.timedelta(minutes=x) for x in range(29, 1529, 30)
    ]

    times = [f"{h:02d}".zfill(4) for h in range(0, 1440, 30)]
    times.insert(0, "1410")
    times.append("0000")
    files_in_day = []
    # times = [f"{h:02d}".zfill(4) for h in range(0, 60, 30)]

    for t in range(len(times)):

        start_time = start_times[t].strftime("%H%M00")
        end_time = end_times[t].strftime("%H%M59")
        time = times[t]

        jday = start_times[t].strftime("%j")
        yr = start_times[t].strftime("%Y")

        url_final = f"https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHH.06/{yr}/{jday}/"  # this is the URL we want to check first
        url_late = f"https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHHL.06/{yr}/{jday}/"  # this is the URL we want to check second
        url_early = f"https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHHE.06/{yr}/{jday}/"  # last URL we want to check
        urls = [url_final, url_late, url_early]
        file_exist_flag = False

        for url in urls:  # loop through all possible URLs to find IMERG file
            if file_exist_flag:
                continue

            result = try_download(
                date=start_times[t],
                url=url,
                target_path=target_path,
                start_time=start_time,
                end_time=end_time,
                time=time,
            )

            if isinstance(result, Path):
                files_in_day.append(result)
                file_exist_flag = True
            else:
                continue

    return files_in_day


if __name__ == "__main__":
    # date = datetime.date(2021, 10, 18)
    date = datetime.date(2012, 7, 11)

    target_path = Path("C:/ACCESS/output_files/_temp")

    files = imerg_half_hourly_request(date=date, target_path=target_path,)

    print(files)
    # hr1 = xr.open_dataset(files[0], group='Grid')
    # rain = hr1["precipitationCal"].values
    # lat = hr1["lat"].values
    # lon = hr1["lon"].values
    # hr1.close()
