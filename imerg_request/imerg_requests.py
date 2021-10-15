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
import numpy as np
import math
import sys

def imerg_half_hourly_request(
    *, date: datetime.date, target_path: Path
    ) -> Path:

    jday = date.strftime('%j')
    year = date.strftime('%Y')
    month = date.strftime('%m')
    day = date.strftime('%d')

    url_final = f'https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHH.06/{year}/{jday}/' # this is the URL we want to check first
    url_late = f'https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHHL.06/{year}/{jday}/' # this is the URL we want to check second
    url_early = f'https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHHE.06/{year}/{jday}/'# last URL we want to check

    times = [f"{h:02d}".zfill(4) for h in range(0, 1440, 30)]
  
    for time in times:
        hour = int(time)/60

        # Might be a nicer way to do this part
        if str(hour).endswith('.5'):
            hourmin = f"{math.floor(hour)}".zfill(2)+"30"
        else:
            hourmin = f"{math.floor(hour)}".zfill(2)+"00"

        filename =f'3B-HHR.MS.MRG.3IMERG.{year}{month}{day}-S{hourmin}00-E002959.{time}.V06B.HDF5'
 
        target = target_path/filename

        if target.exists():
            print(f"File: {target} already exists, skipping")
        else:
            print(f"Getting: {target}")
            print(url_final + filename)
            result = requests.get(url_final+filename)

            try:
                result.raise_for_status()
                f = open(target, 'wb')
                f.write(result.content)
                f.close()
                print(f'contents of URL written to {target}')
            except:
                print(f'requests.get() returned an error code {str(result.status_code)}')


    return filename

if __name__ == "__main__":
    date = datetime.date(2012, 7, 11)

    target_path = Path("C:/ACCESS/output_files/_temp")

    file = imerg_half_hourly_request(
        date=date,
        target_path=target_path,
    )

    print(file)

