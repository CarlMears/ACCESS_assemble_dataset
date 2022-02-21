"""
Code designed to download 30 minute IMERG data for a given day along with the
IMERG data from 23:30 UTC on the previous day and 00:30 UTC on the following
day.

The IMERG data are queried using the CMR EarthData API and subsequently
downloaded. The CMR API is documented at the following URL:
https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html. No login
information is required to query the data however, an EarthData login is
required to download the data. Details on how make a login as well as how to set
up a required .netrc file which works with Python requests is found at the
following URL: https://disc.gsfc.nasa.gov/data-access#python-requests.

Written by AManaster
"""

import datetime
import time
from functools import lru_cache
from pathlib import Path
from typing import Any

import requests


@lru_cache
def get_ids() -> list[str]:
    """Obtain 'concept-ids' for the three half-hourly IMERG products.

    The three products are: final, late, and early. The function results are
    cached since the results are assumed to be the same over an entire run.

    This prevents us from having to hard-code 'concept-ids' for the different
    products since these IDs can potentially change as new IMERG version are
    released.
    """
    COLLECTION_URL = "https://cmr.earthdata.nasa.gov/search/collections"

    params = {"keyword": "imerg"}
    # this search does not like defining the json version as 1.6.4
    headers = {"Accept": "application/vnd.nasa.cmr.umm_results+json"}
    first_response = requests.get(COLLECTION_URL, headers=headers, params=params)
    response_list = first_response.json()

    ids = ["" for _ in range(3)]
    for item in response_list["items"]:
        granule_meta = item["meta"]
        granule_umm = item["umm"]

        # We can expect that, even as versions of half-hourly IMERG change,
        # their native IDs will all still contain 'GPM_3IMERGHH'. Similarly, the
        # conditional statements below operate under the assumption that the
        # 'Final', 'Late', and 'Early' keywords will be present in the different
        # half-hourly IMERG products and will not change from version to
        # version. While it's possible that these things could change, it
        # remains safer than assuming that the 'concept-ids' of these products
        # will remain consistent from version to version.
        if "GPM_3IMERGHH" in granule_meta["native-id"]:
            if "Final" in granule_umm["EntryTitle"]:
                ids[0] = granule_meta["concept-id"]
            elif "Late" in granule_umm["EntryTitle"]:
                ids[1] = granule_meta["concept-id"]
            elif "Early" in granule_umm["EntryTitle"]:
                ids[2] = granule_meta["concept-id"]
            else:
                raise Exception("GPM IMERG Half Hourly Product not recognized")

    return ids


def _parse_umm(granule_umm: dict[str, Any]) -> str:
    # Simple function to parse the data download URL
    # from the .json.umm file

    for url in granule_umm["RelatedUrls"]:
        if url["Type"] == "GET DATA":
            hdf5_url: str = url["URL"]
            return hdf5_url
    else:
        raise Exception("Didn't find a granule download URL in the metadata")


def query_one_day_imerg(date: datetime.date) -> list[str]:
    """Query CMR for one day of IMERG data.

    Return a list of URLs for the daily data we want to download.
    """
    # Base URL for CMR API query
    CMR_URL = "https://cmr.earthdata.nasa.gov/search/granules"

    # Since we want data from the last half hourly file of the day prior and the
    # first half hourly file of the day after the current day
    day_before = date - datetime.timedelta(days=1)
    day_after = date + datetime.timedelta(days=1)

    # Collection IDs for the three IMERG datasets of interest
    id_list = get_ids()

    # check availability of all IMERG half-hourly products.  Should only be
    # relevant for NRT ACCESS applications.
    for id in id_list:
        params = {
            "collection_concept_id": id,
            "temporal": f"{day_before}T23:30:00Z,{day_after}T00:30:00Z",
            "sort_key": "start_date",
            "page_size": "50",
        }
        headers = {"Accept": "application/vnd.nasa.cmr.umm_results+json; version=1.6.4"}
        response = requests.get(CMR_URL, headers=headers, params=params)

        # Sometimes this query will return an empty response body which leads to
        # an error. I believe this is an issue on CMR's end since waiting
        # 30s-60s and trying again will often rectify the issue. If this
        # happens, the following block of code waits 30 seconds and tries again
        # until we receive a successful API response w/ body (i.e., status code
        # 200). I imagine there is a more elegant way to do this since this has
        # the potential to get caught if the response.status_code remains !=
        # 200.
        while response.status_code != 200:
            try:
                response = requests.get(CMR_URL, headers=headers, params=params)
            except requests.HTTPError:
                print(f"API Query error {response.status_code}. Retrying in 30s")
                time.sleep(30)

        response_list = response.json()

        # We want to look for the 'final' versions of IMERG data first
        # The following checks to see if any 'final' data exists for a day.
        # If not, it queries the next type of data (late) to see if a full day has been
        # processed.  If not, the code looks for 'early' data since it will
        # almost always have more data than the 'late' product.
        # This will likely only be relevant for NRT applications.
        if response_list["hits"] == 0:
            print(f"No data for ID {id}; Checking next")
            continue
        elif response_list["hits"] > 0 and response_list["hits"] < 51:
            # Should have maximum 51 'hits' in a day
            print("Some data, but not full day")
            if id != id_list[2]:  # if id does not equal the early ID
                continue
            break
        elif response_list["hits"] > 51:
            # If number of 'hits' is greater than 51 for some reason
            print(
                f"Number of hits exceeds maximum: {response_list['hits']}. "
                "Check downloaded files."
            )
            break
        else:
            break

    files = []
    for item in response_list["items"]:
        if item["meta"]["concept-type"] != "granule":
            raise Exception("Unexpected concept type for granule")
        granule_umm = item["umm"]
        hdf5_file = _parse_umm(granule_umm)
        files.append(hdf5_file)

    return files


def try_download(file_url: str, target_path: Path) -> Path:
    """Download IMERG file to a target directory.

    Return the path to the downloaded file.
    """

    file = file_url.split("/")[-1]
    target = target_path / file

    if target.exists():
        print(f"File: {target} already exists, skipping")
        return target
    else:
        print(f"Getting: {target}")
        print(file_url)

        result = requests.get(file_url)

        try:
            result.raise_for_status()
            target.write_bytes(result.content)
        except requests.HTTPError as e:
            print(f"requests.get() returned an error code {result.status_code}")
            raise e
        else:
            print(f"contents of URL written to {target}")
            return target


def imerg_half_hourly_request(date: datetime.date, target_path: Path) -> list[Path]:
    files = query_one_day_imerg(date=date)

    files_in_day = []
    # loop through all files and download
    for file in files:
        try:
            result = try_download(
                file_url=file,
                target_path=target_path,
            )
        except requests.HTTPError:
            continue
        else:
            files_in_day.append(result)

    return files_in_day


if __name__ == "__main__":
    date = datetime.date(2021, 12, 15)

    target_path = Path.cwd()

    files = imerg_half_hourly_request(date=date, target_path=target_path)
