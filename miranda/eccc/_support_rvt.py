import datetime
import urllib
from pathlib import Path
from typing import Optional, Union

import geopandas as gpd

"https://api.weather.gc.ca/collections/climate-daily/items?datetime=1840-03-01%2000:00:00/2021-06-02%2000:00:00&STN_ID=10761&f=json&limit=1500000&startindex=0"

# TODO: Investigate the API definition: https://api.weather.gc.ca/collections/climate-hourly


def gather_eccc_stations(
    start_date: Optional[datetime.datetime] = None,
    end_date: Optional[datetime.datetime] = None,
    wmo_id: Optional[str] = None,
) -> gpd.GeoDataFrame:

    base_url = "https://api.weather.gc.ca/collections/climate-daily/"

    if not start_date:
        start_date = datetime.datetime(
            year=1840, month=1, day=1, hour=0, minute=0, second=0
        ).strftime("%Y-%m-%d %H:%M:%S")
    start_date_url = str(start_date).replace(" ", "%20")

    if not end_date:
        end_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_date_url = str(end_date).replace(" ", "%20")
    date_range = f"{start_date_url}/{end_date_url}"

    facets = dict(
        f="json",
        datetime=date_range,
        STN_ID=wmo_id,
        # PROVINCE_CODE=province,
        limit=1000,
        startindex=0,
    )
    # if station_id:
    #     facets["STN_ID"] = station_id
    facet_list = list()
    for k, v in facets.items():
        facet_list.append(f"{k}={v}")
    request_facets = f"items?{'&'.join(facet_list)}"
    request_url = urllib.parse.urljoin(base_url, request_facets)  # noqa

    print(request_url)
    # Use geopandas to convert the json output to a GeoDataFrame.
    gdf = gpd.read_file(request_url)

    return gdf


if __name__ == "__main__":
    target_folder = Path().cwd().joinpath("downloaded")
    target_folder.mkdir(exist_ok=True)
    data = gather_eccc_stations(
        start_date="1998-01-01 00:00:00",
        end_date="2010-12-31 00:00:00",
        wmo_id="10761",
    )
    print(data)
