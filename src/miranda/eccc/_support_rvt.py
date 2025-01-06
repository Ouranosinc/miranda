from __future__ import annotations

import datetime
import re
import urllib
from pathlib import Path

import pandas as pd

"https://api.weather.gc.ca/collections/climate-daily/items?datetime=1840-03-01%2000:00:00/2021-06-02%2000:00:00&STN_ID=10761&f=json&limit=1500000&startindex=0"

# TODO: Investigate the API definition: https://api.weather.gc.ca/collections/climate-hourly


# FIXME: This function is a WIP - requires work.
def gather_eccc_stations(
    timestep: str,
    start_date: datetime.datetime | str | None = None,
    end_date: datetime.datetime | str | None = None,
    climate_id: str | None = None,
) -> pd.DataFrame:
    """Collect ECCC station data from the Environment and Climate Change Canada API.

    Parameters
    ----------
    timestep : str
    start_date : datetime.datetime or str, optional
    end_date : datetime.datetime or str, optional
    climate_id : str, optional

    Returns
    -------
    pandas.DataFrame
    """
    if timestep.lower() in ["hourly", "daily"]:
        base_url = f"https://api.weather.gc.ca/collections/climate-{timestep}/"
    else:
        raise ValueError(timestep)

    dates = [start_date, end_date]
    for i, date in enumerate(dates):
        if not date:
            dates[i] = datetime.datetime(
                year=1840, month=1, day=1, hour=0, minute=0, second=0
            ).strftime("%Y-%m-%d %H:%M:%S")
        else:
            if re.match(r"^\d{4}-(0[1-9]|1[0-2])-(0[1-9]|[12]\d|3[01])$", date):
                dates[i] = datetime.datetime.fromisoformat(date).strftime(
                    "%Y-%m-%d %H:%M:%S"
                )
        dates[i] = str(dates[i]).replace(" ", "%20")
    date_range = "/".join(dates)

    facets = dict(
        f="json",
        datetime=date_range,
        CLIMATE_IDENTIFIER=climate_id,
        # PROVINCE_CODE=province,
        limit=20000,
        startindex=0,
    )
    # if station_id:
    #     facets["STN_ID"] = station_id
    facet_list = list()
    for k, v in facets.items():
        facet_list.append(f"{k}={v}")
    request_facets = f"items?{'&'.join(facet_list)}"
    request_url = urllib.parse.urljoin(base_url, request_facets)  # noqa

    # Use geopandas to convert the json output to a DataFrame.
    return pd.read_file(request_url)


if __name__ == "__main__":
    target_folder = Path().cwd().joinpath("downloaded")
    target_folder.mkdir(exist_ok=True)
    data = gather_eccc_stations(
        timestep="hourly",
        start_date="2019-01-01",
        end_date="2020-12-31",
        climate_id="7040815",
    )
