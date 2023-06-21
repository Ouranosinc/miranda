from __future__ import annotations

import os
from urllib.error import HTTPError

import pandas as pd
import xarray as xr


def load_station_metadata(meta: str | os.PathLike) -> xr.Dataset:
    if meta:
        df_inv = pd.read_csv(meta, header=0)
    else:
        try:
            import geopandas as gpd

            station_metadata_url = "https://api.weather.gc.ca/collections/climate-stations/items?f=json&limit=15000000"
            df_inv = gpd.read_file(station_metadata_url)
        except HTTPError as err:
            raise RuntimeError(
                f"Station metadata table unable to be fetched. Considering downloading directly: {err}"
            )
    df_inv["LONGITUDE"] = df_inv.geometry.x
    df_inv["LATITUDE"] = df_inv.geometry.y
    df_inv["ELEVATION"] = df_inv.ELEVATION.astype(float)
    df_inv["CLIMATE_IDENTIFIER"] = df_inv["CLIMATE_IDENTIFIER"].astype(str)

    df_inv = df_inv.drop(["geometry"], axis=1)
    return df_inv.to_xarray()
