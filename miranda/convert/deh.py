"""DEH Hydrograph Conversion module."""
from __future__ import annotations

import json
import logging.config
import os
import re
from pathlib import Path

import pandas as pd
import xarray as xr
from xclim.core.units import units as u

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = ["open_txt"]

# CMOR-like attributes
cmor = json.load(open(Path(__file__).parent / "data" / "deh_cf_attrs.json"))[  # noqa
    "variable_entry"
]

# TODO: Some potentially useful attributes were skipped, because they would be complicated to include in a dataset since they vary per station
meta_patterns = {
    "Station: ": "name",
    "Bassin versant: ": "bv",
    "Coordonnées: (NAD83) ": "coords",
}

data_header_pattern = "Station Date Débit (m³/s) Remarque\n"


def extract_daily(path: os.PathLike | str) -> tuple[dict, pd.DataFrame]:
    """Extract data and metadata from DEH (MELCC) stream flow file."""
    with open(path, encoding="latin1") as fh:
        txt = fh.read()
        txt = re.sub(" +", " ", txt)
        meta, data = txt.split(data_header_pattern)

    m = dict()
    for key in meta_patterns:
        # Various possible separators to take into account
        m[meta_patterns[key]] = (
            meta.split(key)[1].split(" \n")[0].split("\n")[0].split(" Régime")[0]
        )

    d = pd.read_csv(
        path,
        delimiter=r"\s+",
        skiprows=len(meta.splitlines()),
        encoding="latin1",
        converters={0: lambda x: str(x)},  # noqa
        index_col=1,
        parse_dates=True,
        infer_datetime_format=True,
    )
    if len(d["Station"].unique()) == 1:
        m["station"] = d["Station"].unique()[0]
        d = d.drop("Station", axis=1)
    else:
        raise ValueError("Multiple stations detected in the same file.")
    d = d.rename(columns={"Remarque": "Nan", "(m³/s)": "Remarque"})
    d.index.names = ["time"]
    d = d.drop("Nan", axis=1)

    return m, d


def to_cf(meta: dict, data: pd.DataFrame, cf_table: dict) -> xr.Dataset:
    """Return CF-compliant metadata."""
    ds = xr.Dataset()

    ds["q"] = xr.DataArray(data["Débit"], attrs=cf_table["q"])
    ds["flag"] = xr.DataArray(data["Remarque"], attrs=cf_table["flag"])

    ds["name"] = xr.DataArray(meta["name"])
    ds["station_id"] = xr.DataArray(meta["station"])

    ds["area"] = xr.DataArray(
        u.convert(float(meta["bv"].split(" ")[0]), meta["bv"].split(" ")[1], "km²"),
        attrs={"long_name": "drainage area", "units": "km2"},
    )

    def parse_dms(coord):
        deg, minutes, seconds, _ = re.split("[°'\"]", coord)
        if float(deg) > 0:
            return round(
                float(deg) + float(minutes) / 60 + float(seconds) / (60 * 60), 6
            )
        return round(float(deg) - (float(minutes) / 60 + float(seconds) / (60 * 60)), 6)

    coords = meta["coords"].split(" // ")
    ds["lat"] = xr.DataArray(
        parse_dms(coords[0]),
        attrs={
            "standard_name": "latitude",
            "long_name": "latitude",
            "units": "decimal_degrees",
        },
    )
    ds["lon"] = xr.DataArray(
        parse_dms(coords[1]),
        attrs={
            "standard_name": "longitude",
            "long_name": "longitude",
            "units": "decimal_degrees",
        },
    )

    ds.attrs[
        "institution"
    ] = "Ministère de l'Environnement et de la Lutte contre les changements climatiques"
    ds.attrs[
        "source"
    ] = "Hydrometric data <https://www.cehq.gouv.qc.ca/hydrometrie/historique_donnees/index.asp>"
    ds.attrs["redistribution"] = "Redistribution policy unknown. For internal use only."

    return ds


def open_txt(path: str | Path, cf_table: dict | None = cmor) -> xr.Dataset:
    """Extract daily HQ meteorological data and convert to xr.DataArray with CF-Convention attributes."""
    meta, data = extract_daily(path)
    return to_cf(meta, data, cf_table)
