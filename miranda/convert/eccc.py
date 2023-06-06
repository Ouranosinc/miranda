"""Environment and Climate Change Canada Data Conversion module."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from ._data_corrections import dataset_corrections

__all__ = ["convert_canswe"]


def convert_canswe(file: str | Path, output: str | Path):
    """Convert the CanSWE netCDF files to production-ready CF-compliant netCDFs."""
    ds = xr.open_dataset(file)
    ds = ds.set_coords(
        [
            "lon",
            "lat",
            "elevation",
            "source",
            "station_name",
            "station_name_sec",
            "station_name_ter",
            "station_id_sec",
            "station_id_ter",
            "type_mes",
        ]
    )

    def clean_flags(variable):
        values = list(map(bytes.decode, np.sort(pd.unique(variable.values.flatten()))))
        values[0] = "n"
        mean_dict = parse_desc(variable.description)
        meanings = " ".join(np.array([mean_dict[v] for v in values]))
        return dict(flag_values=values, flag_meanings=meanings)

    def parse_desc(desc):
        d = dict(
            map(
                lambda kv: (kv[0].strip(), "_".join(kv[1].replace(">", "").split())),
                map(lambda s: s.split(": "), desc.replace(",", ";").split("; ")),
            )
        )
        d.update(n="nodata")
        return d

    for flag_var in ["data_flag_snd", "data_flag_snw", "qc_flag_snw", "qc_flag_snd"]:
        std_name = "status_flag" if "data" in flag_var else "quality_flag"
        ds[flag_var].attrs.update(standard_name=std_name, **clean_flags(ds[flag_var]))
        ds[flag_var] = ds[flag_var].where(ds[flag_var] != b"", b"n").astype(str)

    ds.snd.attrs["ancillary_variables"] = "data_flag_snd qc_flag_snd"
    ds.snw.attrs["ancillary_variables"] = "data_flag_snw qc_flag_snw"

    ds = dataset_corrections(ds, "ec-canswe")
    ds.attrs["frequency"] = "day"
    date = "-".join(ds.indexes["time"][[0, -1]].strftime("%Y%m"))
    for var in ["snd", "snw"]:
        ds[[var, f"data_flag_{var}", f"qc_flag_{var}"]].to_netcdf(
            Path(output) / f"{var}_CanSWE_day_{date}.nc"
        )
    ds[["sd"]].to_netcdf(Path(output) / f"sd_CanSWE_day_{date}.nc")
