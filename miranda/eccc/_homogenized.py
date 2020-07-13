import calendar
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import xarray as xr

from miranda.utils import eccc_cf_ahccd_metadata


__all__ = ["convert_ahccd_fwf_files"]


def convert_ahccd_fwf_files(
    ff,
    metadata,
    variable,
    cols_specs: Optional[List[Tuple[int, int]]] = None,
    attrs: Optional[dict] = None,
) -> xr.Dataset:

    if attrs is None:
        attrs, _, _, _ = eccc_cf_ahccd_metadata(variable)
    if cols_specs is None:
        _, _, cols_specs, _ = eccc_cf_ahccd_metadata(variable)

    _, _, _, nhead = eccc_cf_ahccd_metadata(variable)

    df = pd.read_fwf(ff, header=nhead, colspecs=cols_specs)
    if "pr" in variable:
        cols = list(df.columns[0:3])
        cols = cols[0::2]
        cols.extend(list(df.columns[4::2]))
    else:
        cols = [c for c in df.columns if "Unnamed" not in c]
    df = df[cols]
    df.replace(attrs["NaN_value"], np.NaN, inplace=True)

    list(df.columns[0:2])
    for i, j in enumerate(["Year", "Month"]):
        df = df.rename(columns={df.columns[i]: j})
    start_date = f"{df['Year'][0]}-{str(df['Month'][0]).zfill(2)}-01"
    # ndays last month
    _, ndays = calendar.monthrange(df["Year"].iloc[-1], df["Month"].iloc[-1])
    end_date = f"{df['Year'].iloc[-1]}-{str(df['Month'].iloc[-1]).zfill(2)}-{str(ndays).zfill(2)}"

    time1 = pd.date_range(start=start_date, end=end_date)
    index = pd.MultiIndex.from_arrays([df["Year"], df["Month"]])
    df.index = index

    cols = [c for c in df.columns if "Year" not in c and "Month" not in c]
    df = df[cols]
    df.columns = np.arange(1, 32)
    ds = df.stack().to_frame()
    ds.index.names = ["Year", "Month", "Day"]
    # find non-valid dates
    for y in time1.year.unique():
        for m in (
            ds[ds.index.get_level_values("Year") == y]
            .index.get_level_values("Month")
            .unique()
        ):
            _, exp_ndays = calendar.monthrange(y, m)
            ndays = (
                (ds.index.get_level_values("Year") == y)
                & (ds.index.get_level_values("Month") == m)
            ).sum()
            if ndays > exp_ndays:
                raise Exception("unknown days present")
    time_ds = pd.DataFrame(
        {
            "year": ds.index.get_level_values("Year"),
            "month": ds.index.get_level_values("Month"),
            "day": ds.index.get_level_values("Day"),
        }
    )
    ds.index = pd.to_datetime(time_ds)
    ds = ds.rename(columns={0: variable})
    ds = ds.to_xarray().rename({"index": "time"})
    dsOut = xr.Dataset(coords={"time": time1})
    dsOut[variable] = ds[variable]
    dsOut[variable].attrs = attrs
    # dsOut
    metadata = metadata.to_xarray().rename({"index": "station"}).drop_vars("station")
    metadata = metadata.assign_coords(
        {"stnid": metadata["stnid"], "station_name": metadata["Station name"]}
    )
    dsOut = dsOut.assign_coords(
        {
            "lon": metadata["long (deg)"],
            "lat": metadata["lat (deg)"],
            "elevation": metadata["elev (m)"],
        }
    )
    metadata = metadata.drop_vars(["long (deg)", "lat (deg)", "elev (m)"])
    metadata = metadata.drop_vars(["stnid", "Station name"])
    if "tas" in variable:
        metadata = metadata.rename({"%Miss": "percMiss"})
    for vv in metadata.data_vars:
        dsOut[vv] = metadata[vv]
    dsOut = dsOut.assign_coords(station=dsOut.stnid)

    return dsOut
