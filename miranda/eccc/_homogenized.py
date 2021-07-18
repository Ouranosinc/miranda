import calendar
import logging.config
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar

from miranda.utils import scripting

from ._utils import ahccd_metadata

logging.config.dictConfig(scripting.LOGGING_CONFIG)

logger = logging.Logger("miranda")

__all__ = ["convert_ahccd", "convert_ahccd_fwf_files"]


def convert_ahccd(
    data_source: Union[str, Path],
    output_dir: Union[str, Path],
    variable: str,
    generation: Optional[int] = None,
):
    code = dict(tasmax="dx", tasmin="dn", tas="dm", pr="dt", prsn="ds", prlp="dr").get(
        variable
    )
    var, col_names, col_spaces, header_row, global_attrs = ahccd_metadata(
        code, generation
    )
    gen = {2: "Second", 3: "Third"}.get(generation)
    if gen == 3 and code in {"dx", "dn", "dm"}:
        meta = "ahccd_gen3_temperature.csv"
        long_name_dict = dict(
            elev="elevation",
            frommonth="from month",
            fromyear="from year",
            prov="province",
            station="station identification number",
            station_name="station name",
            stnid="station identification number",
            joined="joined station (y/n)",
            tomonth="to month",
            toyear="to year",
            pct_miss="%Miss",
        )

    elif gen == 2 and code in {"dt", "ds", "dr"}:
        meta = "ahccd_gen2_precipitation.csv"
        long_name_dict = dict(
            elev="elevation",
            frommonth="from month",
            fromyear="from year",
            prov="province",
            station="station identification number",
            station_name="station name",
            stnid="station identification number",
            stns_joined="joined station (y/n)",
            tomonth="to month",
            toyear="to year",
        )

    else:
        raise NotImplementedError()
    metadata_source = Path().cwd().joinpath("data").joinpath(meta)

    if "tas" in variable and not generation:
        outvar = "temperature"
        metadata = pd.read_excel(metadata_source, header=2)
        metadata.columns = col_names
        cols_specs = col_spaces

    elif "pr" in variable and not generation:
        outvar = "precipitation"
        metadata = pd.read_excel(metadata_source, header=3)
        metadata.columns = col_names
        cols_specs = col_spaces
        for index, row in metadata.iterrows():
            if type(row["stnid"]) == str:
                metadata.loc[index, "stnid"] = metadata.loc[index, "stnid"].replace(
                    " ", ""
                )
    else:
        raise KeyError()

    output_dir.joinpath(variable).mkdir(parents=True, exist_ok=True)

    # Convert station .txt files to netcdf
    for ff in Path(data_source).joinpath(variable).glob("*d*.txt"):
        outfile = output_dir.joinpath(variable, ff.name.replace(".txt", ".nc"))
        if not outfile.exists():
            logger.info(ff.name)

            stid = ff.name.replace(var, "").split(".txt")[0]
            try:
                metadata_st = metadata[metadata["stnid"] == int(stid)]
            except ValueError:
                metadata_st = metadata[metadata["stnid"] == stid]

            if len(metadata_st) == 1:
                ds_out = convert_ahccd_fwf_files(
                    ff, metadata_st, variable, generation, cols_specs, code
                )
                ds_out.attrs = global_attrs[outvar]

                ds_out.to_netcdf(outfile)
            else:
                logger.warning(
                    f"metadata info for station {ff.name} not found : skipping"
                )

    # merge individual stations to single .nc file
    # variable
    ncfiles = list(output_dir.joinpath(variable).glob("*.nc"))
    outfile = output_dir.parent.joinpath(
        "merged_stations", f"ahccd_{gen}_{variable}.nc"
    )

    if not outfile.exists():
        logger.info("merging stations :", variable)
        with ProgressBar():
            ds_ahccd = xr.open_mfdataset(
                ncfiles, concat_dim="station", combine="nested"
            ).load()

            for coord in ds_ahccd.coords:
                # xarray object datatypes mix string and int (e.g. stnid) convert to string for merged nc files
                # Do not apply to datetime object
                if coord != "time" and ds_ahccd[coord].dtype == "O":
                    print(coord)
                    ds_ahccd[coord] = ds_ahccd[coord].astype(str)

            for v in ds_ahccd.data_vars:
                # xarray object datatypes mix string and int (e.g. stnid) convert to string for merged nc files
                # Do not apply to flag timeseries
                if ds_ahccd[v].dtype == "O" and "flag" not in v:
                    logger.info(v)
                    ds_ahccd[v] = ds_ahccd[v].astype(str)

            ds_ahccd[f"{variable}_flag"].attrs[
                "long_name"
            ] = f"{ds_ahccd[f'{variable}'].attrs['long_name']} flag"
            ds_ahccd.lon.attrs["units"] = "degrees_east"
            ds_ahccd.lon.attrs["long_name"] = "longitude"
            ds_ahccd.lat.attrs["units"] = "degrees_north"
            ds_ahccd.lat.attrs["long_name"] = "latitude"

            for ll in long_name_dict:
                ds_ahccd[ll].attrs["long_name"] = long_name_dict[ll]

            outfile.parent.mkdir(parents=True, exist_ok=True)
            ds_ahccd.to_netcdf(outfile, format="NETCDF4_CLASSIC", mode="w")

            del ds_ahccd
    for nc in outfile.parent.glob("*.nc"):
        logger.info(nc)
        ds = xr.open_dataset(nc)
        logger.info(ds)


def convert_ahccd_fwf_files(
    ff,
    metadata,
    variable,
    generation: int = None,
    cols_specs: Optional[List[Tuple[int, int]]] = None,
    attrs: Optional[dict] = None,
) -> xr.Dataset:

    code = dict(tasmax="dx", tasmin="dn", tas="dm", pr="dt", prsn="ds", prlp="dr").get(
        variable
    )

    if attrs is None:
        attrs, _, _, _, _ = ahccd_metadata(code, generation)
    if cols_specs is None:
        _, _, cols_specs, _, _ = ahccd_metadata(code, generation)
    _, _, _, nhead, _ = ahccd_metadata(code, generation)

    df = pd.read_fwf(ff, header=nhead, colspecs=cols_specs)
    if "pr" in variable:
        cols = list(df.columns[0:3])
        cols = cols[0::2]
        cols.extend(list(df.columns[4::2]))
        flags = list(df.columns[5::2])
        dfflags = df[flags]
    else:
        cols = [c for c in df.columns if "Unnamed" not in c]
        flags = [c for c in df.columns if "Unnamed" in c]
        dfflags = df[flags[2:]]

    df = df[cols]
    df.replace(attrs["NaN_value"], np.NaN, inplace=True)

    for i, j in enumerate(["Year", "Month"]):
        df = df.rename(columns={df.columns[i]: j})
    start_date = f"{df['Year'][0]}-{str(df['Month'][0]).zfill(2)}-01"

    _, ndays = calendar.monthrange(df["Year"].iloc[-1], df["Month"].iloc[-1])
    end_date = f"{df['Year'].iloc[-1]}-{str(df['Month'].iloc[-1]).zfill(2)}-{str(ndays).zfill(2)}"
    time1 = pd.date_range(start=start_date, end=end_date)

    index = pd.MultiIndex.from_arrays([df["Year"], df["Month"]])
    df.index = index
    dfflags.index = index
    cols = [c for c in df.columns if "Year" not in c and "Month" not in c]
    df = df[cols]
    df.columns = np.arange(1, 32)
    dfflags.columns = np.arange(1, 32)
    ds = df.stack().to_frame()
    ds = ds.rename(columns={0: variable})
    ds_flag = dfflags.stack().to_frame()
    ds_flag = ds_flag.rename(columns={0: "flag"})
    ds.index.names = ["Year", "Month", "Day"]
    ds_flag.index.names = ["Year", "Month", "Day"]
    ds[f"{variable}_flag"] = ds_flag["flag"]
    del ds_flag

    # find invalid dates
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

    ds = ds.to_xarray().rename({"index": "time"})

    ds_out = xr.Dataset(coords={"time": time1})
    for v in ds.data_vars:
        ds_out[v] = ds[v]

    ds_out[variable].attrs = attrs
    # ds_out
    metadata = metadata.to_xarray().rename({"index": "station"}).drop_vars("station")
    metadata = metadata.assign_coords(
        {
            "stnid": metadata["stnid"].astype(str),
            "station_name": metadata["station_name"],
        }
    )
    # ds_out = ds_out.assign_coords({'lon': metadata['long'], 'lat': metadata['lat'], 'elevation': metadata['elev']})
    #
    ds_out = ds_out.assign_coords(station=metadata.stnid)
    metadata = metadata.drop_vars(["stnid", "station_name"])

    ds_out["lon"] = metadata["long"]
    ds_out["lon"].attrs["units"] = "degrees_east"
    ds_out["lat"] = metadata["lat"]
    ds_out["lat"].attrs["units"] = "degrees_north"
    ds_out["elev"] = metadata["elev"]
    ds_out["elev"].attrs["units"] = "m"

    metadata = metadata.drop_vars(["long", "lat", "elev"])
    for vv in metadata.data_vars:
        if metadata[vv].dtype == "O" and (variable not in vv):
            ds_out[vv] = metadata[vv].astype(str)
        else:
            ds_out[vv] = metadata[vv]
    return ds_out
