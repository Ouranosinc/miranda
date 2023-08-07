"""Adjusted and Homogenized Canadian Clime Data module."""
from __future__ import annotations

import calendar
import logging.config
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from miranda.preprocess._data_definitions import load_json_data_mappings
from miranda.preprocess._treatments import basic_metadata_conversion
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.Logger("miranda")

__all__ = ["convert_ahccd", "convert_ahccd_fwf_file"]


def _ahccd_variable_metadata(
    variable_code: str,
    gen: int,
) -> (dict[str, int | float | str], dict, list[tuple[int, int]], int):
    """

    Parameters
    ----------
    variable_code: {"dm", "dn", "dr", "ds", "dt", "dx"}
    gen: {1, 2, 3}

    Returns
    -------
    dict[str, int or str or float], dict, list[tuple[int, int]], int
    """
    generation = {1: "First", 2: "Second", 3: "Third"}.get(gen)
    if not generation:
        raise NotImplementedError(f"Generation '{gen}' not supported")

    config = load_json_data_mappings("eccc-homogenized")
    metadata = basic_metadata_conversion("eccc-homogenized", config)

    variable_meta = metadata["variables"].get(variable_code)
    if not variable_meta:
        raise NotImplementedError(f"Variable `{variable_code}` not supported.")

    variable_name = variable_meta.get("_variable_name")
    if variable_name:
        variable_meta = {variable_name: variable_meta}
        del variable_meta[variable_name]["_variable_name"]
    else:
        variable_meta = {variable_code: variable_meta}

    header = metadata["Header"]
    # Conditional handling of global attributes based on generation
    for field in [f for f in header if f.startswith("_")]:
        if isinstance(header[field], bool):
            if header[field] and field[1:] == "variable":
                header[field[1:]] = variable_name

        elif isinstance(header[field], dict):
            attr_treatment = header[field]["generation"]
            if field in ["_citation" "_product"]:
                for attribute, value in attr_treatment.items():
                    if attribute == generation:
                        header[field[1:]] = value

        else:
            raise AttributeError(
                f"Attribute treatment configuration for field `{field}` is not properly configured. Verify JSON."
            )

        del header[field]

    return variable_meta, header


def _ahccd_station_metadata(code):
    pass


def _ahccd_column_definitions(
    variable_code: str,
) -> tuple[dict, list[tuple[int, int]], int]:
    config = load_json_data_mappings("eccc-homogenized")
    metadata = basic_metadata_conversion("eccc-homogenized", config)

    variable = metadata["variables"][variable_code]["_variable_name"]
    if variable.startswith("tas"):
        column_names = [
            "No",
            "StnId",
            "Station name",
            "Prov",
            "FromYear",
            "FromMonth",
            "ToYear",
            "ToMonth",
            "%Miss",
            "Lat(deg)",
            "Long(deg)",
            "Elev(m)",
            "Joined",
            "RCS",
        ]
        column_spaces = [(0, 5), (5, 6), (6, 8), (8, 9)]
        ii = 9
        for i in range(1, 32):
            column_spaces.append((ii, ii + 7))
            ii += 7
            column_spaces.append((ii, ii + 1))
            ii += 1
        header_row = 3

    elif variable.startswith("pr"):
        column_names = [
            "Prov",
            "Station name",
            "stnid",
            "beg yr",
            "beg mon",
            "end yr",
            "end mon",
            "lat (deg)",
            "long (deg)",
            "elev (m)",
            "stns joined",
        ]
        column_spaces = [(0, 4), (4, 5), (5, 7), (7, 8)]
        ii = 8
        for i in range(1, 32):
            column_spaces.append((ii, ii + 8))
            ii += 8
            column_spaces.append((ii, ii + 1))
            ii += 1
        header_row = 0

    else:
        raise KeyError

    column_names = {
        col.lower().split("(")[0].replace("%", "pct_").strip().replace(" ", "_"): col
        for col in list(column_names)
    }

    return column_names, column_spaces, header_row


def convert_ahccd_fwf_file(
    ff: Path | str,
    metadata: pd.DataFrame,
    variable: str,
    *,
    generation: int,
) -> xr.Dataset:
    """Convert AHCCD fixed-width files.

    Parameters
    ----------
    ff: str or Path
    metadata: pandas.DataFrame
    variable: str
    generation: int

    Returns
    -------
    xarray.Dataset
    """
    code = dict(tasmax="dx", tasmin="dn", tas="dm", pr="dt", prsn="ds", prlp="dr").get(
        variable
    )

    variable_meta, global_attrs = _ahccd_variable_metadata(code, generation)
    col_names, cols_specs, header = _ahccd_column_definitions(code)

    df = pd.read_fwf(ff, header=header, colspecs=cols_specs)
    if "pr" in variable:
        cols = list(df.columns[0:3])
        cols = cols[0::2]
        cols.extend(list(df.columns[4::2]))
        flags = list(df.columns[5::2])
        dfflags = df[flags]
    elif "tas" in variable:
        cols = [c for c in df.columns if "Unnamed" not in c]
        flags = [c for c in df.columns if "Unnamed" in c]
        dfflags = df[flags[2:]]
    else:
        raise NotImplementedError(f"Variable `{variable}` not supported.")

    df = df[cols]
    df.replace(variable_meta[variable]["NaN_value"], np.NaN, inplace=True)

    for i, j in enumerate(["Year", "Month"]):
        df = df.rename(columns={df.columns[i]: j})
    start_date = f"{df['Year'][0]}-{str(df['Month'][0]).zfill(2)}-01"

    _, ndays = calendar.monthrange(df["Year"].iloc[-1], df["Month"].iloc[-1])
    end_date = f"{df['Year'].iloc[-1]}-{str(df['Month'].iloc[-1]).zfill(2)}-{str(ndays).zfill(2)}"
    time1 = pd.date_range(start=start_date, end=end_date)

    index = pd.MultiIndex.from_arrays([df["Year"], df["Month"]])
    df.index = index
    cols = [c for c in df.columns if "Year" not in c and "Month" not in c]
    df = df[cols]
    df.columns = np.arange(1, 32)
    ds = df.stack().to_frame()
    ds = ds.rename(columns={0: variable})
    ds.index.names = ["Year", "Month", "Day"]

    dfflags.index = index
    dfflags.columns = np.arange(1, 32)
    ds_flag = dfflags.stack().to_frame()
    ds_flag = ds_flag.rename(columns={0: "flag"})
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
            if ndays > np.int(exp_ndays):
                print(f"year {y}, month {m}, ndays={ndays}, exp_ndays={exp_ndays}")
                raise RuntimeError("Unknown days present.")

    time_ds = pd.DataFrame(
        {
            "year": ds.index.get_level_values("Year"),
            "month": ds.index.get_level_values("Month"),
            "day": ds.index.get_level_values("Day"),
        }
    )

    ds.index = pd.to_datetime(time_ds)  # noqa
    ds = ds.to_xarray().rename({"index": "time"})
    ds_out = xr.Dataset(coords={"time": time1})
    for v in ds.data_vars:
        ds_out[v] = ds[v]

    ds_out[variable].attrs = variable_meta[variable]
    metadata = metadata.to_xarray().rename({"index": "station"}).drop_vars("station")
    metadata = metadata.assign_coords(
        {
            "stnid": metadata["stnid"].astype(str),
            "station_name": metadata["station_name"],
        }
    )
    ds_out = ds_out.assign_coords(station=metadata.stnid)
    metadata = metadata.drop_vars(["stnid", "station_name"])

    ds_out[f"{variable}_flag"].attrs["long_name"] = variable_meta[variable]["long_name"]
    ds_out["lon"] = metadata["long"]
    ds_out.lon.attrs["units"] = "degrees_east"
    ds_out.lon.attrs["axis"] = "X"
    ds_out["lat"] = metadata["lat"]
    ds_out.lat.attrs["units"] = "degrees_north"
    ds_out.lat.attrs["axis"] = "Y"
    ds_out["elev"] = metadata["elev"]
    ds_out.elev.attrs["units"] = "m"
    ds_out.elev.attrs["axis"] = "Z"

    metadata = metadata.drop_vars(["long", "lat", "elev"])
    for vv in metadata.data_vars:
        if metadata[vv].dtype == "O" and (variable not in vv):
            ds_out[vv] = metadata[vv].astype(str)
        else:
            ds_out[vv] = metadata[vv]
    return ds_out


def convert_ahccd(
    data_source: str | Path,
    output_dir: str | Path,
    variable: str,
    generation: int,
    merge: bool = False,
) -> None:
    """Convert Adjusted and Homogenized Canadian Climate Dataset files.

    Parameters
    ----------
    data_source: str or Path
    output_dir: str or Path
    variable: str
    generation: int
    merge: bool

    Returns
    -------
    None
    """
    output_dir = Path(output_dir).resolve().joinpath(variable)
    output_dir.mkdir(parents=True, exist_ok=True)

    code = dict(tasmax="dx", tasmin="dn", tas="dm", pr="dt", prsn="ds", prlp="dr").get(
        variable
    )

    var_meta, global_attrs = _ahccd_variable_metadata(code, generation)
    (
        col_names,
        col_spaces,
        header_row,
    ) = _ahccd_column_definitions(code)

    gen = {2: "Second", 3: "Third"}.get(generation)
    if generation == 3 and code in {"dx", "dn", "dm"}:
        station_meta = "ahccd_gen3_temperature.csv"
    elif generation == 2 and code in {"dt", "ds", "dr"}:
        station_meta = "ahccd_gen2_precipitation.csv"

    else:
        raise NotImplementedError(f"Code '{code} for generation {gen}.")
    metadata_source = (
        Path(__file__).resolve().parent.joinpath("configs").joinpath(station_meta)
    )

    if "tas" in variable:
        metadata = pd.read_csv(metadata_source, header=2)
        metadata.columns = col_names.keys()

    elif "pr" in variable:
        metadata = pd.read_csv(metadata_source, header=3)
        metadata.columns = col_names.keys()
        for index, row in metadata.iterrows():
            if isinstance(row["stnid"], str):
                metadata.loc[index, "stnid"] = metadata.loc[index, "stnid"].replace(
                    " ", ""
                )
    else:
        raise KeyError(f"{variable} does not include 'pr' or 'tas'.")

    # Convert station .txt files to netcdf
    for ff in Path(data_source).glob(f"{code}*.txt"):
        outfile = output_dir.joinpath(ff.name.replace(".txt", ".nc"))
        if not outfile.exists():
            logger.info(ff.name)

            station_id = ff.name[2:].split(".txt")[0]
            try:
                metadata_st = metadata[metadata["stnid"] == int(station_id)]
            except ValueError:
                metadata_st = metadata[metadata["stnid"] == station_id]

            if len(metadata_st) == 1:
                ds_out = convert_ahccd_fwf_file(
                    ff, metadata_st, variable, generation=generation
                )
                ds_out.attrs = global_attrs

                ds_out.to_netcdf(outfile, engine="h5netcdf")
            else:
                logger.warning(
                    f"metadata info for station {ff.name} not found : skipping"
                )
    if not merge:
        return

    # merge individual stations to single .nc file
    ncfiles = list(output_dir.glob(f"{code}*.nc"))
    outfile = output_dir.parent.joinpath(
        "merged_stations", f"ahccd_gen{generation}_{variable}.nc"
    )

    if not outfile.exists():
        logger.info("merging stations :", variable)
        ds_ahccd = xr.open_mfdataset(
            ncfiles, concat_dim="station", combine="nested"
        ).load()

        for coord in ds_ahccd.coords:
            # xarray object datatypes mix string and int (e.g. stnid) convert to string for merged nc files
            # Do not apply to datetime object
            if coord != "time" and ds_ahccd[coord].dtype == "O":
                ds_ahccd[coord] = ds_ahccd[coord].astype(str)

        for v in ds_ahccd.data_vars:
            # xarray object datatypes mix string and int (e.g. stnid) convert to string for merged nc files
            # Do not apply to flag timeseries
            if ds_ahccd[v].dtype == "O" and "flag" not in v:
                logger.info(v)
                ds_ahccd[v] = ds_ahccd[v].astype(str)

        outfile.parent.mkdir(parents=True, exist_ok=True)
        ds_ahccd.to_netcdf(outfile, engine="h5netcdf", mode="w")
        del ds_ahccd
    for nc in outfile.parent.glob("*.nc"):
        logger.info(nc)
        ds = xr.open_dataset(nc)
        logger.info(ds)
