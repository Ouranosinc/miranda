"""Adjusted and Homogenized Canadian Clime Data module."""

from __future__ import annotations
import calendar
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from miranda.io import write_dataset
from miranda.io.utils import name_output_file
from miranda.preprocess._metadata import (
    eccc_variable_metadata,
    homogenized_column_definitions,
)
from miranda.treatments import find_project_variable_codes, load_json_data_mappings


logger = logging.Logger("miranda")

__all__ = ["convert_ahccd", "convert_ahccd_fwf_file", "merge_ahccd"]


def convert_ahccd_fwf_file(
    ff: Path | str,
    metadata: pd.DataFrame,
    variable: str,
    *,
    generation: int,
) -> xr.Dataset:
    """
    Convert AHCCD fixed-width files.

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
    configuration = load_json_data_mappings("eccc-ahccd")
    code = find_project_variable_codes(variable, configuration)

    variable_meta, global_attrs = eccc_variable_metadata(code, "eccc-ahccd", generation, configuration)
    column_names, column_spaces, column_dtypes, header = homogenized_column_definitions(code)

    df = pd.read_fwf(ff, header=header, colspecs=column_spaces, dtype=column_dtypes)

    # Handle different variable types
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

    # Extract relevant columns
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
        for m in ds[ds.index.get_level_values("Year") == y].index.get_level_values("Month").unique():
            _, exp_ndays = calendar.monthrange(y, m)
            ndays = ((ds.index.get_level_values("Year") == y) & (ds.index.get_level_values("Month") == m)).sum()
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
    metadata = metadata.assign_coords(dict(station_name=metadata["station_name"]))
    ds_out = ds_out.assign_coords(station=metadata.stnid.astype(str))
    metadata = metadata.drop_vars(["stnid", "station_name"])

    ds_out[f"{variable}_flag"].attrs["long_name"] = variable_meta[variable]["long_name"]

    ds_out["lon"] = metadata["long"]
    ds_out.lon.attrs["units"] = "degrees_east"
    ds_out.lon.attrs["axis"] = "X"
    ds_out["lat"] = metadata["lat"]
    ds_out.lat.attrs["units"] = "degrees_north"
    ds_out.lat.attrs["axis"] = "Y"
    ds_out["elev"] = metadata["elev"]
    ds_out.elev.attrs["units"] = "meters"
    ds_out.elev.attrs["positive"] = "up"
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
    *,
    generation: int,
    merge: bool = False,
    overwrite: bool = False,
) -> None:
    """
    Convert Adjusted and Homogenized Canadian Climate Dataset files.

    Parameters
    ----------
    data_source: str or Path
    output_dir: str or Path
    variable: str
    generation: int
    merge: bool
    overwrite: bool

    Returns
    -------
    None
    """
    configuration = load_json_data_mappings("eccc-ahccd")

    output_dir = Path(output_dir).resolve().joinpath(variable)
    output_dir.mkdir(parents=True, exist_ok=True)

    code = find_project_variable_codes(variable, configuration)
    variable_meta, global_attrs = eccc_variable_metadata(code, "eccc-ahccd", generation, configuration)
    (
        column_names,
        column_spaces,
        column_dtypes,
        header_row,
    ) = homogenized_column_definitions(code)

    gen = {2: "Second", 3: "Third"}.get(generation)
    if generation == 3 and code in {"dx", "dn", "dm"}:
        station_meta = "ahccd_gen3_temperature.csv"
    elif generation == 2 and code in {"dt", "ds", "dr"}:
        station_meta = "ahccd_gen2_precipitation.csv"
    else:
        raise NotImplementedError(f"Code '{code} for generation {gen}.")
    metadata_source = Path(__file__).resolve().parent.joinpath("configs").joinpath(station_meta)

    if "tas" in variable:
        metadata = pd.read_csv(metadata_source, header=2)
        metadata.columns = column_names.keys()

    elif "pr" in variable:
        metadata = pd.read_csv(metadata_source, header=3)
        metadata.columns = column_names.keys()
        for index, row in metadata.iterrows():
            if isinstance(row["stnid"], str):
                metadata.loc[index, "stnid"] = metadata.loc[index, "stnid"].replace(" ", "")
    else:
        raise KeyError(f"{variable} does not include 'pr' or 'tas'.")

    # Convert station .txt files to netcdf
    for ff in Path(data_source).glob(f"{code}*.txt"):
        output_name = ff.name.replace(".txt", ".nc")
        if not output_dir.joinpath(output_name).exists() or overwrite:
            logger.info(ff.name)

            station_id = ff.stem[2:]
            metadata_st = metadata[metadata["stnid"] == station_id]

            if len(metadata_st) == 1:
                ds_out = convert_ahccd_fwf_file(ff, metadata_st, variable, generation=generation)
                ds_out.attrs = global_attrs

                write_dataset(
                    ds_out,
                    output_dir,
                    output_format="netcdf",
                    output_name=output_name,
                    overwrite=overwrite,
                    compute=True,
                )
            else:
                msg = f"Metadata info for station {ff.name} not found: Skipping..."
                logger.warning(msg)
        else:
            msg = f"{output_name} already exists: Skipping..."
            logger.info(msg)
    if merge:
        merge_ahccd(data_source, output_dir, variable)
    return


def merge_ahccd(
    data_source: str | Path,
    output_dir: str | Path | None = None,
    variable: str | None = None,
    overwrite: bool = False,
) -> None:
    """Merge Adjusted and Homogenized Canadian Climate Dataset files."""
    configuration = load_json_data_mappings("eccc-ahccd")

    if variable:
        code = find_project_variable_codes(variable, configuration)
        glob_pattern = f"{code}*.nc"
        output_dir = Path(output_dir).resolve().joinpath(variable)
    else:
        glob_pattern = "*.nc"
        output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Merge individual stations to single .nc file
    ds_ahccd = xr.open_mfdataset(list(data_source.glob(glob_pattern)), concat_dim="station", combine="nested")

    for coord in ds_ahccd.coords:
        # xarray object datatypes mix string and int (e.g. station) convert to string for merged nc files
        # Do not apply to datetime object
        if coord != "time" and ds_ahccd[coord].dtype == "O":
            ds_ahccd[coord] = ds_ahccd[coord].astype(str)

    variables_found = set()
    for v in ds_ahccd.data_vars:
        # xarray object datatypes mix string and int (e.g. station) convert to string for merged nc files
        # Do not apply to flag timeseries
        if ds_ahccd[v].dtype == "O" and "flag" not in v:
            ds_ahccd[v] = ds_ahccd[v].astype(str)
        try:
            variables_found.add(find_project_variable_codes(str(v), configuration))
        except NotImplementedError:
            msg = f"Variable {v} not found in metadata."
            logging.info(msg)
            pass

    # Name output file
    ds_ahccd.attrs["variable"] = ", ".join(variables_found)
    if len(variables_found) > 1:
        variables = "-".join(variables_found)
        msg = f"Many variables found. Merging station and variables files in {data_source}."
        logger.info(msg)
    else:
        variables = variables_found.pop()
    output_name = name_output_file(ds_ahccd, "netcdf", variables)

    try:
        msg = f"Writing merged file to: {output_dir}."
        logger.info(msg)
        write_dataset(
            ds_ahccd,
            output_dir,
            output_format="netcdf",
            output_name=output_name,
            overwrite=overwrite,
            compute=True,
        )
        del ds_ahccd
    except FileExistsError:
        logger.info("Merged file already exists. Use overwrite=`True` to overwrite.")
