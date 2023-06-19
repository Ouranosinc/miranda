"""MELCC (Québec) Weather Stations data conversion module."""
from __future__ import annotations

import datetime as dt
import logging.config
import os
import re
import warnings
from argparse import ArgumentParser
from collections.abc import Sequence
from io import StringIO
from pathlib import Path
from subprocess import run
from typing import Any

import numpy as np
import pandas as pd
import xarray
import xarray as xr
from xclim.core.formatting import update_history
from xclim.core.units import convert_units_to, pint_multiply, str2pint

from miranda import __version__
from miranda.scripting import LOGGING_CONFIG

from ._data_corrections import (
    dataset_corrections,
    load_json_data_mappings,
    metadata_conversion,
)

logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)


frequencies = {
    "a": ("min", "1 min"),  # T
    "b": ("point", ""),  # "Élémentaire"...
    "c": ("12hr", "12 h"),  # 12H
    "f": ("point", ""),  # Instantaneous??
    "h": ("hr", "1 hr"),  # H
    "q": ("day", "1 d"),  # D
    "z": ("4min", "4 min"),  # 4T
}

__all__ = [
    "parse_var_code",
    "list_tables",
    "read_table",
    "read_definitions",
    "read_stations",
    "convert_mdb",
    "convert_snow_table",
    "concat",
    "convert_melcc_obs",
]


def parse_var_code(vcode: str) -> dict[str, Any]:
    """Parse variable code to generate metadata

    Parameters
    ----------
    vcode: str

    Returns
    -------
    dict[str, Any]
    """
    match = re.match(r"(\D*)(\d*)([abcfhqz])", vcode)
    if match is None:
        raise ValueError(f"Failed to parse variable {vcode}.")
    short_code, instrument, freq = match.groups()
    melcc_code = vcode[:2].upper() + vcode[2:-1].lower() + vcode[-1].upper()
    short_code = short_code[:2].upper() + short_code[2:].lower()
    return {
        "freq": frequencies[freq][0],
        "sampling": frequencies[freq][1],
        # PI has different meanings, it's not just a different instrument.
        "var_name": short_code if short_code != "PI" else melcc_code[:-1],
        "instrument": np.int32(instrument),
        "melcc_code": melcc_code,
    }


def list_tables(db_file):
    """List the tables of an MDB file."""
    res = run(["mdb-tables", db_file], capture_output=True, encoding="utf-8")
    if res.returncode != 0:
        raise ValueError(
            f"Calling mdb-tables on {db_file} failed with code {res.returncode}: {res.stderr}"
        )
    return res.stdout.lower().strip().split()


def read_table(db_file: str | os.PathLike, tab: str | os.PathLike) -> xarray.Dataset:
    """Read a MySQL table into an xarray object.

    Parameters
    ----------
    db_file: str or os.PathLike
    tab : str or os.PathLike

    Returns
    -------
    xarray.Dataset
    """
    res = run(
        ["mdb-export", "-T", "%Y-%m-%dT%H:%M:%S", db_file, tab.upper()],
        capture_output=True,
        encoding="utf-8",
        check=True,
    )
    df = (
        pd.read_csv(
            StringIO(res.stdout),
            parse_dates=["DATE"],
            infer_datetime_format=True,
            date_parser=lambda d: pd.to_datetime(d, errors="coerce"),
        )
        .rename(
            columns={
                "NO_SEQ_STATION": "station",
                "DATE": "time",
                "CODE_STATUT_DONNEE": f"{table}_flag",
                "VALEUR_DONNEE": table,
            }
        )
        .drop(columns=["STATUT_APPROBATION"])  # I don't know what this is...
        .set_index(["station", "time"])
    )
    NaT = df.index.get_level_values("time").isnull().sum()
    if NaT > 0:
        logger.warning(
            f"{NaT} values dropped because of unparseable timestamps (probably data from before 1900)."
        )
        df = df[df.index.get_level_values("time").notnull()]
    return df[~df.index.duplicated()].to_xarray()


def read_stations(db_file: str | os.PathLike) -> pd.DataFrame:
    """Read station file using mdbtools.

    Parameters
    ----------
    db_file: str or os.PathLike

    Returns
    -------
    pandas.DataFrame
    """
    res = run(
        ["mdb-export", "-T", "%Y-%m-%dT%H:%M:%S", db_file, "STATIONs"],
        capture_output=True,
        encoding="utf-8",
        check=True,
    )
    df = pd.read_csv(
        StringIO(res.stdout),
        parse_dates=["Date_Ouverture", "Date_Fermeture"],
        infer_datetime_format=True,
        date_parser=lambda d: pd.to_datetime(d, errors="coerce"),
    ).rename(
        columns={
            "NO_STATION_CLIMATO": "station_id",
            "NOM_STATION": "station_name",
            "Latitude": "lat",
            "Longitude": "lon",
            "Altitude": "elevation",
            "Date_Ouverture": "station_opening",
            "Date_Fermeture": "station_closing",
            "Type_Poste": "station_type",
            "CODE_TYPE_POSTE": "station_type",
            "No_Seq_Station": "station",
        },
        errors="ignore",
    )
    ds = df.set_index("station").to_xarray()
    da = ds.set_coords(ds.data_vars.keys()).station
    da.lat.attrs.update(units="degree_north", standard_name="latitude")
    da["lon"] = -da.lon
    da.lon.attrs.update(units="degree_east", standard_name="longitude")
    da.elevation.attrs.update(units="m", standard_name="height")
    da["station_id"] = da["station_id"].astype(str)
    da["station_name"] = da["station_name"].astype(str)
    da["station_type"] = da["station_type"].astype(str)
    da.station_opening.attrs.update(description="Date of station creation.")
    da.station_closing.attrs.update(description="Date of station closure.")
    return da.isel(station=~da.indexes["station"].duplicated())


def read_definitions(dbfile: str):
    """Read variable definition file using mdbtools.

    Parameters
    ----------
    dbfile: str

    Returns
    -------
    pandas.DataFrame
    """
    res = run(
        ["mdb-export", dbfile, "DDs"],
        capture_output=True,
        encoding="utf-8",
        check=True,
    )
    definitions = (
        pd.read_csv(StringIO(res.stdout))
        .rename(
            columns={
                "GROUPE": "melcc_category",
                "NOM_DD": "vcode",
                "PRIORITE": "priority",
                "UNITE": "units",
                "Description": "description_fr",
            }
        )
        .drop(columns=["IND_TRANSFO"])
    )  # TODO: Was ist das?
    definitions["vcode"] = definitions.vcode.str.lower()
    definitions["units"] = definitions.units.fillna("").replace("mm+cm", "mm")
    definitions["priority"] = definitions.priority.astype(np.int32)
    return definitions.set_index("vcode")


def convert_mdb(
    database: str | Path,
    stations: xr.Dataset,
    definitions: xr.Dataset,
    output: str | Path,
    overwrite: bool = True,
) -> dict[tuple[str, str], Path]:
    """Convert microsoft databases of MELCC observation data to xarray objects.

    Parameters
    ----------
    database: str or Path
    stations
    definitions
    output
    overwrite

    Returns
    -------
    dict[tuple[str, str], Path]
    """
    outs = dict()
    tables = list_tables(database)
    for tab in tables:
        if table.startswith("gdb") or tab.startswith("~"):
            continue
        logger.info(f"Parsing {database}:{tab}.")
        meta = parse_var_code(tab)
        code = meta["melcc_code"]
        existing = list(output.glob(f"*_{code}_MELCC_{meta['freq']}_*.nc"))
        if existing and not overwrite:
            if len(existing) > 1:
                raise ValueError(f"Found more than one existing file for table {code}!")
            file = existing[0]
            vname = file.stem.split(code)[0][:-1]
            logger.info(f"File already exists {file}, skipping.")
            outs[(vname, code)] = file
            continue
        raw = read_table(database, tab)
        if np.prod(list(raw.dims.values())) == 0:
            # If any dimension is 0.
            logger.warning("The table is empty.")
            continue
        vv = meta.pop("var_name")

        raw = raw.rename({tab: vv, f"{tab}_flag": f"{vv}_flag"})
        raw.attrs["frequency"] = meta.pop("freq")
        raw[vv].attrs.update(**meta)

        if tab.lower() not in definitions.index:
            warnings.warn(
                f"The {code} variable wasn't defined in the definition table."
            )
        else:
            dd = definitions.loc[tab]
            raw[vv].attrs.update(**dd)

        raw[f"{vv}_flag"] = raw[f"{vv}_flag"].fillna(0).astype(np.int32)
        raw[f"{vv}_flag"].attrs.update(
            standard_name="status_flag",
            flag_values=np.array([0, 1, 3, 5, 7], dtype="int32"),
            flag_meanings="nodata good estimated forced trace",
            flag_meanings_fr="sansdonnée correcte estimée forcée trace",
            long_name=f"Quality flag for {code}.",
        )
        try:
            stat = stations.sel(station=raw.station)
        except KeyError:
            extra_stations = set(raw.station.values) - set(stations.station.values)
            logger.error(
                f"The following station number were not found in the station list : {sorted(list(extra_stations))}"
            )
            raw = xr.merge([raw, xr.Dataset().assign(station=stations)]).sel(
                station=raw.station
            )
        else:
            raw = raw.assign_coords(station=stat)

        ds = dataset_corrections(raw, "melcc-obs")
        new_var_name = list(
            filter(lambda k: not k.endswith("_flag"), ds.data_vars.keys())
        )[0]
        ds.attrs[
            "history"
        ] = f"[{dt.datetime.now():%Y-%m-%d %H:%M:%S}] Conversion from {database.name}:{tab} to netCDF."
        date = "-".join(ds.indexes["time"][[0, -1]].strftime("%Y%m"))
        outs[(new_var_name, code)] = (
            output / f"{new_var_name}_{code}_MELCC_{raw.attrs['frequency']}_{date}.nc"
        )
        ds.to_netcdf(outs[(new_var_name, code)])
    return outs


def convert_melcc_obs(
    metafile: str | Path,
    folder: str | Path,
    output: str | Path | None = None,
    overwrite: bool = True,
) -> dict[tuple[str, str], Path]:
    """Convert MELCC observation data to xarray data objects, returning paths.

    Parameters
    ----------
    metafile: str or Path
    folder: str or Path
    output: str or Path, optional
    overwrite: bool

    Returns
    -------
    dict[str, Path]
    """
    output = Path(output or ".")

    stations = read_stations(metafile)
    definitions = read_definitions(metafile)

    logger.info("Parsing databases.")
    folder = Path(folder)
    outs = {}
    for database in folder.glob("*.mdb"):
        outs.update(convert_mdb(database, stations, definitions, output, overwrite))
    return outs


def concat(
    files: Sequence[str | Path], output_folder: str | Path, overwrite: bool = True
) -> Path:
    """Concatenate converted weather station files.

    Parameters
    ----------
    files: sequence of str or Path
    output_folder: str or Path
    overwrite: bool

    Returns
    -------
    Path
    """
    *vv, _, melcc, freq, _ = Path(files[0]).stem.split("_")
    vv = "_".join(vv)
    logger.info(f"Concatening variables from {len(files)} files ({vv}).")
    # Magic one-liner to parse all date_start and date_end entries from the file names.
    dates_start, dates_end = list(
        zip(*[map(int, Path(file).stem.split("_")[-1].split("-")) for file in files])
    )
    outpath = (
        Path(output_folder)
        / f"{vv}_{melcc}_{freq}_{min(dates_start):06d}-{max(dates_end):06d}.nc"
    )
    if outpath.is_file() and not overwrite:
        logger.info(f"Already done in {outpath}. Skipping.")
        return outpath

    dss = dict()
    for i, file in enumerate(files):
        ds = xr.open_dataset(file, chunks={"time": 1000})
        for crd in ds.coords.values():
            crd.load()
        if list(sorted(ds.data_vars.keys())) != [vv, f"{vv}_flag"]:
            raise ValueError(
                f"Unexpected variables in {file}. Got {ds.data_vars.keys()}, expected {vv} and {vv}_flag."
            )

        priority = ds[vv].attrs["priority"]
        if priority in dss:
            raise ValueError(
                f"Variable {vv} of {file} has the same priority ({priority}) than another file of the same file list."
            )
        dss[priority] = ds

    ds_all = xr.merge(
        [ds.coords.to_dataset() for ds in dss.values()], combine_attrs="drop_conflicts"
    )
    for var in [vv, f"{vv}_flag"]:
        ds_all[var] = xr.concat(
            [ds[var] for ds in dss.values()],
            dim=xr.DataArray(list(dss.keys()), dims=("priority",), name="priority"),
            coords="all",
            combine_attrs="drop_conflicts",
        )
    ds_all = ds_all.sortby("priority")
    ds_all[f"{vv}_flag"] = ds_all[f"{vv}_flag"].where(ds_all[vv].notnull())

    ds_merged = ds_all.ffill("priority").isel(priority=-1, drop=True)
    ds_merged[f"{vv}_flag"] = ds_merged[f"{vv}_flag"].fillna(0).astype(np.int8)
    ds_merged.attrs["history"] = update_history(
        f"Different instruments merged according to their priority - miranda {__version__}",
        new_name=vv,
        **{f"priority={i}": ds for i, ds in dss.items()},
    )
    instruments = [dss[p][vv].melcc_code for p in sorted(dss)]
    ds_merged[vv].attrs[
        "melcc_code"
    ] = "Merged sources in ascending priority : " + " ,".join(map(str, instruments))

    ds_merged.attrs.update(
        source="info-climat-merged",
        title="Station observations of the MELCC - all instruments merged.",
        title_fr="Observations aux stations du MELCC - instruments fusionnés",
    )
    ds_merged.to_netcdf(outpath)
    return outpath


def convert_snow_table(file: str | Path, output: str | Path):
    """Convert snow data given through an Excel file.

    This private data is not included in the MDB files.

    Parameters
    ----------
    file : path
      The excel file with sheets:  "Stations", "Périodes standards" and "Données"
    output : path
      Folder where to put the netCDF files (one for each of snd, sd and snw).
    """
    logging.info("Parsing stations.")
    stations = pd.read_excel(file, sheet_name="Stations")
    stations = stations.rename(
        columns={
            "No": "station",
            "Nom": "station_name",
            "LAT(°)": "lat",
            "LONG(°)": "lon",
            "ALT(m)": "elevation",
            "OUVERTURE": "station_opening",
            "FERMETURE": "station_closing",
        }
    )
    statds = stations.set_index("station").to_xarray()
    stations = statds.set_coords(statds.data_vars.keys()).station
    stations["station_name"] = stations["station_name"].astype(str)
    stations.lat.attrs.update(units="degree_north", standard_name="latitude")
    stations.lon.attrs.update(units="degree_east", standard_name="longitude")
    stations.elevation.attrs.update(units="m", standard_name="height")
    stations.station_opening.attrs.update(description="Date of station creation.")
    stations.station_closing.attrs.update(description="Date of station closure.")

    # Periods
    logging.info("Parsing observation periods.")
    periods = pd.read_excel(
        file, sheet_name="Périodes standards", names=["start", "end", "middle"]
    )
    periods = periods[["start", "end"]].to_xarray()
    periods = (
        periods.to_array()
        .rename(variable="bnds", index="num_period")
        .drop_vars("bnds")
        .rename("period_bnds")
    )
    periods.attrs.update(
        description=(
            "Bounds of the sampling periods of the MELCC. "
            "Observations are taken manually once per period. "
            "The year of these bounds should be ignored."
        )
    )

    # Data
    logging.info("Parsing data.")
    data = pd.read_excel(
        file,
        sheet_name="Données",
        names=[
            "station",
            "time",
            "snd",
            "snd_flag",
            "sd",
            "sd_flag",
            "snw",
            "snw_flag",
        ],
    )
    ds = data.set_index(["station", "time"]).to_xarray()
    ds["station"] = stations.sel(station=ds.station)
    bins = periods.dt.dayofyear
    bins = np.concatenate((bins.isel(bnds=0), bins.isel(bnds=1, num_period=[-1]) + 1))
    ds["period"] = ds.time.copy(data=np.digitize(ds.time.dt.dayofyear, bins))
    ds.period.attrs.update(description="Observational period number.")

    with xr.set_options(keep_attrs=True):
        flag_attrs = dict(
            standard_name="status_flag",
            flag_values=[0, 1, 3, 5, 7],
            flag_meanings="nodata good estimated forced trace",
            flag_meanings_fr="sansdonnée correcte estimée forcée trace",
        )
        ds.snd.attrs.update(
            standard_name="surface_snow_thickness",
            units="cm",
            long_name="Snow depth",
            long_name_fr="Épaisseur de la neige au sol",
            melcc_code="NS000F",
            description_fr="Épaisseur de la neige mesurée (carottier)",
        )
        ds["snd"] = convert_units_to(ds.snd, "m")
        ds["snd_flag"] = ds.snd_flag.fillna(0).astype(int)
        ds.snd_flag.attrs.update(
            long_name="Quality of snow depth measurements.",
            long_name_fr="Qualité de la mesure d'épaisseur de la neige",
            **flag_attrs,
        )
        ds.snw.attrs.update(
            standard_name="surface_snow_amount",
            units="cm",
            long_name="Snow amount",
            long_name_fr="Quantité de neige au sol",
            melcc_code="NSQ000F",
            description_fr="Équivalent en eau de la neige mesurée (carottier)",
            description="Converted from snow water-equivalent using a water density of 1000 kg/m³",
        )
        ds["snw"] = pint_multiply(ds.snw, str2pint("1000 kg m-3"), out_units="kg m^-2")
        ds["snw_flag"] = ds.snd_flag.fillna(0).astype(int)
        ds.snw_flag.attrs.update(
            long_name="Quality of snow amount measurements.",
            long_name_fr="Qualité de la mesure de quantité de neige",
            **flag_attrs,
        )
        # Density given as a percentage of water density
        ds.sd.attrs.update(
            standard_name="surface_snow_density",
            units="%",
            long_name="Snow density",
            long_name_fr="Densité de la neige au sol",
            melcc_code="NSD000F",
            description_fr="Densité de la neige mesurée (carottier)",
        )
        ds["sd"] = pint_multiply(ds.sd, str2pint("1000 kg m-3"), out_units="kg m^-3")
        ds["sd_flag"] = ds.sd_flag.fillna(0).astype(int)
        ds.sd_flag.attrs.update(
            long_name="Quality of snow density measurements.",
            long_name_fr="Qualité de la mesure de densité de la neige",
            **flag_attrs,
        )

    ds.attrs.update(frequency="2sem")
    meta = load_json_data_mappings("melcc-snow")
    ds = metadata_conversion(ds, "melcc-snow", meta)
    date = "-".join(ds.indexes["time"][[0, -1]].strftime("%Y%m"))
    # Save
    logging.info("Saving to files.")
    for vv in ["sd", "snd", "snw"]:
        ds[[vv, f"{vv}_flag"]].to_netcdf(
            output / f"{vv}_{ds[vv].melcc_code}_MELCC_2sem_{date}.nc"
        )


if __name__ == "__main__":
    argparser = ArgumentParser(
        description="Convert MS Access data to netCDF. Requires mdbtools to be installed. "
        "By default, a single netCDF is created for each data table."
    )
    argparser.add_argument(
        "-v", "--verbose", help="Increase verbosity of the script.", action="store_true"
    )
    argparser.add_argument(
        "-c",
        "--concat",
        help="Merge all variables with the same frequency in individual datasets, according to their priority.",
        action="store_true",
    )
    argparser.add_argument(
        "-o", "--output", help="Output folder where to put the netCDFs.", default="."
    )
    argparser.add_argument(
        "--raw-output",
        help="Output folder where to put the non-merged netCDFs.",
    )
    argparser.add_argument(
        "-s",
        "--skip-existing",
        help="Do not overwrite existing files.",
        action="store_true",
    )
    argparser.add_argument(
        "metafile",
        help="Path to the MDB file containing the station and variable information.",
    )
    argparser.add_argument(
        "folder", help="Path to a folder including many *.mdb files to convert."
    )
    args = argparser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    outputs = convert_melcc_obs(
        args.metafile,
        args.folder,
        output=args.raw_output or args.output,
        overwrite=not args.skip_existing,
    )

    if args.concat:
        new_vars = {}
        for (new_var, table), path in outputs.items():
            new_vars.setdefault(new_var, []).append(path)

        for new_var, fichiers in new_vars.items():
            try:
                concat(fichiers, args.output, overwrite=not args.skip_existing)
            except Exception as err:
                logger.error(err)
