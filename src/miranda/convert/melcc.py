"""MELCC (Québec) Weather Stations data conversion module."""

from __future__ import annotations
import datetime as dt
import logging
import os
import re
import subprocess  # noqa: S404
import warnings
from argparse import ArgumentParser
from collections import defaultdict
from collections.abc import Sequence
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray
import xarray as xr
from xclim.core.formatting import update_history
from xclim.core.units import convert_units_to, pint_multiply, str2pint

from miranda import __version__
from miranda.convert.corrections import dataset_corrections
from miranda.treatments import metadata_conversion
from miranda.treatments.utils import load_json_data_mappings


logger = logging.getLogger("miranda.convert.melcc")


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
    "concat",
    "convert_mdb",
    "convert_melcc_obs",
    "convert_snow_table",
    "list_tables",
    "parse_var_code",
    "read_definitions",
    "read_stations",
    "read_table",
]


def parse_var_code(vcode: str) -> dict[str, Any]:
    """
    Parse variable code to generate metadata.

    Parameters
    ----------
    vcode : str
        The variable code.

    Returns
    -------
    dict[str, Any]
        The metadata dictionary.
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


def _validate_db_file(db_file: str) -> str:
    """
    Validate the database file and ensure that input is trustworthy.

    Parameters
    ----------
    db_file : str
        The database file.

    Returns
    -------
    str
        The database file.
    """
    if len(db_file) > 1:
        raise ValueError("Only one database file can be processed at a time.")
    if not Path(db_file).is_file():
        raise FileNotFoundError(f"File {db_file} not found.")
    return db_file


def list_tables(db_file: str | os.PathLike[str]) -> list[str]:
    """
    List the tables of an MDB file.

    Parameters
    ----------
    db_file : str or os.PathLike
        The database file.

    Returns
    -------
    list of str
        The list of tables.
    """
    try:
        res = subprocess.run(  # noqa: S603
            ["mdb-tables", _validate_db_file(db_file)],
            capture_output=True,
            encoding="utf-8",
        )
    except subprocess.CalledProcessError as err:
        raise ValueError(f"Calling mdb-tables on {db_file} failed with code {err.returncode}: {err.stderr}") from err
    if res.returncode != 0:
        raise ValueError(f"Calling mdb-tables on {db_file} failed with code {res.returncode}: {res.stderr}")
    return res.stdout.lower().strip().split()


def read_table(db_file: str | os.PathLike[str], table: str | os.PathLike) -> xarray.Dataset:
    """
    Read a MySQL table into an xarray object.

    Parameters
    ----------
    db_file : str or os.PathLike
        The database file.
    table : str or os.PathLike
        The table to read.

    Returns
    -------
    xarray.Dataset
        An xarray Dataset with the table data.
    """
    try:
        res = subprocess.run(  # noqa: S603
            [
                "mdb-export",
                "-T",
                "%Y-%m-%dT%H:%M:%S",
                _validate_db_file(db_file),
                table.upper(),
            ],
            capture_output=True,
            encoding="utf-8",
            check=True,
        )
    except subprocess.CalledProcessError as err:
        msg = f"Calling mdb-export on {db_file} failed with code {err.returncode}: {err.stderr}"
        raise ValueError(msg) from err
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
    nat = df.index.get_level_values("time").isnull().sum()
    if nat > 0:
        msg = f"{nat} values dropped because of unparsable timestamps (probably data from before 1900)."

        logger.warning(msg)
        df = df[df.index.get_level_values("time").notnull()]
    return df[~df.index.duplicated()].to_xarray()


def read_stations(db_file: str | os.PathLike) -> pd.DataFrame:
    """
    Read station file using mdbtools.

    Parameters
    ----------
    db_file : str or os.PathLike
        The database file.

    Returns
    -------
    pandas.DataFrame
        A Pandas DataFrame with the station information.
    """
    try:
        res = subprocess.run(  # noqa: S603
            [
                "mdb-export",
                "-T",
                "%Y-%m-%dT%H:%M:%S",
                _validate_db_file(db_file),
                "STATIONs",
            ],
            capture_output=True,
            encoding="utf-8",
            check=True,
        )
    except subprocess.CalledProcessError as err:
        msg = f"Calling mdb-export on {db_file} failed with code {err.returncode}: {err.stderr}"
        raise ValueError(msg) from err

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


def read_definitions(db_file: str) -> pd.DataFrame:
    """
    Read variable definition file using mdbtools.

    Parameters
    ----------
    db_file : str
        The database file.

    Returns
    -------
    pandas.DataFrame
        The variable definitions.
    """
    try:
        res = subprocess.run(  # noqa: S603
            ["mdb-export", _validate_db_file(db_file), "DDs"],
            capture_output=True,
            encoding="utf-8",
            check=True,
        )
    except subprocess.CalledProcessError as err:
        msg = f"Calling mdb-export on {db_file} failed with code {err.returncode}: {err.stderr}"
        raise ValueError(msg) from err
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
    """
    Convert microsoft databases of MELCCFP observation data to xarray objects.

    Parameters
    ----------
    database : str or Path
        The database file.
    stations : xr.Dataset
        The station list.
    definitions : xr.Dataset
        The variable definitions.
    output : str or Path
        The output folder.
    overwrite : bool
        Whether to overwrite existing files. Default: True.

    Returns
    -------
    dict[tuple[str, str], Path]
        The converted files.
    """
    outs = {}
    tables = list_tables(database)
    for table in tables:
        if table.startswith("gdb") or table.startswith("~"):
            continue
        msg = f"Parsing {database}:{table}."
        logger.info(msg)
        meta = parse_var_code(table)
        code = meta["melcc_code"]
        existing = list(output.glob(f"*_{code}_MELCC_{meta['freq']}_*.nc"))
        if existing and not overwrite:
            if len(existing) > 1:
                raise ValueError(f"Found more than one existing file for table {code}!")
            file = existing[0]
            vname = file.stem.split(code)[0][:-1]
            msg = f"File already exists {file}, skipping."
            logger.info(msg)
            outs[(vname, code)] = file
            continue
        raw = read_table(database, table)
        if np.prod(list(raw.dims.values())) == 0:
            # If any dimension is 0.
            logger.warning("The table is empty.")
            continue
        vv = meta.pop("var_name")

        raw = raw.rename({table: vv, f"{table}_flag": f"{vv}_flag"})
        raw.attrs["frequency"] = meta.pop("freq")
        raw[vv].attrs.update(**meta)

        if table.lower() not in definitions.index:
            warnings.warn(f"The {code} variable wasn't defined in the definition table.", stacklevel=2)
        else:
            dd = definitions.loc[table]
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
            msg = f"The following station number were not found in the station list : {sorted(list(extra_stations))}"

            logger.error(msg)
            raw = xr.merge([raw, xr.Dataset().assign(station=stations)]).sel(station=raw.station)
        else:
            raw = raw.assign_coords(station=stat)

        ds = dataset_corrections(raw, "melcc-obs")
        new_var_name = list(filter(lambda k: not k.endswith("_flag"), ds.data_vars.keys()))[0]
        ds.attrs["history"] = f"[{dt.datetime.now():%Y-%m-%d %H:%M:%S}] Conversion from {database.name}:{table} to netCDF."
        date = "-".join(ds.indexes["time"][[0, -1]].strftime("%Y%m"))
        outs[(new_var_name, code)] = output / f"{new_var_name}_{code}_MELCC_{raw.attrs['frequency']}_{date}.nc"
        ds.to_netcdf(outs[(new_var_name, code)])
    return outs


def convert_melcc_obs(
    metafile: str | Path,
    folder: str | Path,
    output: str | Path | None = None,
    overwrite: bool = True,
) -> dict[tuple[str, str], Path]:
    """
    Convert MELCCFP observation data to xarray data objects, returning paths.

    Parameters
    ----------
    metafile : str or Path
        The metadata file.
    folder : str or Path
        The folder containing the MDB files.
    output : str or Path, optional
        The output folder. Default: None.
    overwrite : bool
        Whether to overwrite existing files. Default: True.

    Returns
    -------
    dict[(str, str), Path]
        The converted files.
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
    files: Sequence[str | os.PathLike[str]],
    output_folder: str | os.PathLike[str],
    overwrite: bool = True,
) -> Path:
    """
    Concatenate converted weather station files.

    Parameters
    ----------
    files : Sequence of str or os.PathLike
        The files to concatenate.
    output_folder : str or os.PathLike
        The output folder.
    overwrite : bool
        Whether to overwrite existing files. Default: True.

    Returns
    -------
    Path
        The output path.
    """
    *vv, _, melcc, freq, _ = Path(files[0]).stem.split("_")
    vv = "_".join(vv)
    msg = f"Concatenating variables from {len(files)} files ({vv})."
    logger.info(msg)
    # Magic one-liner to parse all date_start and date_end entries from the file names.
    dates_start, dates_end = list(
        zip(
            *[map(int, Path(file).stem.split("_")[-1].split("-")) for file in files],
            strict=False,
        )
    )
    outpath = Path(output_folder) / f"{vv}_{melcc}_{freq}_{min(dates_start):06d}-{max(dates_end):06d}.nc"
    if outpath.is_file() and not overwrite:
        msg = f"Already done in {outpath}. Skipping."
        logger.info(msg)
        return outpath

    dss = dict()
    for file in files:
        ds = xr.open_dataset(file, chunks={"time": 1000})
        for crd in ds.coords.values():
            crd.load()
        if list(sorted(ds.data_vars.keys())) != [vv, f"{vv}_flag"]:
            raise ValueError(f"Unexpected variables in {file}. Got {ds.data_vars.keys()}, expected {vv} and {vv}_flag.")

        priority = ds[vv].attrs["priority"]
        if priority in dss:
            raise ValueError(f"Variable {vv} of {file} has the same priority ({priority}) than another file of the same file list.")
        dss[priority] = ds

    ds_all = xr.merge([ds.coords.to_dataset() for ds in dss.values()], combine_attrs="drop_conflicts")
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
    ds_merged[vv].attrs["melcc_code"] = "Merged sources in ascending priority : " + " ,".join(map(str, instruments))

    ds_merged.attrs.update(
        source="info-climat-merged",
        title="Station observations of the MELCC - all instruments merged.",
        title_fr="Observations aux stations du MELCC - instruments fusionnés",
    )
    ds_merged.to_netcdf(outpath)
    return outpath


def convert_snow_table(file: str | os.PathLike[str] | Path, output: str | os.PathLike[str] | Path) -> None:
    """
    Convert snow data given through an Excel file.

    This private data is not included in the MDB files.

    Parameters
    ----------
    file : str or os.PathLike or Path
        The Excel file with sheets: "Stations", "Périodes standards", and "Données".
    output : str or os.PathLike or Path
        Folder where to put the netCDF files (one for each of snd, sd and snw).
    """
    logger.info("Parsing stations.")
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
    logger.info("Parsing observation periods.")
    periods = pd.read_excel(file, sheet_name="Périodes standards", names=["start", "end", "middle"])
    periods = periods[["start", "end"]].to_xarray()
    periods = periods.to_array().rename(variable="bnds", index="num_period").drop_vars("bnds").rename("period_bnds")
    periods.attrs.update(
        description=(
            "Bounds of the sampling periods of the MELCC. "
            "Observations are taken manually once per period. "
            "The year of these bounds should be ignored."
        )
    )

    # Data
    logger.info("Parsing data.")
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
    meta = load_json_data_mappings("melcc")
    ds = metadata_conversion(ds, "melcc-snow", meta)
    date = "-".join(ds.indexes["time"][[0, -1]].strftime("%Y%m"))
    # Save
    logger.info("Saving to files.")
    for vv in ["sd", "snd", "snw"]:
        ds[[vv, f"{vv}_flag"]].to_netcdf(output / f"{vv}_{ds[vv].melcc_code}_MELCC_2sem_{date}.nc")


if __name__ == "__main__":
    argparser = ArgumentParser(
        description="Convert MS Access data to netCDF. Requires mdbtools to be installed. By default, a single netCDF is created for each data table."
    )
    argparser.add_argument("-v", "--verbose", help="Increase verbosity of the script.", action="store_true")
    argparser.add_argument(
        "-c",
        "--concat",
        help="Merge all variables with the same frequency in individual datasets, according to their priority.",
        action="store_true",
    )
    argparser.add_argument("-o", "--output", help="Output folder where to put the netCDFs.", default=".")
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
    argparser.add_argument("folder", help="Path to a folder including many *.mdb files to convert.")
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
        new_vars = defaultdict(list)
        for (new_var, _), path in outputs.items():
            new_vars[new_var].append(path)

        for fichiers in new_vars.values():
            try:
                concat(fichiers, args.output, overwrite=not args.skip_existing)
            except Exception as err:  # noqa: PERF203,BLE001
                logger.error(err)
