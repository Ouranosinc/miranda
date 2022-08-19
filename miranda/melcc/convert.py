import datetime as dt
import json
import logging
import re
from argparse import ArgumentParser
from copy import deepcopy
from functools import reduce
from io import StringIO
from pathlib import Path
from subprocess import run

import numpy as np
import pandas as pd
import xarray as xr
from pint.errors import DimensionalityError
from xclim.core.units import amount2rate, convert_units_to, pint_multiply, str2pint

logger = logging.getLogger("miranda.melcc")


frequencies = {
    "A": ("min", "1 min"),  # T
    "B": ("point", None),  # "Élémentaire"...
    "C": ("12hr", "12 h"),  # 12H
    "F": ("point", None),  # Instantaneous??
    "H": ("hr", "1 hr"),  # H
    "Q": ("day", "1 d"),  # D
    "Z": ("4min", "4 min"),  # 4T
}


with (Path(__file__).parent / "attrs.json").open() as f:
    melcc_attrs = json.load(f)


def parse_var_code(vcode, definitions=None):
    vcode = vcode[:2] + vcode[2:-1].lower() + vcode[-1]
    match = re.match(r"([^\d]*)(\d*)([ABCFHQZ])", vcode)
    if match is None:
        raise ValueError(f"Failed to parse variable {vcode}.")
    short_code, instrument, freq = match.groups()
    meta = {
        "freq": frequencies[freq][0],
        "sampling": frequencies[freq][1],
        "var_name": vcode[:-1],
    }
    if vcode[:-1] in melcc_attrs["variable_entry"]:
        attrs = deepcopy(melcc_attrs["variable_entry"][vcode[:-1]])
    else:
        attrs = deepcopy(melcc_attrs["variable_entry"][short_code])
    meta["attrs"] = {k: v for k, v in attrs.items() if not k.startswith("_")}
    meta["attrs"].update(instrument=np.int32(instrument), melcc_code=vcode)
    meta.update({k[1:]: v for k, v in attrs.items() if k.startswith("_")})
    if definitions is not None and vcode.lower() in definitions:
        meta["attrs"]["melcc_description"] = definitions.loc[vcode]
    return meta


def list_tables(dbfile):
    """List the tables of a MDB file."""
    res = run(["mdb-tables", dbfile], capture_output=True, encoding="utf-8")
    if res.returncode != 0:
        raise ValueError(
            f"Calling mdb-tables on {dbfile} failed with code {res.returncode}: {res.stderr}"
        )
    return res.stdout.upper().strip().split()


def read_table(dbfile, table, **kwargs):
    res = run(
        ["mdb-export", "-T", "%Y-%m-%dT%H:%M:%S", dbfile, table],
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
        .rename(columns={"NO_SEQ_STATION": "station", "DATE": "time"})
        .drop(columns=["CODE_STATUT_DONNEE", "STATUT_APPROBATION"])
        .set_index(["station", "time"])
        .VALEUR_DONNEE
    )
    NaT = df.index.get_level_values("time").isnull().sum()
    if NaT > 0:
        logger.warning(
            f"{NaT} values dropped because of unparsable timestamps (probably data from before 1900)."
        )
        df = df[df.index.get_level_values("time").notnull()]
    return df[~df.index.duplicated()].to_xarray()


def convert_melcc(da, vcode, stations=None, definitions=None):
    data = parse_var_code(vcode, definitions)
    da.attrs.update(data["attrs"])
    if "scale" in data:
        da = pint_multiply(da, str2pint(data["scale"]))
    if data.get("to_rate"):
        da = pint_multiply(da, 1 / str2pint(data["sampling"]))
    if "convert_to" in data:
        try:
            da = convert_units_to(da, data["convert_to"])
        except DimensionalityError:
            logger.warn(f"Converting {da.name} from amount to rate.")
            da = convert_units_to(amount2rate(da), data["convert_to"])
    if "nc_name" in data:
        da = da.rename(data["nc_name"])
    else:
        da = da.rename(data["var_name"])
    if stations is not None:
        try:
            stat = stations.sel(station=da.station)
        except KeyError:
            extra_stations = set(da.station.values) - set(stations.station.values)
            logger.error(
                f"The following station number were not found in the station list : {sorted(list(extra_stations))}"
            )
            with open("unknown_stations.txt", "a") as f:
                f.write("\n".join(map(str, extra_stations)))
                f.write("\n")
            da = xr.merge([da, xr.Dataset().assign(station=stations)])[da.name].sel(
                station=da.station
            )
        else:
            da = da.assign_coords(station=stat)
    da.attrs.update(data["attrs"])
    return da, data


def read_stations(file):
    df = pd.read_excel(file, skiprows=4)
    df["OUVERTURE"] = pd.to_datetime(df.OUVERTURE, errors="coerce")
    df = df.rename(
        columns={
            "NO_STATION": "station_id",
            "NOM_STATION": "station_name",
            "LAT(°)": "lat",
            "LONG(°)": "lon",
            "ALT(m)": "elevation",
            "OUVERTURE": "station_opening",
            "FERMETURE": "station_closing",
            "TYPE_STATION": "station_type",
            "NO_SEQ_STATION": "station",
        }
    )
    ds = df.set_index("station").to_xarray()
    da = ds.set_coords(ds.data_vars.keys()).station
    da.lat.attrs.update(units="degree_north", standard_name="latitude")
    da.lon.attrs.update(units="degree_east", standard_name="longitude")
    da.elevation.attrs.update(units="m", standard_name="height")
    da["station_id"] = da["station_id"].astype(str)
    da["station_name"] = da["station_name"].astype(str)
    da["station_type"] = da["station_type"].astype(str)
    da.station_opening.attrs.update(description="Date of station creation.")
    da.station_closing.attrs.update(description="Date of station closure.")
    return da


def read_definitions(file):
    definitions = pd.read_excel(
        file, skiprows=4, names=["num", "type", "vcode", "description"]
    )
    definitions["vcode"] = definitions.vcode.str.lower()
    return definitions.set_index("vcode").description


def convert_snow_levels(file):
    # Stations
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
    stations.lat.attrs.update(units="degree_north", standard_name="latitude")
    stations.lon.attrs.update(units="degree_east", standard_name="longitude")
    stations.elevation.attrs.update(units="m", standard_name="height")
    stations.station_opening.attrs.update(description="Date of station creation.")
    stations.station_closing.attrs.update(description="Date of station closure.")

    # Periods
    periods = pd.read_excel(
        file, sheet_name="Périodes standards", names=["start", "end", "middle"]
    )
    periods = (
        periods.to_xarray()
        .to_array()
        .rename(variable="bnds", index="period")
        .drop_vars("bnds")
        .rename("period_bnds")
    )

    # Data
    data = pd.read_excel(
        file,
        sheet_name="Données",
        names=[
            "station",
            "time",
            "snd",
            "snd_flag",
            "snow_density",
            "snow_density_flag",
            "snw",
            "snw_flag",
        ],
    )
    ds = data.set_index(["station", "time"]).to_xarray()
    ds["station"] = stations.sel(station=ds.station)


if __name__ == "__main__":
    argparser = ArgumentParser(
        description="Convert MS Access data to netCDF. Requires mdbtools to be installed. By default, a single netCDF is created for each data table."
    )
    argparser.add_argument(
        "-v", "--verbose", help="Increase verbosity of the script.", action="store_true"
    )
    argparser.add_argument(
        "-c",
        "--concat",
        help="Merge all variables with the same frequency in individual datasets.",
        action="store_true",
    )
    argparser.add_argument(
        "-o", "--output", help="Output folder where to put the netCDFs.", default="."
    )
    argparser.add_argument(
        "-d",
        "--definitions",
        help="Excel file containing the variable code definitions. Optional if --concat is not passed, unused otherwise.",
    )
    argparser.add_argument(
        "stations", help="Path to the excel file including the station information."
    )
    argparser.add_argument(
        "folder", help="Path to a folder including many *.mdb files to convert."
    )
    args = argparser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    definitions = None
    if args.definitions is not None:
        # a pd.Series  : index is variable code and value is french description
        definitions = read_definitions(args.definitions)

    output = Path(args.output)

    stations = read_stations(args.stations)

    logger.info("Parsing databases.")
    folder = Path(args.folder)
    dss = {}
    for database in folder.glob("*.mdb"):
        tables = list_tables(str(database))
        for table in tables:
            logger.info(f"Parsing {database}:{table}.")
            raw = read_table(database, table).rename(table)
            da, data = convert_melcc(
                raw, table, stations=stations, definitions=definitions
            )
            if args.concat:
                dss.setdefault(data["freq"], {}).setdefault(da.name, []).append(da)
            else:
                ds = da.to_dataset()
                ds.attrs.update(melcc_attrs["header"])
                ds.attrs[
                    "history"
                ] = f"[{dt.datetime.now():%Y-%m-%d %H:%M:%S}] Conversion from {database.name}:{table} to netCDF."
                da.to_netcdf(
                    output / f"MELCC_{data['freq']}_{da.name}_{da.melcc_code}.nc"
                )

    if args.concat:
        logger.info("Concatening all variables.")
        for freq, variables in dss.items():
            ds = xr.Dataset()
            for var_name, das in variables.items():
                if len(das) > 1:
                    n_conflicts = reduce(np.logical_and, das).sum().item()
                    if n_conflicts > 0:
                        logger.warning(
                            f"{n_conflicts} samples have multiple measurements for {var_name}/{freq}."
                        )
                    da = xr.merge(
                        sorted(das, key=lambda da: da.instrument),
                        compat="override",
                        combine_attrs="drop_conflicts",
                    )[var_name]
                else:
                    da = das[0]
                ds = ds.assign({var_name: da})
            ds.attrs.update(melcc_attrs["header"])
            ds.attrs[
                "history"
            ] = f"[{dt.datetime.now():%Y-%m-%d %H:%M:%S}] Conversion from MS Access files to netCDF."
            ds.to_netcdf(output / f"MELCC_{data['freq']}.nc")
