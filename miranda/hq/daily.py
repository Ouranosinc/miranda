import csv
import datetime as dt
import json
import re
from pathlib import Path
from typing import Optional
from typing import Tuple
from typing import Union

import numpy as np
import pandas as pd
import pint
import xarray as xr


def units():
    """Return a unit environment."""
    u = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)
    null = pint.Context("none")
    u.add_context(null)

    u.define("fraction = []")
    u.define("percent = 1e-2 fraction = pct")
    u.define("degC = kelvin; offset: 273.15 = celsius = C")
    u.define("d = day")
    u.define("h = hour")  # Not the Planck constant...
    u.define("[speed] = [length] / [time]")

    u.define("[precipitation] = [mass] / [length] ** 2 / [time]")
    u.define("mmday = 1000 kg / meter ** 2 / day")

    u.define("[discharge] = [length] ** 3 / [time]")
    u.define("cms = meter ** 3 / second")

    hydro = pint.Context("hydro")
    hydro.add_transformation(
        "[mass] / [length]**2",
        "[length]",
        lambda ureg, x: x / (1000 * ureg.kg / ureg.m ** 3),
    )
    hydro.add_transformation(
        "[mass] / [length]**2 / [time]",
        "[length] / [time]",
        lambda ureg, x: x / (1000 * ureg.kg / ureg.m ** 3),
    )
    hydro.add_transformation(
        "[length] / [time]",
        "[mass] / [length]**2 / [time]",
        lambda ureg, x: x * (1000 * ureg.kg / ureg.m ** 3),
    )
    u.add_context(hydro)
    u.enable_contexts(hydro)
    return u


unites = units()

# CMOR-like attributes
cmor = json.load(open(Path(__file__).parent / "cf_attrs.json"))["variable_entry"]

fp = r"[-+]?\d*,\d+|\d+"

section_patterns = r"(\w+) :\n"

meta_patterns = {
    "Installation": {
        "nom": "Nom;(.*)",
        "type": "Type;(.*)",
        "ouverture": "Ouverture;(.+)",
        "fermeture": "Fermeture;(.+)",
        "x": rf"XCOORD \(degrés\.décimales\);({fp})",
        "y": rf"YCOORD \(degrés\.décimales\);({fp})",
        "z": rf"ZCOORD \(mètres\);({fp})",
    },
    "Données": {
        "variable": "Type de donnée;(.*)",
        "fuseau": "Fuseau horaire;(.*)",
        "pas": "Pas de temps;(.*)",
        "mesure": "Type de mesure;(.*)",
        "unité": "Unité;(.*)",
    },
}

csv_patterns = ("Propriétaires", "Correspondances")

data_header_pattern = "Dateheure;Valeur;Qualite;Statut"

converters = {
    "ouverture": np.datetime64,
    "fermeture": np.datetime64,
    "x": lambda x: float(x.replace(",", ".")),
    "y": lambda x: float(x.replace(",", ".")),
    "z": lambda x: float(x.replace(",", ".")),
    "Dateheure": lambda x: dt.datetime.strptime(x, "%Y-%m-%d %H:%M"),
    "Valeur": lambda x: float(x.replace(",", ".")) if x == "" else np.nan,
    "Qualite": lambda x: int(x.split("-")[0]),
    "Statut": lambda x: int(x.split("-")[0]),
    "unité": lambda x: cf_units.get(x, x),
    "pas": lambda x: cf_frequency.get(x, x),
}


def guess_variable(meta, cf_table: Optional[dict]) -> str:
    """Return the corresponding CMOR variable."""
    if cf_table is None:
        cf_table = cmor

    v = meta["variable"]

    corr = {
        "Pluie (mm)": "prlp",
        "Neige (mm)": "prsn",
        "Épaisseur de neige": "snd",
        "Vitesse du vent 10 mètres": "sfcWind",
        "Direction du vent 10 mètres": "sfcWindAz",
        "Humidité relative 2 mètres": "hurs",
    }

    name = str()
    if v in corr:
        name = corr[v]
    else:
        if v == "Température":
            if meta["mesure"] == "Maximum":
                name = "tasmax"
            elif meta["mesure"] == "Minimum":
                name = "tasmin"

    if meta["pas"] != cf_table[name]["frequency"]:
        raise ValueError("Unexpected frequency.")

    return name


cf_units = {"°C": "celsius", "mm": "mmday"}
cf_frequency = {"Fin du pas journalier": "day", "Instantanée du pas horaire": "1h"}
cf_attrs_names = {"x": "lon", "y": "lat", "z": "elevation", "nom": "site"}


def extract_daily(path) -> Tuple[dict, pd.DataFrame]:
    """Extract data and metadata from HQ meteo file."""

    with open(path, encoding="latin1") as fh:
        txt = fh.read()
        meta, data = re.split(data_header_pattern, txt, maxsplit=2)

    sections = iter(re.split(section_patterns, meta)[1:])

    m = dict()
    for sec in sections:
        if sec in meta_patterns:
            content = next(sections)
            for key, pat in meta_patterns[sec].items():
                match = re.search(pat, content)
                if match:
                    m[key] = match.groups()[0]

        elif sec in csv_patterns:
            content = next(sections).strip()
            m[sec] = list(
                csv.reader(content.splitlines(), dialect="unix", delimiter=";")
            )

    d = pd.read_csv(
        path,
        delimiter=";",
        skiprows=len(meta.splitlines()),
        encoding="latin1",
        index_col=0,
        parse_dates=True,
        infer_datetime_format=True,
        decimal=",",
    )

    return m, d


def to_cf(
    meta: dict, data: pd.DataFrame, cf_table: Optional[dict] = None
) -> xr.DataArray:
    """Return CF-compliant metadata."""

    if dict is None:
        cf_table = cmor

    # Convert meta values
    m = dict()
    for key, val in meta.items():
        m[key] = converters.get(key, lambda x: x)(val)

    # Get default variable attributes
    name = guess_variable(m, cf_table)
    attrs = cf_table[name]

    # Add custom HQ attributes
    for key, val in cf_attrs_names.items():
        if key in m:
            attrs[val] = m[key]

    x = data["Valeur"].values.astype(float)

    # Convert units
    if attrs["units"] != m["unité"]:
        x = unites.convert(x, m["unité"], units2pint(attrs["units"]))

    coords = {k: attrs.pop(k, np.nan) for k in ["lon", "lat", "elevation", "site"]}
    coords["time"] = data.index.values
    cf_corrected = xr.DataArray(
        data=x, dims="time", coords=coords, name=name, attrs=attrs
    )
    return cf_corrected


def open_csv(path: Union[str, Path], cf_table: Optional[dict] = None) -> xr.DataArray:
    """Extract daily HQ meteo data and convert to xr.DataArray with CF-Convention attributes."""
    meta, data = extract_daily(path)
    return to_cf(meta, data, cf_table)


def units2pint(value: str) -> pint.Quantity:
    """Return the pint Unit for the DataArray units.

    Parameters
    ----------
    value : str
      Unit expression.

    Returns
    -------
    pint.Quantity
      Pint compatible units.

    """

    def _transform(s):
        """Convert a CF-unit string to a pint expression."""
        return re.subn(r"\^?(-?\d)", r"**\g<1>", s)[0]

    value = value.replace("%", "pct")
    try:  # Pint compatible
        return unites.parse_expression(value).units
    except (
        pint.UndefinedUnitError,
        pint.DimensionalityError,
    ):  # Convert from CF-units to pint-compatible
        return unites.parse_expression(_transform(value)).units


# if __name__ == "__main__":
#     fns = list(Path("/home/david/projects/ForestFires/data/HQ/LG_Tmin").glob("*.csv"))[
#         :4
#     ]
#     das = [open_csv(fn) for fn in fns]
#     da = xr.concat(das, dim="site")
