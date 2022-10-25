import re

import pandas as pd
import pint
import xarray as xr
from pint import Unit
from xclim.core import calendar

KiB = int(pow(2, 10))
MiB = int(pow(2, 20))
GiB = int(pow(2, 30))

u = pint.UnitRegistry(autoconvert_offset_to_baseunit=True, on_redefinition="ignore")

# General purpose units
null = pint.Context("none")
u.add_context(null)
u.define("fraction = []")
u.define("percent = 1e-2 fraction = pct")
u.define("degC = kelvin; offset: 273.15 = celsius = C")
u.define("d = day")
u.define("h = hour")  # Not the Planck constant...
u.define("[speed] = [length] / [time]")

u.define("[precipitation] = [mass] / [length] ** 2 / [time]")
u.define("mmday = 1 kg / meter ** 2 / day")
u.define("[discharge] = [length] ** 3 / [time]")
u.define("cms = meter ** 3 / second")

# Allows transformation from mm day to kg / m2 / s
hq = pint.Context("hq")

hq.add_transformation(
    "[mass] / [length]**2",
    "[length]",
    lambda ureg, x: x / (1000 * ureg.kg / ureg.m**3),
)
hq.add_transformation(
    "[mass] / [length]**2 / [time]",
    "[length] / [time]",
    lambda ureg, x: x / (1000 * ureg.kg / ureg.m**3),
)
hq.add_transformation(
    "[length] / [time]",
    "[mass] / [length]**2 / [time]",
    lambda ureg, x: x * (1000 * ureg.kg / ureg.m**3),
)
u.add_context(hq)
u.enable_contexts(hq)


def units2pint(value: str) -> Unit:
    """Return the pint Unit for the DataArray units.

    Parameters
    ----------
    value : str
      Unit expression.

    Returns
    -------
    pint.Unit
      Pint compatible units.

    """

    def _transform(s):
        """Convert a CF-unit string to a pint expression."""
        return re.subn(r"\^?(-?\d)", r"**\g<1>", s)[0]

    value = value.replace("%", "pct")
    try:  # Pint compatible
        return u.parse_expression(value).units
    except (
        pint.UndefinedUnitError,
        pint.DimensionalityError,
    ):  # Convert from CF-units to pint-compatible
        return u.parse_expression(_transform(value)).units


def get_time_frequency(d: xr.Dataset):
    """Try to understand the datasets frequency.

    If it can't be inferred with :py:func:`xarray.infer_freq` it tries to:
    - look for a "freq" attrs in the global or time variable attributes.
    - infer monthly frequency if all time steps are between 27 and 32 days

    returns the offset a list of (multiplicator, base) and it's meaning (single word)
    """
    freq = xr.infer_freq(d.time)

    # Hacky workaround for irregular Monthly data
    if freq is None or (
        1 < int(calendar.parse_offset(freq)[0]) < 32 and freq.endswith("D")
    ):
        if 'freq' in d.attrs:
            freq = d.attrs['freq']
        elif 'freq' in d.time.attrs:
            freq = d.time.attrs['freq']
        elif (
            (d.diff("time") < pd.Timedelta(32, "D"))
            & (d.diff("time") > pd.Timedelta(27, "D"))
        ).all():
            freq = "1M"
        else:
            raise TypeError()

    offset = [int(calendar.parse_offset(freq)[0]), calendar.parse_offset(freq)[1]]

    time_units = {
        "s": "second",
        "T": "minute",
        "h": "hour",
        "D": "day",
        "M": "month",
        "W": "week",
        "A": "year",
    }
    if offset[1] in ["S", "H"]:
        offset[1] = offset[1].lower()
    offset_meaning = time_units[offset[1]]
    return offset, offset_meaning
