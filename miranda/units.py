import pint
import re

u = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)

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
    lambda ureg, x: x / (1000 * ureg.kg / ureg.m ** 3),
)
hq.add_transformation(
    "[mass] / [length]**2 / [time]",
    "[length] / [time]",
    lambda ureg, x: x / (1000 * ureg.kg / ureg.m ** 3),
)
hq.add_transformation(
    "[length] / [time]",
    "[mass] / [length]**2 / [time]",
    lambda ureg, x: x * (1000 * ureg.kg / ureg.m ** 3),
)
u.add_context(hq)
u.enable_contexts(hq)


def units2pint(value: str) -> u.Unit:
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
        return u.parse_expression(value).units
    except (
        pint.UndefinedUnitError,
        pint.DimensionalityError,
    ):  # Convert from CF-units to pint-compatible
        return u.parse_expression(_transform(value)).units
