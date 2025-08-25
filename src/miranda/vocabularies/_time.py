__all__ = ["xarray_frequencies_to_cmip6like"]

# Manually map xarray frequencies to CMIP6/CMIP5 controlled vocabulary.
# see: https://github.com/ES-DOC/pyessv-archive
xarray_frequencies_to_cmip6like = {
    "h": "hr",
    "H": "hr",
    "D": "day",
    "W": "sem",
    "M": "mon",
    "Q": "qtr",  # TODO does this make sense? does not exist in cmip6 CV
    "A": "yr",
    "Y": "yr",
}
