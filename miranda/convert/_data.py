__all__ = [
    "nasa_ag_variables",
    "nrcan_variables",
    "project_institutes",
    "sc_earth_variables",
    "wfdei_gem_capa_variables",
    "xarray_frequencies_to_cmip6",
]

nrcan_variables = ["tasmin", "tasmax", "pr"]
nasa_ag_variables = ["prate", "rhstmax", "srad", "tavg", "tmax", "tmin", "wndpsd"]
wfdei_gem_capa_variables = ["huss", "pr", "ps", "rlds", "rsds", "sfcWind", "tas"]
sc_earth_variables = ["prcp", "tdew", "tmean", "trange", "wind"]


project_institutes = {
    "cfsr": "ncar",
    "era5": "ecmwf",
    "era5-single-levels": "ecmwf",
    "era5-land": "ecmwf",
    "merra2": "nasa",
    "nrcan-gridded-10km": "nrcan",
    "wfdei-gem-capa": "usask",
}

# map xarray freq to CMIP6 controlled vocabulary.
# see: https://github.com/WCRP-CMIP/CMIP6_CVs/blob/master/CMIP6_frequency.json
xarray_frequencies_to_cmip6 = {
    "H": "hr",
    "D": "day",
    "M": "mon",
    "A": "yr",
    "Q": "qtr",  # TODO does this make sense? does not exist in cmip6 CV
}
