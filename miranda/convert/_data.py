__all__ = [
    "era5_variables",
    "nasa_ag_variables",
    "nrcan_variables",
    "project_institutes",
    "sc_earth_variables",
    "wfdei_gem_capa_variables",
    "xarray_frequencies_to_cmip6like",
]


era5_variables = [
    "d2m",
    "pev",
    "sde",
    "sd",
    "sf",
    "t2m",
    "tp",
    "u10",
    "v10",
]
nrcan_variables = ["tasmin", "tasmax", "pr"]
nasa_ag_variables = ["prate", "rhstmax", "srad", "tavg", "tmax", "tmin", "wndpsd"]
sc_earth_variables = ["prcp", "tdew", "tmean", "trange", "wind"]
wfdei_gem_capa_variables = ["huss", "pr", "ps", "rlds", "rsds", "sfcWind", "tas"]


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
xarray_frequencies_to_cmip6like = {
    "H": "hr",
    "D": "day",
    "M": "mon",
    "A": "yr",
    "Q": "qtr",  # TODO does this make sense? does not exist in cmip6 CV
}
