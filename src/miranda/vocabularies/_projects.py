__all__ = ["project_institute"]


institutes = {
    "ncar": ["cfsr"],
    "ecmwf": [
        "era5",
        "era5-land",
        "era5-land-monthly-means",
        "era5-monthly",
        "era5-pressure-levels",
        "era5-pressure-levels-preliminary-back-extension",
        "era5-pressure-monthly-means-levels-preliminary-back-extension",
        "era5-single-levels",
        "era5-single-levels-monthly-means",
        "era5-single-levels-monthly-means-preliminary-back-extension",
        "era5-single-levels-preliminary-back-extension",
    ],
    "nasa": ["merra2", "NEX-GDDP-CMIP6"],
    "nrcan": ["nrcan-gridded-10km"],
    "usask": ["wfdei-gem-capa"],
    "eccc": ["rdrs-v21"],
}


def project_institute(project: str, raise_on_error: bool = False) -> str:
    for institute, projects in institutes.items():
        if project in projects:
            return institute
    if raise_on_error:
        raise ValueError("Project does not have associated institute.")
    return "Unknown"
