import json
import logging.config
import os
from pathlib import Path
from typing import List, Mapping, Union

from miranda.scripting import LOGGING_CONFIG
from miranda.storage import report_file_size

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "era5_variables",
    "gather_agcfsr",
    "gather_agmerra",
    "gather_era5_land",
    "gather_era5_pressure_levels",
    "gather_era5_single_levels",
    "gather_nrcan_gridded_obs",
    "gather_sc_earth",
    "gather_wfdei_gem_capa",
    "nasa_ag_variables",
    "nrcan_variables",
    "reanalysis_project_institutes",
    "sc_earth_variables",
    "wfdei_gem_capa_variables",
    "xarray_frequencies_to_cmip6like",
]

data_folder = Path(__file__).parent / "data"
era5_variables = json.load(open(data_folder / "ecmwf_cf_attrs.json"))[
    "variable_entry"
].keys()
nrcan_variables = ["tasmin", "tasmax", "pr"]
nasa_ag_variables = json.load(open(data_folder / "nasa_cf_attrs.json"))[
    "variable_entry"
].keys()
sc_earth_variables = ["prcp", "tdew", "tmean", "trange", "wind"]
wfdei_gem_capa_variables = json.load(open(data_folder / "usask_cf_attrs.json"))[
    "variable_entry"
].keys()

reanalysis_project_institutes = {
    "cfsr": "ncar",
    "era5": "ecmwf",
    "era5-land": "ecmwf",
    "era5-land-monthly-means": "ecmwf",
    "era5-monthly": "ecmwf",
    "era5-pressure-levels": "ecmwf",
    "era5-pressure-levels-preliminary-back-extension": "ecmwf",
    "era5-pressure-monthly-means-levels-preliminary-back-extension": "ecmwf",
    "era5-single-levels": "ecmwf",
    "era5-single-levels-monthly-means": "ecmwf",
    "era5-single-levels-monthly-means-preliminary-back-extension": "ecmwf",
    "era5-single-levels-preliminary-back-extension": "ecmwf",
    "merra2": "nasa",
    "nrcan-gridded-10km": "nrcan",
    "wfdei-gem-capa": "usask",
}


# Manually map xarray frequencies to CMIP6/CMIP5 controlled vocabulary.
# see: https://github.com/ES-DOC/pyessv-archive
xarray_frequencies_to_cmip6like = {
    "H": "hr",
    "D": "day",
    "W": "sem",
    "M": "mon",
    "Q": "qtr",  # TODO does this make sense? does not exist in cmip6 CV
    "A": "yr",
    "Y": "yr",
}


def _gather(
    name: str,
    variables: Mapping[str, List[str]],
    source: Union[str, os.PathLike],
    back_extension: bool,
    monthly_means: bool,
) -> Mapping[str, List[Path]]:
    source = Path(source)
    name = (
        f"{name}"
        f"{'-monthly-means' if monthly_means else ''}"
        f"{'-preliminary-back-extension' if back_extension else ''}"
    )
    logging.info(f"Gathering {name} files from: {source.as_posix()}")
    in_files = list()
    for v in variables:
        in_files.extend(list(sorted(source.rglob(f"{v}*{name}*.nc"))))
    logging.info(
        f"Found {len(in_files)} files, totalling {report_file_size(in_files)}."
    )
    return {name: in_files}


def gather_era5_single_levels(
    path: Union[str, os.PathLike],
    back_extension: bool = False,
    monthly_means: bool = False,
) -> Mapping[str, List[Path]]:
    # ERA5-Single-Levels source data
    return _gather(
        "era5-single-levels",
        era5_variables,
        source=path,
        back_extension=back_extension,
        monthly_means=monthly_means,
    )


def gather_era5_pressure_levels(
    path: Union[str, os.PathLike],
    back_extension: bool = False,
    monthly_means: bool = False,
) -> Mapping[str, List[Path]]:
    # ERA5-Pressure-Levels source data
    return _gather(
        "era5-pressure-levels",
        era5_variables,
        source=path,
        back_extension=back_extension,
        monthly_means=monthly_means,
    )


def gather_era5_land(path: Union[str, os.PathLike]) -> Mapping[str, List[Path]]:
    # ERA5-Land source data
    return _gather("era5-land", era5_variables, source=path, back_extension=False)


def gather_era5_land_sea_mask(path: Union[str, os.PathLike]) -> Mapping[str, Path]:
    try:
        land_sea_mask = dict(lsm=next(Path(path).glob("sftlf*era5*.nc")))
    except StopIteration:
        logging.error("No land_sea_mask found for ERA5.")
        raise FileNotFoundError()
    return land_sea_mask


def gather_agmerra(path: Union[str, os.PathLike]) -> Mapping[str, List[Path]]:
    # agMERRA source data
    source_agmerra = Path(path)
    logging.info(f"Gathering agMERRA from: {source_agmerra.as_posix()}")
    in_files_agmerra = list()
    for v in nasa_ag_variables:
        in_files_agmerra.extend(
            list(sorted(source_agmerra.rglob(f"AgMERRA_*_{v}.nc4")))
        )
    logging.info(
        f"Found {len(in_files_agmerra)} files, totalling {report_file_size(in_files_agmerra)}."
    )
    return dict(cfsr=in_files_agmerra)


def gather_agcfsr(path: Union[str, os.PathLike]) -> Mapping[str, List[Path]]:
    # agCFSR source data
    source_agcfsr = Path(path)
    logging.info(f"Gathering CFSR from: {source_agcfsr.as_posix()}")
    in_files_agcfsr = list()
    for v in nasa_ag_variables:
        in_files_agcfsr.extend(list(sorted(source_agcfsr.rglob(f"AgCFSR_*_{v}.nc4"))))
    logging.info(
        f"Found {len(in_files_agcfsr)} files, totalling {report_file_size(in_files_agcfsr)}."
    )
    return dict(cfsr=in_files_agcfsr)


def gather_nrcan_gridded_obs(path: Union[str, os.PathLike]) -> Mapping[str, List[Path]]:
    # NRCan Gridded Obs source data
    source_nrcan = Path(path)
    logging.info(f"Gathering NRCAN Gridded Obs from {source_nrcan.as_posix()}")
    in_files_nrcan = list()
    for v in nrcan_variables:
        in_files_nrcan.extend(list(sorted(source_nrcan.joinpath(v).glob(f"*{v}_*.nc"))))
    logging.info(
        f"Found {len(in_files_nrcan)} files, totalling {report_file_size(in_files_nrcan)}."
    )
    return dict(nrcan=in_files_nrcan)


def gather_wfdei_gem_capa(path: Union[str, os.PathLike]) -> Mapping[str, List[Path]]:
    # WFDEI-GEM-CaPa source data
    source_wfdei = Path(path)
    logging.info(f"Gathering WFDEI-GEM_CaPa from: {source_wfdei.as_posix()}")
    in_files_wfdei = list()
    for v in wfdei_gem_capa_variables:
        in_files_wfdei.extend(list(sorted(source_wfdei.rglob(f"{v}_*.nc"))))
    logging.info(
        f"Found {len(in_files_wfdei)} files, totalling {report_file_size(in_files_wfdei)}."
    )
    return {"wfdei-gem-capa": in_files_wfdei}


def gather_sc_earth(path: Union[str, os.PathLike]) -> Mapping[str, List[Path]]:
    # SC-Earth source data
    source_sc_earth = Path(path)
    logging.info(f"Gathering SC-Earth from: {source_sc_earth.as_posix()}")
    in_files_sc_earth = list()
    for v in sc_earth_variables:
        in_files_sc_earth.extend(
            list(sorted(source_sc_earth.rglob(f"SC-Earth_{v}_*.nc")))
        )
    logging.info(
        f"Found {len(in_files_sc_earth)} files, totalling {report_file_size(in_files_sc_earth)}."
    )
    return {"wfdei-gem-capa": in_files_sc_earth}
