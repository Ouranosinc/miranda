import datetime
import json
import logging.config
import os
from pathlib import Path
from typing import Dict, List, Optional, Union

from miranda.scripting import LOGGING_CONFIG
from miranda.storage import report_file_size

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "era5_variables",
    "eccc_rdrs_variables",
    "gather_agcfsr",
    "gather_agmerra",
    "gather_ecmwf",
    "gather_grnch",
    "gather_nex",
    "gather_nrcan_gridded_obs",
    "gather_raw_rdrs_by_years",
    "gather_rdrs",
    "gather_sc_earth",
    "gather_wfdei_gem_capa",
    "nasa_ag_variables",
    "nrcan_variables",
    "project_institutes",
    "sc_earth_variables",
    "wfdei_gem_capa_variables",
    "xarray_frequencies_to_cmip6like",
]

_data_folder = Path(__file__).parent / "data"

eccc_rdrs_variables = dict()
eccc_rdrs_variables["raw"] = [
    v
    for v in json.load(open(_data_folder / "eccc_rdrs_cf_attrs.json"))[
        "variables"
    ].keys()
]
eccc_rdrs_variables["cf"] = [
    attrs["_cf_variable_name"]
    for attrs in json.load(open(_data_folder / "eccc_rdrs_cf_attrs.json"))[
        "variables"
    ].values()
]

era5_variables = json.load(open(_data_folder / "ecmwf_cf_attrs.json"))[
    "variables"
].keys()
grnch_variables = ["T", "Tmin", "Tmax", "P"]
nrcan_variables = ["tasmin", "tasmax", "pr"]
nasa_ag_variables = json.load(open(_data_folder / "nasa_cf_attrs.json"))[
    "variables"
].keys()
sc_earth_variables = ["prcp", "tdew", "tmean", "trange", "wind"]
wfdei_gem_capa_variables = json.load(open(_data_folder / "usask_cf_attrs.json"))[
    "variables"
].keys()

project_institutes = {
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
    "rdrs-v21": "eccc",
    "NEX-GDDP-CMIP6": "nasa",
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
    variables: List[str],
    source: Union[str, os.PathLike],
    glob_pattern: str,
    suffix: Optional[str] = None,
    recursive: Optional[bool] = False,
) -> Dict[str, List[Path]]:
    source = Path(source).expanduser()
    logging.info(f"Gathering {name} files from: {source.as_posix()}")
    in_files = list()
    for variable in variables:
        if suffix:
            pattern = glob_pattern.format(variable=variable, name=name, suffix=suffix)
        else:
            pattern = glob_pattern.format(variable=variable)
        if recursive:
            in_files.extend(list(sorted(source.rglob(pattern))))
        else:
            in_files.extend(list(sorted(source.glob(pattern))))
    logging.info(
        f"Found {len(in_files)} files, totalling {report_file_size(in_files)}."
    )
    return {name: in_files}


def gather_ecmwf(
    project: str,
    path: Union[str, os.PathLike],
    back_extension: bool = False,
    monthly_means: bool = False,
) -> Dict[str, List[Path]]:
    """

    Parameters
    ----------
    project : {"era5-single-levels", "era5-pressure-levels", "era5-land"}
    path : str or os.PathLike
    back_extension : bool
    monthly_means : bool

    Returns
    -------
    dict(str, list[pathlib.Path])
    """
    # ERA5-Single-Levels source data
    name = (
        f"{project}"
        f"{'-monthly-means' if monthly_means else ''}"
        f"{'-preliminary-back-extension' if back_extension else ''}"
    )
    glob_pattern = "".join(["{variable}", f"_*_{name}_*.nc"])

    return _gather(name, era5_variables, source=path, glob_pattern=glob_pattern)


def gather_agmerra(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    """Gather agMERRA source data.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict(str, list[pathlib.Path])
    """
    return _gather(
        "merra", nasa_ag_variables, source=path, glob_pattern="AgMERRA_*_{variable}.nc4"
    )


def gather_agcfsr(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    """Gather agCFSR source data.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict(str, list[pathlib.Path])
    """
    return _gather(
        "cfsr", nasa_ag_variables, source=path, glob_pattern="AgCFSR_*_{variable}.nc4"
    )


def gather_nrcan_gridded_obs(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    """Gather NRCan Gridded Observations source data.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict(str, list[pathlib.Path])
    """
    return _gather(
        "nrcan", nrcan_variables, source=path, glob_pattern="*{variable}_*.nc"
    )


def gather_wfdei_gem_capa(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    """Gather WFDEI-GEM-CaPa source data.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict(str, list[pathlib.Path])
    """
    return _gather(
        "wfdei-gem-capa",
        wfdei_gem_capa_variables,
        source=path,
        glob_pattern="{variable}_*.nc",
    )


def gather_sc_earth(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    """Gather SC-Earth source data

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict(str, list[pathlib.Path])
    """
    return _gather(
        "sc-earth",
        sc_earth_variables,
        source=path,
        glob_pattern="SC-Earth_{variable}_*.nc",
    )


def gather_rdrs(
    name: str, path: Union[str, os.PathLike], suffix: str, key: str
) -> Dict[str, Dict[str, List[Path]]]:
    """Gather RDRS processed source data.

    Parameters
    ----------
    name : str
    path : str or os.PathLike
    suffix : str
    key : str  one of 'raw' or 'cf' indicating which variable name dictionary to search for

    Returns
    -------
    dict(str, list[pathlib.Path])
    """
    if isinstance(path, str):
        path = Path(path).expanduser()

    files = dict({name: dict()})
    for vv in eccc_rdrs_variables[key]:
        tmp = _gather(
            name,
            [vv],
            source=path.joinpath(vv),
            glob_pattern="{variable}_*_{name}_*.{suffix}",
            suffix=suffix,
            recursive=True,
        )
        files[name][vv] = tmp[name]
    return files


def gather_raw_rdrs_by_years(
    path: Union[str, os.PathLike]
) -> Dict[str, Dict[str, List[Path]]]:
    """Gather raw RDRS files for preprocessing.

    Parameters
    ----------
    path: str or os.PathLike

    Returns
    -------
    dict(str, dict(str, list[Path])) or None
    """
    # Time stamps starts at noon and flow into subsequent months
    # Need full year plus previous december in order to easily produce complete hourly frequency monthly files
    path = Path(path)
    year_sets = dict()
    for year in range(1950, datetime.datetime.now().year + 1):
        files = sorted(list(path.glob(f"{year - 1}12*.nc")))
        if files:
            files = [files[-1]]
        files.extend(sorted(list(path.glob(f"{year}*.nc"))))
        year_sets[str(year)] = files
    return {"rdrs-v21": year_sets}


def gather_grnch(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # GRNCH-ETS source data
    source_grnch = Path(path)
    logging.info(f"Gathering GRNCH from: {source_grnch.as_posix()}")
    in_files_grnch = list()
    for v in grnch_variables:
        for yyyy in range(1970, 2020):
            in_files_grnch.extend(list(source_grnch.rglob(f"{v}_{yyyy}.nc")))
    logging.info(
        f"Found {len(in_files_grnch)} files, totalling {report_file_size(in_files_grnch)}."
    )
    return dict(cfsr=sorted(in_files_grnch))


def gather_nex(
    path: Union[str, os.PathLike],
) -> Dict[str, List[Path]]:
    """Put all files that should be contained in one dataset in one entry of the dictionnary.

    Parameters
    ----------
    path : str or os.PathLike
    back_extension : bool
    monthly_means : bool

    Returns
    -------
    dict(str, list[pathlib.Path])
    """

    source = Path(path)
    datasets = source.glob("*/*/*/*/*/*/*/*/*/")

    out_dict = {}
    # separate files by datasets
    for dataset in datasets:
        in_files = list()
        in_files.extend(list(sorted(dataset.glob("*.nc"))))
        out_dict[str(dataset)] = in_files
    return out_dict
