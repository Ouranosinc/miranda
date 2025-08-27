from __future__ import annotations
import datetime
import json
import logging
import os
import re
from pathlib import Path

from miranda.storage import report_file_size


logger = logging.getLogger("miranda.convert.data_definitions")

__all__ = [
    "eccc_rdrs_variables",
    "era5_variables",
    "gather_agcfsr",
    "gather_agmerra",
    "gather_eccc_rdrs",
    "gather_ecmwf",
    "gather_emdna",
    "gather_grnch",
    "gather_nex",
    "gather_nrcan_gridded_obs",
    "gather_raw_rdrs_by_years",
    "gather_sc_earth",
    "gather_wfdei_gem_capa",
    "nasa_ag_variables",
    "nrcan_variables",
    "sc_earth_variables",
    "wfdei_gem_capa_variables",
]

_data_folder = Path(__file__).resolve().parent / "data"


eccc_rdrs_variables = {}
eccc_rdrs_variables["raw"] = [v for v in json.load(_data_folder.joinpath("eccc_rdrs_cf_attrs.json").open("r", encoding="utf-8"))["variables"].keys()]
eccc_rdrs_variables["cf"] = [
    attrs["_cf_variable_name"]
    for attrs in json.load(_data_folder.joinpath("eccc_rdrs_cf_attrs.json").open("r", encoding="utf-8"))["variables"].values()
    if "_cf_variable_name" in attrs
]

era5_variables = json.load(_data_folder.joinpath("ecmwf_cf_attrs.json").open("r", encoding="utf-8"))["variables"].keys()
grnch_variables = ["T", "Tmin", "Tmax", "P"]
nrcan_variables = ["tasmin", "tasmax", "pr"]
nasa_ag_variables = json.load(_data_folder.joinpath("nasa_cf_attrs.json").open("r", encoding="utf-8"))["variables"].keys()
sc_earth_variables = ["prcp", "tdew", "tmean", "trange", "wind"]
wfdei_gem_capa_variables = json.load(_data_folder.joinpath("usask_cf_attrs.json").open())["variables"].keys()


def _gather(
    name: str,
    variables: list[str],
    source: str | os.PathLike,
    glob_pattern: str,
    suffix: str | None = None,
    recursive: bool | None = False,
) -> dict[str, list[Path]]:
    source = Path(source).expanduser()
    msg = f"Gathering {name} files from: {source.as_posix()}"
    logger.info(msg)
    in_files = []
    for variable in variables:
        if suffix:
            pattern = glob_pattern.format(variable=variable, name=name, suffix=suffix)
        else:
            pattern = glob_pattern.format(variable=variable)
        if recursive:
            in_files.extend(list(sorted(source.rglob(pattern))))
        else:
            in_files.extend(list(sorted(source.glob(pattern))))
    msg = f"Found {len(in_files)} files, totalling {report_file_size(in_files)}."

    logger.info(msg)
    return {name: in_files}


def gather_ecmwf(
    project: str,
    path: str | os.PathLike,
    back_extension: bool = False,
    monthly_means: bool = False,
) -> dict[str, list[Path]]:
    """
    Gather ECMWF source data.

    Parameters
    ----------
    project : {"era5-single-levels", "era5-pressure-levels", "era5-land"}
    path : str or os.PathLike
    back_extension : bool
    monthly_means : bool

    Returns
    -------
    dict[str, list[pathlib.Path]]
    """
    name = f"{project}{'-monthly-means' if monthly_means else ''}{'-preliminary-back-extension' if back_extension else ''}"
    glob_pattern = "".join(["{variable}", f"_*_{name}_*.nc"])

    return _gather(name, era5_variables, source=path, glob_pattern=glob_pattern)


def gather_agmerra(path: str | os.PathLike) -> dict[str, list[Path]]:
    """
    Gather agMERRA source data.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict[str, list[pathlib.Path]]
    """
    return _gather("merra", nasa_ag_variables, source=path, glob_pattern="AgMERRA_*_{variable}.nc4")


def gather_agcfsr(path: str | os.PathLike) -> dict[str, list[Path]]:
    """
    Gather agCFSR source data.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict[str, list[pathlib.Path]]
    """
    return _gather("cfsr", nasa_ag_variables, source=path, glob_pattern="AgCFSR_*_{variable}.nc4")


def gather_nrcan_gridded_obs(path: str | os.PathLike) -> dict[str, list[Path]]:
    """
    Gather NRCan Gridded Observations source data.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict(str, list[pathlib.Path])
    """
    return _gather("nrcan", nrcan_variables, source=path, glob_pattern="*{variable}_*.nc")


def gather_wfdei_gem_capa(path: str | os.PathLike) -> dict[str, list[Path]]:
    """
    Gather WFDEI-GEM-CaPa source data.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict[str, list[pathlib.Path]]
    """
    return _gather(
        "wfdei-gem-capa",
        wfdei_gem_capa_variables,
        source=path,
        glob_pattern="{variable}_*.nc",
    )


def gather_sc_earth(path: str | os.PathLike) -> dict[str, list[Path]]:
    """
    Gather SC-Earth source data

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict[str, list[pathlib.Path]]
    """
    return _gather(
        "sc-earth",
        sc_earth_variables,
        source=path,
        glob_pattern="SC-Earth_{variable}_*.nc",
    )


def gather_eccc_rdrs(name: str, path: str | os.PathLike, suffix: str, key: str) -> dict[str, dict[str, list[Path]]]:
    """
    Gather RDRS processed source data.

    Parameters
    ----------
    name : str
        The variable to gather.
    path : str or os.PathLike
        The location of the source data.
    suffix : str
        The filename suffix.
    key : {"raw", "cf"}
        Indicating which variable name dictionary to search for.

    Returns
    -------
    dict[str, list[pathlib.Path]]
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
            recursive=False,
        )
        files[name][vv] = tmp[name]
    return files


def gather_raw_rdrs_by_years(
    path: str | os.PathLike,
    project: str,
) -> dict[str, dict[str, list[Path]]]:
    """
    Gather raw RDRS files for preprocessing.

    Parameters
    ----------
    path: str or os.PathLike
    project: str

    Returns
    -------
    dict[str, dict[str, list[pathlib.Path]]
    """
    # Time stamps starts at noon and flow into subsequent months
    # Need full year plus previous december in order to easily produce complete hourly frequency monthly files
    path = Path(path)
    year_sets = dict()
    for year in range(1950, datetime.datetime.now().year + 1):
        dec_prev_year_files = []
        this_year_files = []

        for file in path.glob("*.nc"):
            match = re.search(r"(\d{10})", file.name)  # search for 10 digits (YYYYMMDDHH)
            if match:
                date_str = match.group(1)
                dt = datetime.datetime.strptime(date_str, "%Y%m%d%H")
                if dt.year == year - 1 and dt.month == 12:
                    dec_prev_year_files.append(file)
                elif dt.year == year:
                    this_year_files.append(file)

        # if there are files from the previous December, get the last one
        dec_prev_year_files.sort()
        if dec_prev_year_files:
            files = [dec_prev_year_files[-1]]
        else:
            files = []

        this_year_files.sort()
        files.extend(this_year_files)
        year_sets[str(year)] = files

    return {project: year_sets}


def gather_grnch(path: str | os.PathLike) -> dict[str, list[Path]]:
    """
    Gather raw ETS-GRNCH files for preprocessing.

    Parameters
    ----------
    path: str or os.PathLike

    Returns
    -------
    dict(str, dict(str, list[Path])) or None
    """
    # GRNCH-ETS source data
    source_grnch = Path(path)
    msg = f"Gathering GRNCH from: {source_grnch.as_posix()}"
    logger.info(msg)
    in_files_grnch = list()
    for v in grnch_variables:
        for yyyy in range(1970, 2020):
            in_files_grnch.extend(list(source_grnch.rglob(f"{v}_{yyyy}.nc")))
    msg = f"Found {len(in_files_grnch)} files, totalling {report_file_size(in_files_grnch)}."

    logger.info(msg)
    return dict(cfsr=sorted(in_files_grnch))


def gather_nex(
    path: str | os.PathLike,
) -> dict[str, list[Path]]:
    """
    Gather raw NEX files for preprocessing.

    Put all files that should be contained in one dataset in one entry of the dictionary.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict[str, list[pathlib.Path]]
    """
    source = Path(path)
    datasets = source.glob("*/*/*/*/*/*/*/*/*/")

    out_dict = dict()
    # separate files by datasets
    for dataset in datasets:
        in_files = list()
        in_files.extend(list(sorted(dataset.glob("*.nc"))))
        out_dict[str(dataset)] = in_files
    return out_dict


def gather_emdna(
    path: str | os.PathLike,
) -> dict[str, list[Path]]:
    """
    Gather raw EMDNA files for preprocessing.

    Put all files with the same member together.

    Parameters
    ----------
    path : str or os.PathLike

    Returns
    -------
    dict[str, list[pathlib.Path]]
    """
    source = Path(path)
    member_dict = {}
    # 100 members
    members = [f"{i:03d}" for i in range(1, 101)]
    for member in members:
        member_dict[member] = list(sorted(source.glob(f"EMDNA_estimate/*/EMDNA_*.{member}.nc4")))

    # OI
    member_dict["OI"] = list(sorted(source.glob("OI_estimate/*.nc4")))

    return member_dict
