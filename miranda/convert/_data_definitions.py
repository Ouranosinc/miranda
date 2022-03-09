import logging.config
import os
from pathlib import Path
from typing import Dict, List, Union

from miranda.scripting import LOGGING_CONFIG
from miranda.storage import report_file_size

from ._data import (
    era5_variables,
    nasa_ag_variables,
    nrcan_variables,
    sc_earth_variables,
    wfdei_gem_capa_variables,
)

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "gather_agcfsr",
    "gather_agmerra",
    "gather_era5_land",
    "gather_era5_single_levels",
    "gather_nrcan_gridded_obs",
    "gather_sc_earth",
    "gather_wfdei_gem_capa",
]


def gather_era5_single_levels(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # ERA5 source data
    source_era5 = Path(path)
    logging.info("Gathering ERA5 from %s" % source_era5.as_posix())
    infiles_era5 = list()
    for v in era5_variables:
        infiles_era5.extend(list(sorted(source_era5.rglob(f"{v}_*.nc"))))
    logging.info(
        f"Found {len(infiles_era5)} files, totalling {report_file_size(infiles_era5)}."
    )
    return {"era5-single-levels": infiles_era5}


def gather_era5_land_sea_mask(path: Union[str, os.PathLike]) -> Dict:
    try:
        land_sea_mask = dict(lsm=next(Path(path).glob("sftlf*.nc")))
    except StopIteration:
        logging.error("No land_sea_mask found for ERA5.")
        raise FileNotFoundError()
    return land_sea_mask


def gather_era5_land(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # ERA5-Land source data
    source_era5l = Path(path)
    logging.info("Gathering ERA5-Land from %s" % source_era5l.as_posix())
    infiles_era5l = list()
    for v in era5_variables:
        infiles_era5l.extend(list(sorted(source_era5l.rglob(f"{v}_*.nc"))))
    logging.info(
        f"Found {len(infiles_era5l)} files, totalling {report_file_size(infiles_era5l)}."
    )
    return {"era5-land": infiles_era5l}


def gather_agmerra(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # agMERRA source data
    source_agmerra = Path(path)
    logging.info("Gathering agMERRA from %s" % source_agmerra.as_posix())
    infiles_agmerra = list()
    for v in nasa_ag_variables:
        infiles_agmerra.extend(list(sorted(source_agmerra.rglob(f"AgMERRA_*_{v}.nc4"))))
    logging.info(
        f"Found {len(infiles_agmerra)} files, totalling {report_file_size(infiles_agmerra)}."
    )
    return dict(cfsr=infiles_agmerra)


def gather_agcfsr(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # agCFSR source data
    source_agcfsr = Path(path)
    logging.info("Gathering CFSR from %s" % source_agcfsr.as_posix())
    infiles_agcfsr = list()
    for v in nasa_ag_variables:
        infiles_agcfsr.extend(list(sorted(source_agcfsr.rglob(f"AgCFSR_*_{v}.nc4"))))
    logging.info(
        f"Found {len(infiles_agcfsr)} files, totalling {report_file_size(infiles_agcfsr)}."
    )
    return dict(cfsr=infiles_agcfsr)


def gather_nrcan_gridded_obs(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # NRCan Gridded Obs source data
    source_nrcan = Path(path)
    logging.info("Gathering NRCAN Gridded Obs from %s" % source_nrcan.as_posix())
    infiles_nrcan = list()
    for v in nrcan_variables:
        infiles_nrcan.extend(list(sorted(source_nrcan.joinpath(v).glob(f"*{v}_*.nc"))))
    logging.info(
        f"Found {len(infiles_nrcan)} files, totalling {report_file_size(infiles_nrcan)}."
    )
    return dict(nrcan=infiles_nrcan)


def gather_wfdei_gem_capa(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # WFDEI-GEM-CaPa source data
    source_wfdei = Path(path)
    logging.info("Gathering WFDEI-GEM_CaPa from %s" % source_wfdei.as_posix())
    infiles_wfdei = list()
    for v in wfdei_gem_capa_variables:
        infiles_wfdei.extend(list(sorted(source_wfdei.rglob(f"{v}_*.nc"))))
    logging.info(
        f"Found {len(infiles_wfdei)} files, totalling {report_file_size(infiles_wfdei)}."
    )
    return {"wfdei-gem-capa": infiles_wfdei}


def gather_sc_earth(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # SC_Earth source data
    source_sc_earth = Path(path)
    logging.info("Gathering SC-Earth from %s" % source_sc_earth.as_posix())
    infiles_sc_earth = list()
    for v in sc_earth_variables:
        infiles_sc_earth.extend(
            list(sorted(source_sc_earth.rglob(f"SC-Earth_{v}_*.nc")))
        )
    logging.info(
        f"Found {len(infiles_sc_earth)} files, totalling {report_file_size(infiles_sc_earth)}."
    )
    return {"wfdei-gem-capa": infiles_sc_earth}
