import logging.config
import os
from pathlib import Path
from typing import Dict, List, Union

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "gather_agcfsr",
    "gather_agmerra",
    "gather_era5",
    "gather_era5_land",
    "gather_nrcan_gridded_obs",
    "gather_sc_earth",
    "gather_wfdei_gem_capa",
]

ERA5_VARIABLES = [
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


NASA_AG_VARIABLES = ["prate", "rhstmax", "srad", "tavg", "tmax", "tmin", "wndpsd"]

WFDEI_GEM_CAPA_VARIABLES = ["huss", "pr", "ps", "rlds", "rsds", "sfcWind", "tas"]

SC_EARTH_VARIABLES = ["prcp", "tdew", "tmean", "trange", "wind"]


def gather_era5(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # ERA5 source data
    source_era5 = Path(path)
    logging.info("Gathering ERA5 from %s" % source_era5.as_posix())
    infiles_era5 = list()
    for v in ERA5_VARIABLES:
        infiles_era5.extend(list(sorted(source_era5.rglob(f"{v}_*.nc"))))
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
    for v in ERA5_VARIABLES:
        infiles_era5l.extend(list(sorted(source_era5l.rglob(f"{v}_*.nc"))))
    return {"era5-land": infiles_era5l}


def gather_agmerra(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # agMERRA source data
    source_agmerra = Path(path)
    logging.info("Gathering agMERRA from %s" % source_agmerra.as_posix())
    infiles_agmerra = list()
    for v in NASA_AG_VARIABLES:
        infiles_agmerra.extend(list(sorted(source_agmerra.rglob(f"AgMERRA_*_{v}.nc4"))))
    return dict(cfsr=infiles_agmerra)


def gather_agcfsr(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # agCFSR source data
    source_agcfsr = Path(path)
    logging.info("Gathering CFSR from %s" % source_agcfsr.as_posix())
    infiles_agcfsr = list()
    for v in NASA_AG_VARIABLES:
        infiles_agcfsr.extend(list(sorted(source_agcfsr.rglob(f"AgCFSR_*_{v}.nc4"))))
    return dict(cfsr=infiles_agcfsr)


def gather_nrcan_gridded_obs(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # NRCan Gridded Obs source data
    source_nrcan = Path(path)
    logging.info("Gathering NRCAN Gridded Obs from %s" % source_nrcan.as_posix())
    infiles_nrcan = list(sorted(source_nrcan.joinpath("tasmax").glob("*tasmax_*.nc")))
    infiles_nrcan.extend(
        list(sorted(source_nrcan.joinpath("tasmin").glob("*tasmin_*.nc")))
    )
    infiles_nrcan.extend(list(sorted(source_nrcan.joinpath("pr").glob("*pr*.nc"))))
    return dict(nrcan=infiles_nrcan)


def gather_wfdei_gem_capa(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # WFDEI-GEM-CaPa source data
    source_wfdei = Path(path)
    logging.info("Gathering WFDEI-GEM_CaPa from %s" % source_wfdei.as_posix())
    infiles_wfdei = list()
    for v in WFDEI_GEM_CAPA_VARIABLES:
        infiles_wfdei.extend(list(sorted(source_wfdei.rglob(f"{v}_*.nc"))))
    return {"wfdei-gem-capa": infiles_wfdei}


def gather_sc_earth(path: Union[str, os.PathLike]) -> Dict[str, List[Path]]:
    # SC_Earth source data
    source_sc_earth = Path(path)
    logging.info("Gathering SC-Earth from %s" % source_sc_earth.as_posix())
    infiles_sc_earth = list()
    for v in SC_EARTH_VARIABLES:
        infiles_sc_earth.extend(
            list(sorted(source_sc_earth.rglob(f"SC-Earth_{v}_*.nc")))
        )
    return {"wfdei-gem-capa": infiles_sc_earth}
