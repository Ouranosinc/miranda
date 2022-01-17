import logging.config
import os
from pathlib import Path
from typing import Dict, Union

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "gather_cfsr",
    "gather_era5",
    "gather_era5_land",
    "gather_merra2",
    "gather_nrcan_gridded_obs",
    "gather_wfdei_gem_capa",
]


def gather_era5(path: Union[str, os.PathLike]) -> (list, Dict):
    # ERA5 source data
    source_era5 = Path(path)
    logging.info("Gathering ERA5 from %s" % source_era5.as_posix())
    infiles_era5 = list(sorted(source_era5.rglob("tas_*.nc")))
    infiles_era5.extend(list(sorted(source_era5.rglob("pr_*.nc"))))
    land_sea_mask = dict(lsm=next(source_era5.glob("sftlf*.nc")))
    return infiles_era5, land_sea_mask


def gather_era5_land(path: Union[str, os.PathLike]) -> (list, None):
    # ERA5-Land source data
    source_era5l = Path(path)
    logging.info("Gathering ERA5-Land from %s" % source_era5l.as_posix())
    infiles_era5l = list(sorted(source_era5l.rglob("*tas_*.nc")))
    infiles_era5l.extend(list(sorted(source_era5l.rglob("*pr_*.nc"))))
    return infiles_era5l, None


def gather_merra2(path: Union[str, os.PathLike]) -> (list, Dict):
    # MERRA2 source data
    source_merra2 = Path(path)
    logging.info("Gathering MERRA2 from %s" % source_merra2.as_posix())
    infiles_merra2 = list(sorted(source_merra2.rglob("tas_*.nc")))
    infiles_merra2.extend(list(sorted(source_merra2.rglob("prbc*.nc"))))
    land_sea_mask = dict(
        sftlf=next(source_merra2.joinpath("sftlf").glob("sftlf*.nc")),
    )
    return infiles_merra2, land_sea_mask


def gather_cfsr(path: Union[str, os.PathLike]) -> (list, Dict):
    # CFSR source data
    source_cfsr = Path(path)
    logging.info("Gathering CFSR from %s" % source_cfsr.as_posix())
    infiles_cfsr = list(sorted(source_cfsr.joinpath("tas").glob("tas_*.nc")))
    infiles_cfsr.extend(
        list(sorted(source_cfsr.joinpath("tasmax").glob("tasmax_*.nc")))
    )
    infiles_cfsr.extend(
        list(sorted(source_cfsr.joinpath("tasmin").glob("tasmin_*.nc")))
    )
    infiles_cfsr.extend(list(sorted(source_cfsr.rglob("pr*.nc"))))
    land_sea_mask = dict(
        sftlf=next(source_cfsr.joinpath("sftlf").glob("sftlf*.nc")),
    )
    return infiles_cfsr, land_sea_mask


def gather_nrcan_gridded_obs(path: Union[str, os.PathLike]) -> (list, None):
    # NRCan Gridded Obs source data
    source_nrcan = Path(path)
    logging.info("Gathering NRCAN Gridded Obs from %s" % source_nrcan.as_posix())
    infiles_nrcan = list(sorted(source_nrcan.joinpath("tasmax").glob("*tasmax_*.nc")))
    infiles_nrcan.extend(
        list(sorted(source_nrcan.joinpath("tasmin").glob("*tasmin_*.nc")))
    )
    infiles_nrcan.extend(list(sorted(source_nrcan.joinpath("pr").glob("*pr*.nc"))))


def gather_wfdei_gem_capa(path: Union[str, os.PathLike]) -> (list, None):
    # WFDEI-GEM-CaPa source data
    source_wfdei = Path(path)
    logging.info("Gathering WFDEI-GEM_CaPa from %s" % source_wfdei.as_posix())
    infiles_wfdei = list(sorted(source_wfdei.glob("tas_*.nc")))
    infiles_wfdei.extend(list(sorted(source_wfdei.glob("pr_*.nc"))))
