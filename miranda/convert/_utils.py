import json
import logging.config
import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import dask.config
import numpy as np
import regionmask
import xarray
import xarray as xr
from clisops.core import subset
from dask.distributed import Client
from xclim.indices import tas

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


dask.config.set(local_directory=f"{Path(__file__).parent}/dask_workers/")


PROJECT_INSTITUTES = {
    "era5": "ecmwf",
    "era5-land": "ecmwf",
    "merra2": "nasa",
    "cfsr": "ncar",
    "nrcan-gridded-10km": "nrcan",
    "wfdei-gem-capa": "usask",
}

HOURLY_ACCUMULATED_VARIABLES = dict()
HOURLY_ACCUMULATED_VARIABLES["era5-land"] = [
    "e",
    "evabs",
    "evaow",
    "evatc",
    "evavt",
    "pev",
    "ro",
    "sf",
    "slhf",
    "smlt",
    "sro",
    "sshf",
    "ssr",
    "ssrd",
    "ssro",
    "str",
    "strd",
    "tp",
]


# Needed pre-processing function
def _drop_those_time_bnds(dataset: xr.Dataset):
    if "time_bnds" in dataset.variables:
        return dataset.drop_vars(["time_bnds"])
    return dataset


def reanalysis_processing(
    data: Dict[str, List[Union[str, os.PathLike]]],
    output_folder: Union[str, os.PathLike],
    variables: Sequence[str],
    time_freq: str,
    domains: Optional[Sequence[str]] = None,
    start: Optional[str] = None,
    end: Optional[str] = None,
    land_sea_mask: Optional[Dict] = None,
    land_sea_percentage: int = 50,
) -> None:
    """

    Parameters
    ----------
    data: Dict[str, List[str]]
    output_folder: Union[str, os.PathLike]
    variables: Sequence[str]
    time_freq: {"hourly", "daily"}
    domains: Sequence[str]
      Allowed options: {"QC", "CAN", "AMNO", None}
    start: str, optional
    end: str, optional
    land_sea_mask: Dict, optional
    land_sea_percentage: int

    Returns
    -------
    None
    """
    out_files = Path(output_folder)
    if domains is None:
        domains = ["raw"]

    for domain in domains:
        output_folder = out_files.joinpath(domain)
        output_folder.mkdir(exist_ok=True)
        for project, infiles in data.items():
            if domain != "raw":
                logging.info(f"Processing {project} data for domain {domain}.")
            else:
                logging.info(f"Processing {project} data.")
            for var in variables:
                institute = PROJECT_INSTITUTES[project]
                file_name = f"{var}_{time_freq}_{institute}_{project}_{domain.lower()}"

                # Select only for variable of interest
                multi_files = [x for x in infiles if var in str(x)]
                with Client():
                    logging.info("Resampling variable `%s`." % var)

                    # Subsetting operations
                    subset_geo = True
                    if domain.upper() == "QC":
                        lons = np.array([-79.76, -57.10])
                        lats = np.array([44.99, 62.59])
                    elif domain.upper() == "CAN":
                        lons = np.array([-141.02, -52.60])
                        lats = np.array([41.68, 83.14])
                    elif domain.upper() == "AMNO":
                        lons = np.array([-180, -10])
                        lats = np.array([10, 90])
                    elif domain.upper() == "RAW":
                        subset_geo = False
                        if start and end:
                            subset_time = True
                        else:
                            subset_time = False
                    else:
                        raise NotImplementedError()

                    if subset_geo:
                        ds = subset.subset_bbox(
                            xr.open_mfdataset(
                                multi_files,
                                combine="by_coords",
                                concat_dim=["time"],
                                chunks={"time": "auto"},
                                engine="netcdf4",
                                preprocess=_drop_those_time_bnds,
                            ),
                            lon_bnds=lons
                            if project not in {"era5", "cfsr"}
                            else lons + 360,
                            lat_bnds=lats,
                            start_date=start,
                            end_date=end,
                        )
                    elif subset_time:
                        ds = subset.subset_time(
                            xr.open_mfdataset(
                                multi_files,
                                combine="by_coords",
                                concat_dim=["time"],
                                chunks={"time": "auto"},
                                engine="netcdf4",
                                preprocess=_drop_those_time_bnds,
                            ),
                            start_date=start,
                            end_date=end,
                        )
                    else:
                        ds = xr.open_mfdataset(
                            multi_files,
                            combine="by_coords",
                            concat_dim=["time"],
                            chunks={"time": "auto"},
                            engine="netcdf4",
                            preprocess=_drop_those_time_bnds,
                        )

                    ds = variable_conversion(ds, project=project)

                    # Land-Sea mask operations
                    if project in land_sea_mask.keys():
                        logging.info(
                            f"Masking variable with land-sea mask at {land_sea_percentage} % cutoff."
                        )
                        land_sea_mask, lsm_file = land_sea_mask[project]
                        lsm_raw = xr.open_dataset(lsm_file)
                        try:
                            lsm_raw = lsm_raw.rename(
                                {"longitude": "lon", "latitude": "lat"}
                            )
                        except ValueError:
                            raise

                        lsm = subset.subset_bbox(
                            lsm_raw,
                            lon_bnds=lons
                            if project not in ["era5", "cfsr"]
                            else lons + 360,
                            lat_bnds=lats,
                        ).load()
                        lsm = lsm.where(
                            lsm[land_sea_mask] > float(land_sea_percentage) / 100
                        )
                        if project == "era5":
                            ds = ds.where(
                                lsm[land_sea_mask].isel(time=0, drop=True).notnull()
                            )
                            try:
                                ds = ds.rename({"longitude": "lon", "latitude": "lat"})
                            except ValueError:
                                raise
                        elif project in ["merra2", "cfsr"]:
                            ds = ds.where(lsm[land_sea_mask].notnull())
                        ds.attrs["land_sea_cutoff"] = f"{land_sea_percentage} %"

                    if time_freq.lower() == "daily":
                        daily_aggregation(ds, project, var, file_name, output_folder)
                    else:
                        logging.info("Writing out fixed files for %s." % file_name)
                        years, datasets = zip(*ds.groupby("time.year"))  # noqa
                        out_filenames = [
                            output_folder.joinpath(f"{file_name}_{year}.nc")
                            for year in years
                        ]
                        xr.save_mfdataset(datasets, out_filenames)


def variable_conversion(ds: xarray.Dataset, project: str) -> xarray.Dataset:
    """Convert variables to CF-compliant format"""

    # Add and update existing metadata fields
    def _metadata_conversion(d: xarray.Dataset, p: str) -> xarray.Dataset:
        if p in ["era5", "era5-single-levels", "era5-land"]:
            metadata_definition = json.load(
                open(Path(__file__).parent / "ecmwf_cf_attrs.json")
            )
        else:
            raise NotImplementedError()

        # Add global attributes
        d.attrs.update(metadata_definition["Header"])
        descriptions = metadata_definition["variable_entry"]

        # Add variable metadata
        for v in d.data_vars:
            d[v].attrs.update(descriptions[v])

        # Rename data variables
        for v in d.data_vars:
            try:
                d = d.rename({descriptions[v]["original_variable"]: v})
            except (ValueError, IndexError):
                pass
        return d

    # For converting variable units
    def _units_conversion(d: xarray.Dataset, p: str) -> xarray.Dataset:

        if p in HOURLY_ACCUMULATED_VARIABLES.keys():
            if any([v in d.data_vars for v in HOURLY_ACCUMULATED_VARIABLES[p]]):
                d["time"] = d.time - np.timedelta64(1, "h")

        for v in d.data_vars:
            if hasattr(d[v].attrs, "_conversion"):
                d[v] = d[v] * d[v].attrs["_conversion"]
            del d[v].attrs["_conversion"]
        return d

    ds = _metadata_conversion(ds, project)
    ds = _units_conversion(ds, project)

    return ds


def daily_aggregation(
    ds, project: str, variable: str, input_file: Union[str, os.PathLike], output_folder
) -> None:
    # Daily variable aggregation operations
    input_file = Path(input_file)
    output_folder = Path(output_folder)

    logging.info("Creating daily upscaled reanalyses.")
    if variable == "tas":
        # Some looping to deal with memory consumption issues
        for v, func in {
            "tasmax": "max",
            "tasmin": "min",
            "tas": "mean",
        }.items():
            if project in ["cfsr", "nrcan"]:
                v_desired = v
            else:
                v_desired = "tas"

            input_file_parts = input_file.stem.split("_")
            input_file_parts[1] = "daily"
            output = output_folder.joinpath(f"{'_'.join(input_file_parts)}")

            if any([f for f in output.glob("*")]):
                logging.info("Files for `%s` exist. Continuing..." % output.name)
                return

            ds_out = xr.Dataset()
            if v == "tas" and not hasattr(ds, "tas"):
                ds_out[v] = tas(tasmax=ds.tasmax, tasmin=ds.tasmin)
            else:
                # Thanks for the help, xclim contributors
                r = ds[v_desired].resample(time="D", keep_attrs=True)
                ds_out[v] = getattr(r, func)(dim="time", keep_attrs=True)

            logging.info("Masking %s with AR6 regions." % v)
            mask = regionmask.defined_regions.ar6.all.mask(ds_out.lon, ds_out.lat)
            ds_out = ds_out.assign_coords(region=mask)

            logging.info("Writing out daily converted %s." % output)
            years, datasets = zip(*ds_out.groupby("time.year"))  # noqa
            out_filenames = [
                output_folder.joinpath(f"{output}_{year}.nc") for year in years
            ]
            xr.save_mfdataset(datasets, out_filenames)

            del ds_out

    if variable in ["pr"]:
        input_file_parts = input_file.stem.split("_")
        input_file_parts[1] = "daily"
        output = output_folder.joinpath(f"{'_'.join(input_file_parts)}")

        if any([f for f in output.glob("*")]):
            logging.info("Files for `%s` exist. Continuing..." % output.name)
            return

        ds_out = xr.Dataset()
        logging.info("Converting precipitation units")
        if project in HOURLY_ACCUMULATED_VARIABLES.keys():
            ds_out["pr"] = ds.pr.resample(time="D").max(dim="time", keep_attrs=True)
        else:
            ds_out["pr"] = ds.pr.resample(time="D").mean(dim="time", keep_attrs=True)

        logging.info("Masking %s with AR6 regions." % variable)
        mask = regionmask.defined_regions.ar6.all.mask(ds_out.lon, ds_out.lat)
        ds_out = ds_out.assign_coords(region=mask)

        logging.info("Writing out daily converted %s." % output)
        years, datasets = zip(*ds_out.groupby("time.year"))  # noqa
        out_filenames = [
            output_folder.joinpath(f"{output}_{year}.nc") for year in years
        ]
        xr.save_mfdataset(datasets, out_filenames)

        del ds_out

    return
