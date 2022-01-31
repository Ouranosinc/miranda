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
from xclim.indices import tas

from miranda.gis.subset import subsetting_domains
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


dask.config.set(local_directory=f"{Path(__file__).parent}/dask_workers/")


PROJECT_INSTITUTES = {
    "cfsr": "ncar",
    "era5": "ecmwf",
    "era5-land": "ecmwf",
    "merra2": "nasa",
    "nrcan-gridded-10km": "nrcan",
    "wfdei-gem-capa": "usask",
}
PROJECT_LONS_FROM_0TO360 = ["era5", "cfsr"]

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
    aggregate: Union[str, bool] = False,
    domains: Optional[Sequence[str]] = None,
    start: Optional[str] = None,
    end: Optional[str] = None,
) -> None:
    """

    Parameters
    ----------
    data: Dict[str, List[str]]
    output_folder: Union[str, os.PathLike]
    variables: Sequence[str]
    aggregate: {"daily", False}
    domains: Sequence[str]
      Allowed options: {"QC", "CAN", "AMNO", None}
    start: str, optional
    end: str, optional

    Returns
    -------
    None
    """
    out_files = Path(output_folder)
    if domains is None:
        domains = [None]

    for domain in domains:
        if domain is not None:
            output_folder = out_files.joinpath(domain)  # noqa
        else:
            output_folder = output_folder
        output_folder.mkdir(exist_ok=True)

        for project, in_files in data.items():
            if domain is not None:
                logging.info(f"Processing {project} data for domain {domain}.")
            else:
                logging.info(f"Processing {project} data.")
            for var in variables:
                # Select only for variable of interest
                multi_files = sorted(x for x in in_files if var in str(x))
                logging.info("Resampling variable `%s`." % var)

                if aggregate:
                    time_freq = aggregate
                else:
                    time_freq = multi_files[0].split("_")[1]

                institute = PROJECT_INSTITUTES[project]
                file_name = "_".join([var, time_freq, institute, project])
                if domain is not None:
                    file_name = f"{file_name}_{domain}"

                # Subsetting operations
                ds = None
                subset_time = False
                if domain is None:
                    if start and end:
                        subset_time = True

                elif domain is not None:
                    if domain.upper() == "AMNO":
                        domain = "NAM"
                    region = subsetting_domains(domain)
                    lon_values = region[1], region[3]
                    lat_values = region[0], region[2]

                    ds = subset.subset_bbox(
                        xr.open_mfdataset(
                            multi_files,
                            chunks={"time": "auto"},
                            engine="netcdf4",
                            preprocess=_drop_those_time_bnds,
                        ),
                        lon_bnds=lon_values
                        if project not in {"era5", "cfsr"}
                        else lon_values + 360,
                        lat_bnds=lat_values,
                        start_date=start,
                        end_date=end,
                    )
                if subset_time:
                    ds = subset.subset_time(
                        xr.open_mfdataset(
                            multi_files,
                            chunks={"time": "auto"},
                            engine="netcdf4",
                            preprocess=_drop_those_time_bnds,
                        ),
                        start_date=start,
                        end_date=end,
                    )
                elif not any([subset_time, domain, ds]):
                    ds = xr.open_mfdataset(
                        multi_files,
                        chunks={"time": "auto"},
                        engine="netcdf4",
                        preprocess=_drop_those_time_bnds,
                    )

                if project in HOURLY_ACCUMULATED_VARIABLES.keys():
                    if all([v in HOURLY_ACCUMULATED_VARIABLES[project] for v in ds.data_vars]):
                        ds = deaccumulate(ds, project)
                        #ds["time"] = ds.time - np.timedelta64(1, "h")

                ds = variable_conversion(ds, project=project)

                if time_freq.lower() == "daily":
                    daily_aggregation(ds, project, file_name, output_folder)
                else:
                    logging.info("Writing out fixed files for %s." % file_name)
                    # TODO: Resample for year*months
                    years, datasets = zip(*ds.groupby("time.year"))  # noqa
                    out_filenames = [
                        output_folder.joinpath(f"{file_name}_{year}.nc")
                        for year in years
                    ]
                    xr.save_mfdataset(datasets, out_filenames)

def deaccumulate(ds: xarray.Dataset, project: str):

    ds


def variable_conversion(ds: xarray.Dataset, project: str) -> xarray.Dataset:
    """Convert variables to CF-compliant format"""

    # Add and update existing metadata fields
    def _metadata_conversion(d: xarray.Dataset, p: str) -> xarray.Dataset:
        if p in ["era5", "era5-single-levels", "era5-land"]:
            metadata_definition = json.load(
                open(Path(__file__).parent.parent / "ecmwf" / "ecmwf_cf_attrs.json")
            )
        else:
            raise NotImplementedError()

        # Add global attributes
        d.attrs.update(metadata_definition["Header"])
        d.attrs.update(dict(project=p))
        descriptions = metadata_definition["variable_entry"]

        # Add variable metadata
        for v in d.data_vars:
            d[v].attrs.update(descriptions[v])

        # Rename data variables
        for v in d.data_vars:
            try:
                cf_name = descriptions[v]["_cf_variable_name"]
                d = d.rename({v: cf_name})
                d[cf_name].attrs.update(dict(original_variable=v))
                del d[cf_name].attrs["_cf_variable_name"]
            except (ValueError, IndexError):
                pass
        return d

    # For renaming lat and lon dims
    def _dims_conversion(d: xarray.Dataset):
        for orig, new in dict(longitude="lon", latitude="lat").items():
            try:
                d = d.rename({orig: new})
            except KeyError:
                pass
        return d

    # For converting variable units
    def _units_conversion(d: xarray.Dataset, p: str) -> xarray.Dataset:

        # if p in HOURLY_ACCUMULATED_VARIABLES.keys():
        #     if any([v in d.data_vars for v in HOURLY_ACCUMULATED_VARIABLES[p]]):
        #         d["time"] = d.time - np.timedelta64(1, "h")

        # TODO: We need to de-accumulate this variable as well

        for v in d.data_vars:
            if hasattr(d[v].attrs, "_conversion"):
                d[v] = d[v] * d[v].attrs["_conversion"]
                del d[v].attrs["_conversion"]
        return d

    ds = _metadata_conversion(ds, project)
    ds = _dims_conversion(ds)
    ds = _units_conversion(ds, project)

    return ds


def daily_aggregation(
    ds, project: str, input_file: Union[str, os.PathLike], output_folder
) -> None:
    # Daily variable aggregation operations
    input_file = Path(input_file)
    output_folder = Path(output_folder)
    logging.info("Creating daily upscaled reanalyses.")

    for variable in ds.data_vars:
        input_file_parts = input_file.stem.split("_")
        input_file_parts[0] = variable
        input_file_parts[1] = "daily"
        output = output_folder.joinpath(f"{'_'.join(input_file_parts)}")

        if any([f for f in output.glob("*")]):
            logging.info("Files for `%s` exist. Continuing..." % output.name)
            return

        if variable == "tas":
            # Some looping to deal with memory consumption issues
            # TODO: Add cell methods for tasmax and tasmin
            for v, func in {
                "tasmax": "max",
                "tasmin": "min",
                "tas": "mean",
            }.items():
                if project in ["cfsr", "nrcan"]:
                    v_desired = v
                else:
                    v_desired = "tas"

                input_file_parts[0] = v
                output = output_folder.joinpath(f"{'_'.join(input_file_parts)}")

                ds_out = xr.Dataset()
                ds_out.attrs = ds.attrs.copy()
                if v == "tas" and not hasattr(ds, "tas"):
                    ds_out[v] = tas(tasmax=ds.tasmax, tasmin=ds.tasmin)
                else:
                    # Thanks for the help, xclim contributors
                    r = ds[v_desired].resample(time="D", keep_attrs=True)
                    ds_out[v] = getattr(r, func)(dim="time", keep_attrs=True)

                logging.info("Writing out daily converted %s." % output)
                years, datasets = zip(*ds_out.groupby("time.year"))  # noqa
                out_filenames = [
                    output_folder.joinpath(f"{output}_{year}.nc") for year in years
                ]
                xr.save_mfdataset(datasets, out_filenames)

                del ds_out
            return

        if variable in ["pr"]:
            ds_out = xr.Dataset()
            ds_out.attrs = ds.attrs.copy()
            logging.info("Converting precipitation units")
            if project in HOURLY_ACCUMULATED_VARIABLES.keys():
                ds_out["pr"] = ds.pr.resample(time="D").max(dim="time", keep_attrs=True)
            else:
                ds_out["pr"] = ds.pr.resample(time="D").mean(
                    dim="time", keep_attrs=True
                )
        else:
            continue

        logging.info("Writing out daily converted %s." % output)
        years, datasets = zip(*ds_out.groupby("time.year"))  # noqa
        out_filenames = [
            output_folder.joinpath(f"{output}_{year}.nc") for year in years
        ]
        xr.save_mfdataset(datasets, out_filenames)

        del ds_out
    return


def add_ar6_regions(ds: xarray.Dataset) -> xarray.Dataset:
    """Add the IPCC AR6 Regions to dataset.

    Parameters
    ----------
    ds : xarray.Dataset

    Returns
    -------
    xarray.Dataset
    """
    mask = regionmask.defined_regions.ar6.all.mask(ds.lon, ds.lat)
    ds = ds.assign_coords(region=mask)
    return ds


def threshold_land_sea_mask(
    ds: Union[xr.Dataset, str, os.PathLike],
    *,
    land_sea_mask: Dict[str, Union[os.PathLike, str]],
    land_sea_percentage: int = 50,
    output_folder: Optional[Union[str, os.PathLike]] = None,
) -> Optional[Path]:
    """Land-Sea mask operations.

    Parameters
    ----------
    ds: Union[xr.Dataset, str, os.PathLike]
    land_sea_mask: dict
    land_sea_percentage: int
    output_folder: str or os.PathLike, optional

    Returns
    -------
    Path
    """
    file_name = ""
    if isinstance(ds, (str, os.PathLike)):
        if output_folder is not None:
            output_folder = Path(output_folder)
            file_name = f"{Path(ds).stem}_land-sea-masked.nc"
        ds = xr.open_dataset(ds)

    if output_folder is not None and file_name == "":
        logging.warning(
            "Cannot generate filenames from xarray.Dataset objects. Consider writing NetCDF manually."
        )

    try:
        project = ds.attrs["project"]
    except KeyError:
        raise ValueError("No 'project' found for given dataset.")

    if project in PROJECT_LONS_FROM_0TO360:
        add_lon_values = 360
    else:
        add_lon_values = 0

    if project in land_sea_mask.keys():
        logging.info(
            f"Masking variable with land-sea mask at {land_sea_percentage} % cutoff."
        )
        land_sea_mask_variable, lsm_file = land_sea_mask[project]
        lsm_raw = xr.open_dataset(lsm_file)
        try:
            lsm_raw = lsm_raw.rename({"longitude": "lon", "latitude": "lat"})
        except ValueError:
            raise

        lon_bounds = np.array([ds.lon.min(), ds.lon.max()]) + add_lon_values
        lat_bounds = np.array([ds.lat.min(), ds.lat.max()])

        lsm = subset.subset_bbox(
            lsm_raw,
            lon_bnds=lon_bounds,
            lat_bnds=lat_bounds,
        ).load()
        lsm = lsm.where(lsm[land_sea_mask_variable] > float(land_sea_percentage) / 100)
        if project == "era5":
            ds = ds.where(lsm[land_sea_mask].isel(time=0, drop=True).notnull())
            try:
                ds = ds.rename({"longitude": "lon", "latitude": "lat"})
            except ValueError:
                raise
        elif project in ["merra2", "cfsr"]:
            ds = ds.where(lsm[land_sea_mask].notnull())

        ds.attrs["land_sea_cutoff"] = f"{land_sea_percentage} %"

        if len(file_name) > 0:
            out = output_folder / file_name
            ds.to_netcdf(out)
            return out
        else:
            return ds
    else:
        logging.warning("Project was not found.")
        raise RuntimeError()
