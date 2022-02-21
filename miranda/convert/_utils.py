import datetime
import json
import logging.config
import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import dask.config
import netCDF4
import numpy as np
import regionmask
import xarray
import xarray as xr
import zarr
from clisops.core import subset
from dask import compute
from dask.diagnostics import ProgressBar
from xarray import Dataset
from xclim.core import calendar, units
from xclim.indices import tas

from miranda.gis.subset import subsetting_domains
from miranda.scripting import LOGGING_CONFIG
from miranda.utils import chunk_iterables

from ._data import project_institutes, xarray_frequencies_to_cmip6

logging.config.dictConfig(LOGGING_CONFIG)


dask.config.set(local_directory=f"{Path(__file__).parent}/dask_workers/")

LATLON_COORDINATE_PRECISION = dict()
LATLON_COORDINATE_PRECISION["era5-land"] = 4

VERSION = datetime.datetime.now().strftime("%Y.%m.%d")


# Needed pre-processing function
def _drop_those_time_bnds(dataset: xr.Dataset) -> xr.Dataset:
    if "time_bnds" in dataset.variables:
        return dataset.drop_vars(["time_bnds"])
    return dataset


def get_chunks_on_disk(ncfile: Union[os.PathLike, str]) -> dict:
    ds = netCDF4.Dataset(ncfile)
    chunks = dict()
    for v in ds.variables:
        chunks[v] = dict()
        for ii, dim in enumerate(ds[v].dimensions):
            chunks[v][dim] = ds[v].chunking()[ii]
    return chunks


def reanalysis_processing(
    data: Dict[str, List[Union[str, os.PathLike]]],
    output_folder: Union[str, os.PathLike],
    variables: Sequence[str],
    aggregate: Union[str, bool] = False,
    domains: Union[str, List[str]] = "_DEFAULT",
    start: Optional[str] = None,
    end: Optional[str] = None,
    target_chunks: Optional[dict] = None,
    output_format: str = "netcdf",
    overwrite: bool = False,
    engine: str = "h5netcdf"
    # n_workers: int = 4,
    # **dask_kwargs,
) -> None:
    """

    Parameters
    ----------
    data: Dict[str, List[str]]
    output_folder: Union[str, os.PathLike]
    variables: Sequence[str]
    aggregate: {"day", None}
    domains: {"QC", "CAN", "AMNO", "GLOBAL"}
    start: str, optional
    end: str, optional
    target_chunks: dict, optional
    output_format: {"netcdf", "zarr"}
    overwrite: bool
    engine: {"netcdf4", "h5netcdf"}

    Returns
    -------
    None
    """
    if output_format == "netcdf":
        suffix = ".nc"
    elif output_format == "zarr":
        suffix = ".zarr"
    else:
        raise NotImplementedError(f"`output_format`: '{output_format}")

    with ProgressBar(), dask.config.set(**{"array.slicing.split_large_chunks": False}):
        out_files = Path(output_folder)
        if isinstance(domains, str):
            domains = [domains]

        for domain in domains:
            if domain == "_DEFAULT":
                logging.warning("No domain specified. proceeding with 'not-specified'")
                output_folder = output_folder
                domain = "not-specified"
            elif isinstance(domain, str):
                output_folder = out_files.joinpath(domain)  # noqa
            else:
                raise NotImplementedError()

            output_folder.mkdir(exist_ok=True)

            for project, in_files in data.items():
                if domain != "not_specified":
                    logging.info(f"Processing {project} data for domain {domain}.")
                else:
                    logging.info(f"Processing {project} data.")
                for var in variables:
                    # Select only for variable of interest
                    multi_files = sorted(x for x in in_files if f"{var}_" in str(x))

                    if multi_files:
                        chunks = get_chunks_on_disk(multi_files[0])
                        chunks = chunks[var]

                        if target_chunks is None:
                            output_chunks = dict()
                            mappings = dict(longitude="lon", latitude="lat")
                            for k, v in chunks.items():
                                if k in mappings.keys():
                                    output_chunks[mappings[k]] = v
                                else:
                                    output_chunks[k] = v

                            logging.warning(
                                "No target_chunks set."
                                f" Proceeding with following found chunks: {output_chunks}"
                            )
                        else:
                            output_chunks = target_chunks

                        logging.info(f"Resampling variable `{var}`.")

                        if aggregate:
                            time_freq = aggregate
                        else:
                            parse_freq = calendar.parse_offset(
                                xr.infer_freq(xr.open_dataset(multi_files[0]).time)
                            )
                            time_freq = f"{parse_freq[0]}{xarray_frequencies_to_cmip6[parse_freq[1]]}"

                        institute = project_institutes[project]
                        file_name = "_".join([var, time_freq, institute, project])
                        if domain != "not-specified":
                            file_name = f"{file_name}_{domain}"

                        xr_kwargs = dict(
                            chunks=chunks,
                            engine=engine,
                            preprocess=_drop_those_time_bnds,
                            parallel=True,
                        )

                        # Subsetting operations
                        if domain.lower() in ["global", "not-specified"]:
                            if start or end:
                                ds = subset.subset_time(
                                    xr.open_mfdataset(multi_files, **xr_kwargs),
                                    start_date=start,
                                    end_date=end,
                                )
                            else:
                                ds = xr.open_mfdataset(multi_files, **xr_kwargs)
                        else:
                            region = subsetting_domains(domain)
                            lon_values = np.array([region[1], region[3]])
                            lat_values = np.array([region[0], region[2]])

                            ds = subset.subset_bbox(
                                xr.open_mfdataset(multi_files, **xr_kwargs),
                                lon_bnds=lon_values,
                                lat_bnds=lat_values,
                                start_date=start,
                                end_date=end,
                            )

                        ds.attrs.update(dict(frequency=time_freq, domain=domain))
                        ds = variable_conversion(
                            ds, project=project, output_format=output_format
                        )

                        if time_freq.lower() == "day":
                            dataset = daily_aggregation(ds)
                            freq = "YS"
                        else:
                            out_variable = (
                                list(ds.data_vars)[0]
                                if len(list(ds.data_vars)) == 1
                                else None
                            )
                            dataset = {out_variable: ds}
                            freq = "MS"

                        if len(dataset) == 0:
                            logging.warning(
                                f"Daily aggregation methods for variable `{var}` are not supported. "
                                "Continuing..."
                            )

                        for key in dataset.keys():
                            ds = dataset[key]
                            out_variable = (
                                list(ds.data_vars)[0]
                                if len(list(ds.data_vars)) == 1
                                else None
                            )
                            file_name1 = file_name.replace(
                                f"{var}_", f"{out_variable}_"
                            )
                            logging.info(f"Writing out fixed files for {file_name1}.")

                            years, datasets = zip(*ds.resample(time=freq))

                            if freq == "MS":
                                format_str = "%Y-%m"
                            else:
                                format_str = "%Y"

                            out_filenames = [
                                output_folder.joinpath(
                                    f"{file_name1}_{xr.DataArray(year).dt.strftime(format_str).values}{suffix}"
                                )
                                for year in years
                            ]

                            jobs = list()

                            if overwrite:
                                logging.warning(
                                    f"Removing existing {output_format} files for {var}."
                                )
                            for i, d in enumerate(datasets):
                                if out_filenames[i].exists() and overwrite:
                                    if (
                                        out_filenames[i].is_dir()
                                        and output_format == "zarr"
                                    ):
                                        shutil.rmtree(out_filenames[i])
                                    else:
                                        out_filenames[i].unlink()

                                if not out_filenames[i].exists():
                                    jobs.append(
                                        delayed_write(
                                            d,
                                            out_filenames[i],
                                            output_chunks,
                                            output_format,
                                        )
                                    )
                            if len(jobs) == 0:
                                logging.warning(
                                    f"All output files for `{var}` currently exist."
                                    " To overwrite them, set `overwrite=True`. Continuing..."
                                )
                            else:
                                chunked_jobs = chunk_iterables(jobs, 10)
                                for chunk in chunked_jobs:
                                    compute(chunk)
                    else:
                        logging.info(f"No files found for variable {var}.")


def delayed_write(
    ds: xarray.Dataset, outfile: Path, target_chunks: dict, output_format: str
):
    # Set correct chunks in encoding options
    kwargs = dict()
    kwargs["encoding"] = dict()
    for name, da in ds.data_vars.items():
        chunks = list()
        for dim in da.dims:
            if dim in target_chunks.keys():
                chunks.append(target_chunks[str(dim)])
            else:
                chunks.append(len(da[dim]))

        if output_format == "netcdf":
            kwargs["encoding"][name] = {
                "chunksizes": chunks,
                "zlib": True,
            }
            kwargs["compute"] = False
        elif output_format == "zarr":
            ds = ds.chunk(target_chunks)
            kwargs["encoding"][name] = {
                "chunks": chunks,
                "compressor": zarr.Blosc(),
            }
            kwargs["compute"] = False
    if kwargs["encoding"]:
        kwargs["encoding"]["time"] = {"dtype": "int32"}

    return getattr(ds, f"to_{output_format}")(outfile, **kwargs)


def variable_conversion(
    ds: xarray.Dataset, project: str, output_format: str
) -> xarray.Dataset:
    """Convert variables to CF-compliant format"""

    def _correct_units_names(d: xarray.Dataset, p: str):
        key = "_corrected_units"
        for vv in d.data_vars:
            if p in metadata_definition["variable_entry"][vv][key].keys():
                # ds_units = d[vv].attrs["units"]  # FIXME: Is this needed anywhere ?
                d[vv].attrs["units"] = metadata_definition["variable_entry"][vv][key][
                    project
                ]
        return d

    if project in ["era5", "era5-single-levels", "era5-land"]:
        metadata_definition = json.load(
            open(Path(__file__).parent.parent / "ecmwf" / "ecmwf_cf_attrs.json")
        )
    else:
        raise NotImplementedError()

    # for de-accumulation or conversion to flux
    def _transform(d: xarray.Dataset, p: str):
        key = "_transformation"
        d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
        for vv in d.data_vars:
            if p in metadata_definition["variable_entry"][vv][key].keys():
                if metadata_definition["variable_entry"][vv][key][p] == "deaccumulate":
                    freq = xr.infer_freq(ds.time)
                    try:
                        offset = (
                            float(calendar.parse_offset(freq)[0])
                            if calendar.parse_offset(freq)[0] != ""
                            else 1.0
                        )
                    except TypeError:
                        logging.error(
                            f"Unable to parse the time frequency for variable `{vv}`. "
                            "Verify data integrity before retrying."
                        )
                        raise

                    # accumulated hourly to hourly flux (de-accumulation)
                    with xr.set_options(keep_attrs=True):
                        out = d[vv].diff(dim="time")
                        out = d[vv].where(
                            d[vv].time.dt.hour == int(offset), out.broadcast_like(d[vv])
                        )
                        out = units.amount2rate(out)
                    d_out[out.name] = out
                elif metadata_definition["variable_entry"][vv][key][p] == "amount2rate":
                    out = units.amount2rate(
                        d[vv],
                        out_units=metadata_definition["variable_entry"][vv]["units"],
                    )
                    d_out[out.name] = out
                else:
                    raise NotImplementedError(
                        f"Unknown transformation "
                        f"{metadata_definition['variable_entry'][vv][key][p]}"
                    )

            else:
                d_out[vv] = d[vv]
        return d_out

    # For converting variable units to standard workflow units
    def _units_cf_conversion(d: xarray.Dataset) -> xarray.Dataset:
        descriptions = metadata_definition["variable_entry"]
        for v in d.data_vars:
            d[v] = units.convert_units_to(d[v], descriptions[v]["units"])
        return d

    # Add and update existing metadata fields
    def _metadata_conversion(d: xarray.Dataset, p: str, o: str) -> xarray.Dataset:

        # Add global attributes
        d.attrs.update(metadata_definition["Header"])
        d.attrs.update(dict(project=p, output_format=o))

        # Date-based versioning
        d.attrs.update(dict(version=VERSION))

        history = (
            f"{d.attrs['history']}\n[{datetime.datetime.now()}] Converted from original data to {o}"
            " with modified metadata for CF-like compliance."
        )
        d.attrs.update(dict(history=history))

        descriptions = metadata_definition["variable_entry"]

        # Add variable metadata
        for v in d.data_vars:
            descriptions[v].pop("_corrected_units")
            descriptions[v].pop("_transformation")
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
        sort_dims = []
        for orig, new in dict(longitude="lon", latitude="lat").items():
            try:

                d = d.rename({orig: new})
                if new == "lon" and np.any(d.lon > 180):
                    lon1 = d.lon.where(d.lon <= 180.0, d.lon - 360.0)
                    d[new] = lon1
                sort_dims.append(new)
            except KeyError:
                pass
            if project in LATLON_COORDINATE_PRECISION.keys():
                d[new] = d[new].round(LATLON_COORDINATE_PRECISION[project])
        if sort_dims:
            d = d.sortby(sort_dims)
        return d

    ds = _correct_units_names(ds, project)
    ds = _transform(ds, project)
    ds = _units_cf_conversion(ds)
    ds = _metadata_conversion(ds, project, output_format)
    ds = _dims_conversion(ds)

    return ds


def daily_aggregation(ds) -> Dict[str, Dataset]:
    logging.info("Creating daily upscaled reanalyses.")

    daily_dataset = dict()
    for variable in ds.data_vars:
        if variable in ["tas", "tdps"]:
            # Some looping to deal with memory consumption issues
            for v, func in {
                f"{variable}max": "max",
                f"{variable}min": "min",
                f"{variable}": "mean",
            }.items():
                ds_out = xr.Dataset()
                ds_out.attrs = ds.attrs.copy()
                ds_out.attrs["frequency"] = "day"

                method = (
                    f"time: {func}{'imum' if func != 'mean' else ''} (interval: 1 day)"
                )
                ds_out.attrs["cell_methods"] = method

                if v == "tas" and not hasattr(ds, "tas"):
                    ds_out[v] = tas(tasmax=ds.tasmax, tasmin=ds.tasmin)
                else:
                    # Thanks for the help, xclim contributors
                    r = ds[variable].resample(time="D")
                    ds_out[v] = getattr(r, func)(dim="time", keep_attrs=True)

                daily_dataset[v] = ds_out
                del ds_out

        elif variable in ["pr", "prsn", "snd", "snw"]:
            ds_out = xr.Dataset()
            ds_out.attrs = ds.attrs.copy()
            ds_out.attrs["frequency"] = "day"
            ds_out.attrs["cell_methods"] = "time: mean (interval: 1 day)"
            logging.info(f"Converting {variable} to daily time step (daily mean).")
            ds_out[variable] = (
                ds[variable].resample(time="D").mean(dim="time", keep_attrs=True)
            )

            daily_dataset[variable] = ds_out
            del ds_out
        else:
            continue

    return daily_dataset


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

        lon_bounds = np.array([ds.lon.min(), ds.lon.max()])
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
        return ds
    else:
        logging.warning("Project was not found.")
        raise RuntimeError()
