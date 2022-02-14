import json
import logging.config
import os
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

logging.config.dictConfig(LOGGING_CONFIG)


dask.config.set(local_directory=f"{Path(__file__).parent}/dask_workers/")

# map xarray freq to CMIP6 controlled vocabulary.
# see: https://github.com/WCRP-CMIP/CMIP6_CVs/blob/master/CMIP6_frequency.json
XR_FREQ_TO_CMIP6 = {
    "H": "hr",
    "D": "day",
    "M": "mon",
    "A": "yr",
    "Q": "qtr",  # TODO does this make sense? does not exist in cmip6 CV
}

PROJECT_INSTITUTES = {
    "cfsr": "ncar",
    "era5": "ecmwf",
    "era5-land": "ecmwf",
    "merra2": "nasa",
    "nrcan-gridded-10km": "nrcan",
    "sc-earth": "usask",
    "wfdei-gem-capa": "usask",
}
PROJECT_LONS_FROM_0TO360 = ["era5", "cfsr"]

LATLON_COORDINATE_TOLERANCE = dict()
LATLON_COORDINATE_TOLERANCE["era5-land"] = 4


# Needed pre-processing function
def _drop_those_time_bnds(dataset: xr.Dataset) -> xr.Dataset:
    if "time_bnds" in dataset.variables:
        return dataset.drop_vars(["time_bnds"])
    return dataset


def get_chunks_on_disk(ncfile: Union[os.PathLike, str]) -> dict:
    ds = netCDF4.Dataset(ncfile)
    chunks = {}
    for v in ds.variables:
        chunks[v] = {}
        for ii, dim in enumerate(ds[v].dimensions):
            chunks[v][dim] = ds[v].chunking()[ii]
    return chunks


def reanalysis_processing(
    data: Dict[str, List[Union[str, os.PathLike]]],
    output_folder: Union[str, os.PathLike],
    variables: Sequence[str],
    aggregate: Union[str, bool] = False,
    domains: Optional[Sequence[str]] = None,
    start: Optional[str] = None,
    end: Optional[str] = None,
    target_chunks: Optional[dict] = None,
    output_format: str = "netcdf",
) -> None:
    """

    Parameters
    ----------
    data: Dict[str, List[str]]
    output_folder: Union[str, os.PathLike]
    variables: Sequence[str]
    aggregate: {"day", None}
    domains: {"QC", "CAN", "AMNO", None}
    start: str, optional
    end: str, optional
    target_chunks: dict, optional
    output_format: {"netcdf", "zarr"}

    Returns
    -------
    None
    """
    with ProgressBar(), dask.config.set(**{"array.slicing.split_large_chunks": False}):
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
                    multi_files = sorted(x for x in in_files if f"{var}_" in str(x))
                    if multi_files:
                        chunks = get_chunks_on_disk(multi_files[0])
                        try:
                            chunks = chunks[var]
                        except KeyError:
                            chunks = chunks["sd"]  # era5 'sde' file has 'sd' variable??
                        logging.info("Resampling variable `%s`." % var)

                        if aggregate:
                            time_freq = aggregate
                        else:
                            parse_freq = calendar.parse_offset(
                                xr.infer_freq(xr.open_dataset(multi_files[0]).time)
                            )
                            time_freq = (
                                f"{parse_freq[0]}{XR_FREQ_TO_CMIP6[parse_freq[1]]}"
                            )

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
                                    chunks=chunks,
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
                                    chunks=chunks,
                                    engine="netcdf4",
                                    preprocess=_drop_those_time_bnds,
                                ),
                                start_date=start,
                                end_date=end,
                            )
                        elif not any([subset_time, domain, ds]):
                            ds = xr.open_mfdataset(
                                multi_files,
                                chunks=chunks,
                                engine="netcdf4",
                                preprocess=_drop_those_time_bnds,
                            )

                        ds = variable_conversion(ds, project=project)

                        if time_freq.lower() == "day":
                            dataset = daily_aggregation(
                                ds, project
                            )  # file_name, output_folder)
                            freq = "YS"
                        else:
                            outvar = (
                                list(ds.data_vars)[0]
                                if len(list(ds.data_vars)) == 1
                                else None
                            )
                            dataset = {outvar: ds}
                            freq = "MS"

                        for key in dataset.keys():
                            ds = dataset[key]
                            outvar = (
                                list(ds.data_vars)[0]
                                if len(list(ds.data_vars)) == 1
                                else None
                            )
                            file_name1 = file_name.replace(f"{var}_", f"{outvar}_")
                            logging.info("Writing out fixed files for %s." % file_name1)

                            years, datasets = zip(*ds.resample(time=freq))

                            if freq == "MS":
                                format_str = "%Y-%m"
                            else:
                                format_str = "%Y"
                            if output_format == "netcdf":
                                suffix = ".nc"
                            elif output_format == "zarr":
                                suffix = ".zarr"
                            out_filenames = [
                                output_folder.joinpath(
                                    f"{file_name1}_{xr.DataArray(year).dt.strftime(format_str).values}{suffix}"
                                )
                                for year in years
                            ]

                            jobs = []
                            for ii, d in enumerate(datasets):
                                jobs.append(
                                    delayed_write(
                                        d,
                                        out_filenames[ii],
                                        target_chunks,
                                        output_format,
                                    )
                                )
                            compute(jobs)
                    else:
                        logging.info(f"No files found for variable {var}.")


def delayed_write(
    ds: xarray.Dataset, outfile: Path, target_chunks: dict, output_format: str
):
    # Set correct chunks in encoding options
    kwargs = dict()
    kwargs["encoding"] = dict()
    for name, da in ds.data_vars.items():
        chunks = []
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


def variable_conversion(ds: xarray.Dataset, project: str) -> xarray.Dataset:
    """Convert variables to CF-compliant format"""

    # Add and update existing metadata fields
    def _metadata_conversion(d: xarray.Dataset, p: str) -> xarray.Dataset:

        # Add global attributes
        d.attrs.update(metadata_definition["Header"])
        d.attrs.update(dict(project=p))
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
            if project in LATLON_COORDINATE_TOLERANCE.keys():
                d[new] = d[new].round(LATLON_COORDINATE_TOLERANCE[project])
        if sort_dims:
            d = d.sortby(sort_dims)
        return d

    # For converting variable units to standard workflow units
    def _units_cf_conversion(d: xarray.Dataset, p: str) -> xarray.Dataset:
        descriptions = metadata_definition["variable_entry"]
        for v in d.data_vars:
            d[v] = units.convert_units_to(d[v], descriptions[v]["units"])
        return d

    # for de-accumulation or conversion to flux
    def _transform(d: xarray.Dataset, p: str):
        key = "_transformation"
        d_out = xr.Dataset(coords=d.coords, attrs=d.attrs)
        for vv in d.data_vars:
            if p in metadata_definition["variable_entry"][vv][key].keys():
                if metadata_definition["variable_entry"][vv][key][p] == "deaccumulate":
                    freq = xr.infer_freq(ds.time)
                    offset = (
                        float(calendar.parse_offset(freq)[0])
                        if calendar.parse_offset(freq)[0] != ""
                        else 1.0
                    )

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

    ds = _correct_units_names(ds, project)
    ds = _transform(ds, project)
    ds = _units_cf_conversion(ds, project)
    ds = _metadata_conversion(ds, project)
    ds = _dims_conversion(ds)

    return ds


def daily_aggregation(ds, project: str) -> Dict[str, Dataset]:
    # Daily variable aggregation operations
    # input_file = Path(input_file)
    # output_folder = Path(output_folder)
    logging.info("Creating daily upscaled reanalyses.")

    daily_dataset = dict()

    for variable in ds.data_vars:
        # input_file_parts = input_file.stem.split("_")
        # input_file_parts[0] = variable
        # input_file_parts[1] = "day"
        # output = output_folder.joinpath(f"{'_'.join(input_file_parts)}")

        # if any([f for f in output.glob("*")]):
        #     logging.info("Files for `%s` exist. Continuing..." % output.name)
        #     return

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

                # input_file_parts[0] = v
                # output = output_folder.joinpath(f"{'_'.join(input_file_parts)}")

                ds_out = xr.Dataset()
                ds_out.attrs = ds.attrs.copy()
                ds_out.attrs["frequency"] = "day"
                if v == "tas" and not hasattr(ds, "tas"):
                    ds_out[v] = tas(tasmax=ds.tasmax, tasmin=ds.tasmin)
                else:
                    # Thanks for the help, xclim contributors
                    r = ds[v_desired].resample(time="D", keep_attrs=True)
                    ds_out[v] = getattr(r, func)(dim="time", keep_attrs=True)

                daily_dataset[v] = ds_out
                del ds_out

        elif variable in ["pr", "prsn", "snd", "snw"]:
            ds_out = xr.Dataset()
            ds_out.attrs = ds.attrs.copy()
            ds_out.attrs["frequency"] = "day"
            logging.info(f"Converting {variable} to daily time step (daily mean)")
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