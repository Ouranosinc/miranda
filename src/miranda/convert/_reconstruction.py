from __future__ import annotations
import logging
import os
from collections.abc import Sequence
from pathlib import Path

import dask.config
import xarray as xr
from dask import compute
from dask.diagnostics import ProgressBar
from xclim.core import calendar

from miranda.gis import subset_domain
from miranda.io.utils import delayed_write, get_chunks_on_disk
from miranda.utils import chunk_iterables
from miranda.vocabularies import project_institute, xarray_frequencies_to_cmip6like

from ._aggregation import aggregate as aggregate_func
from .corrections import dataset_corrections


logger = logging.getLogger("miranda.convert.reconstruction")


dask.config.set(local_directory=f"{Path(__file__).parent}/dask_workers/")


__all__ = ["reanalysis_processing"]


# Needed pre-processing function
def _drop_those_time_bnds(dataset: xr.Dataset) -> xr.Dataset:
    if "time_bnds" in dataset.variables:
        return dataset.drop_vars(["time_bnds"])
    return dataset


def reanalysis_processing(
    data: dict[str, list[str | os.PathLike]],
    output_folder: str | os.PathLike,
    variables: Sequence[str],
    aggregate: str | bool = False,
    domains: str | list[str] = "_DEFAULT",
    start: str | None = None,
    end: str | None = None,
    target_chunks: dict | None = None,
    output_format: str = "netcdf",
    overwrite: bool = False,
    engine: str = "h5netcdf",
    n_workers: int = 4,
    **dask_kwargs,
) -> None:
    """
    Reanalysis processing.

    Parameters
    ----------
    data : dict[str, list[str]]
    output_folder : str or os.PathLike
    variables : Sequence[str]
    aggregate : {"day", None}
    domains : {"QC", "CAN", "AMNO", "NAM", "GLOBAL"}
    start : str, optional
    end : str, optional
    target_chunks : dict, optional
    output_format : {"netcdf", "zarr"}
    overwrite : bool
    engine : {"netcdf4", "h5netcdf"}
    n_workers : int

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

    with (
        ProgressBar(),
        dask.config.set(
            **{"array.slicing.split_large_chunks": False},
            n_workers=n_workers,
            **dask_kwargs,
        ),
    ):
        out_files = Path(output_folder)
        if isinstance(domains, str):
            domains = [domains]

        for domain in domains:
            if domain == "_DEFAULT":
                logger.warning("No domain specified. proceeding with 'not-specified'.")
                output_folder = output_folder
                domain = "not-specified"
            elif isinstance(domain, str):
                output_folder = out_files.joinpath(domain)  # noqa
            else:
                raise NotImplementedError()

            output_folder.mkdir(exist_ok=True)

            for project, in_files in data.items():
                msg = f"Processing {project} data{f' for domain {domain}' if domain != 'not_specified' else ''}."
                logger.info(msg)

                for var in variables:
                    # Select only for variable of interest
                    multi_files = sorted(str(x) for x in in_files if f"{var}_" in str(x))

                    if multi_files:
                        all_chunks = get_chunks_on_disk(multi_files[0])
                        chunks = all_chunks[var]

                        if target_chunks is None:
                            output_chunks = dict()
                            mappings = dict(longitude="lon", latitude="lat")
                            for k, v in chunks.items():
                                if k in mappings.keys():
                                    output_chunks[mappings[k]] = v
                                else:
                                    output_chunks[k] = v

                            msg = f"No `target_chunks` set. Proceeding with following found chunks: {output_chunks}."
                            logger.warning(msg)
                        else:
                            output_chunks = target_chunks

                        msg = f"Resampling variable `{var}`."
                        logger.info(msg)

                        if aggregate:
                            time_freq = aggregate
                        else:
                            parse_freq = calendar.parse_offset(xr.infer_freq(xr.open_dataset(multi_files[0]).time))
                            time_freq = f"{parse_freq[0]}{xarray_frequencies_to_cmip6like[parse_freq[1]]}"

                        institute = project_institute(project)
                        file_name = f"{var}_{time_freq}_{institute}_{project}"
                        if domain != "not-specified":
                            file_name = f"{file_name}_{domain}"
                        if not chunks:
                            chunks = dict(time=24 * 10, lon=50, lat=50)
                            print(chunks)
                        xr_kwargs = dict(
                            chunks=chunks,
                            engine=engine,
                            preprocess=_drop_those_time_bnds,
                            parallel=True,
                        )

                        # Subsetting operations
                        if domain.lower() in ["global", "not-specified"]:
                            if start or end:
                                ds = xr.open_mfdataset(multi_files, **xr_kwargs).sel(time=slice(start, end))
                            else:
                                ds = xr.open_mfdataset(multi_files, **xr_kwargs)
                        else:
                            ds = subset_domain(
                                xr.open_mfdataset(multi_files, **xr_kwargs),
                                domain,
                                start_date=start,
                                end_date=end,
                            )

                        ds.attrs.update(dict(frequency=time_freq, domain=domain))
                        ds = dataset_corrections(ds, project=project)

                        if time_freq.lower() == "day":
                            dataset = aggregate_func(ds, freq="day")
                            freq = "YS"
                        else:
                            out_variable = list(ds.data_vars)[0] if len(list(ds.data_vars)) == 1 else None
                            dataset = {out_variable: ds}
                            freq = "MS"

                        if len(dataset) == 0:
                            msg = f"Daily aggregation methods for variable `{var}` are not supported. Continuing..."
                            logger.warning(msg)

                        for key in dataset.keys():
                            ds = dataset[key]

                            # TODO: What do we do about multivariable files. Are they even allowed?
                            out_variable = list(ds.data_vars)[0] if len(list(ds.data_vars)) == 1 else None
                            file_name1 = file_name.replace(f"{var}_", f"{out_variable}_")

                            msg = f"Writing out fixed files for {file_name1}."
                            logger.info(msg)
                            years, datasets = zip(*ds.resample(time=freq), strict=False)
                            if freq == "MS":
                                format_str = "%Y-%m"
                                iterable_chunks = 36
                            else:
                                format_str = "%Y"
                                iterable_chunks = 10

                            out_filenames = [
                                output_folder.joinpath(f"{file_name1}_{xr.DataArray(year).dt.strftime(format_str).values}{suffix}") for year in years
                            ]

                            jobs = list()
                            if output_format != "zarr" and overwrite:
                                msg = f"Removing existing {output_format} files for {var}."

                                logger.warning(msg)
                            for i, d in enumerate(datasets):
                                if out_filenames[i].exists() and out_filenames[i].is_file() and overwrite:
                                    out_filenames[i].unlink()

                                if not out_filenames[i].exists() or (out_filenames[i].is_dir() and overwrite):
                                    jobs.append(
                                        delayed_write(
                                            d,
                                            out_filenames[i],
                                            output_format,
                                            overwrite,
                                            target_chunks=output_chunks,
                                        )
                                    )

                            if len(jobs) == 0:
                                msg = f"All output files for `{var}` currently exist. To overwrite them, set `overwrite=True`. Continuing..."
                                logger.warning(msg)
                            else:
                                chunked_jobs = chunk_iterables(jobs, iterable_chunks)
                                logger.info(msg)
                                iterations = 0
                                for chunk in chunked_jobs:
                                    iterations += 1
                                    msg = f"Processing iteration {iterations} for variable `{var}`."
                                    logger.info(msg)
                                    compute(chunk)
                    else:
                        msg = f"No files found for variable {var}."
                        logger.info(msg)
