"""Dataset corrections submodule."""
from __future__ import annotations

import datetime
import os
from functools import partial
from pathlib import Path
from typing import Callable, Iterator, Sequence

import xarray as xr

from miranda.convert import (
    dimensions_compliance,
    metadata_conversion,
    threshold_mask,
    variable_conversion,
)
from miranda.convert._data_definitions import load_json_data_mappings
from miranda.convert._treatments import (
    cf_units_conversion,
    clip_values,
    conservative_regrid,
    correct_unit_names,
    ensure_correct_time_frequency,
    invert_value_sign,
    offset_time_dimension,
    preprocessing_corrections,
    transform_values,
)
from miranda.convert.utils import find_version_hash
from miranda.gis import subset_domain


def dataset_corrections(ds: xr.Dataset, project: str) -> xr.Dataset:
    """Convert variables to CF-compliant format"""
    metadata_definition = load_json_data_mappings(project)

    ds = correct_unit_names(ds, project, metadata_definition)
    ds = transform_values(ds, project, metadata_definition)
    ds = invert_value_sign(ds, project, metadata_definition)
    ds = cf_units_conversion(ds, metadata_definition)
    ds = clip_values(ds, project, metadata_definition)

    ds = dimensions_compliance(ds, project, metadata_definition)
    ds = ensure_correct_time_frequency(ds, project, metadata_definition)
    ds = offset_time_dimension(ds, project, metadata_definition)

    ds = variable_conversion(ds, project, metadata_definition)

    ds = metadata_conversion(ds, project, metadata_definition)

    ds.attrs["history"] = (
        f"{datetime.datetime.now()}: "
        f"Variables converted from original files using miranda.convert.{dataset_corrections.__name__}. "
        f"{ds.attrs.get('history')}".strip()
    )

    return ds


def dataset_conversion(
    input_files: (
        str
        | os.PathLike
        | Sequence[str | os.PathLike]
        | Iterator[os.PathLike]
        | xr.Dataset
    ),
    project: str,
    domain: str | None = None,
    mask: xr.Dataset | xr.DataArray | None = None,
    mask_cutoff: float | bool = False,
    regrid: bool = False,
    add_version_hashes: bool = True,
    preprocess: Callable | str | None = "auto",
    **xr_kwargs,
) -> xr.Dataset | xr.DataArray:
    """Convert an existing Xarray-compatible dataset to another format with variable corrections applied.

    Parameters
    ----------
    input_files : str or os.PathLike or Sequence[str or os.PathLike] or Iterator[os.PathLike] or xr.Dataset
        Files or objects to be converted.
        If sent a list or GeneratorType, will open with :py:func:`xarray.open_mfdataset` and concatenate files.
    project : {"cordex", "cmip5", "cmip6", "ets-grnch", "isimip-ft", "pcic-candcs-u6", "converted"}
        Project name for decoding/handling purposes.
    domain: {"global", "nam", "can", "qc", "mtl"}, optional
        Domain to perform subsetting for. Default: None.
    mask : Optional[Union[xr.Dataset, xr.DataArray]]
        DataArray or single data_variable dataset containing mask.
    mask_cutoff : float or bool
        If land_sea_mask supplied, the threshold above which to mask with land_sea_mask. Default: False.
    regrid : bool
        Performing regridding with xesmf. Default: False.
    add_version_hashes : bool
        If True, version name and sha256sum of source file(s) will be added as a field among the global attributes.
    preprocess : callable or str, optional
        Preprocessing functions to perform over each Dataset.
        Default: "auto" - Run preprocessing fixes based on supplied fields from metadata definition.
        Callable - Runs function over Dataset (single) or supplied to `preprocess` (multifile dataset).
    **xr_kwargs
        Arguments passed directly to xarray.

    Returns
    -------
    xr.Dataset or xr.DataArray
    """
    if isinstance(input_files, xr.Dataset):
        ds = input_files
    else:
        if isinstance(input_files, (str, os.PathLike)):
            if Path(input_files).is_dir():
                files = []
                files.extend([f for f in Path(input_files).glob("*.nc")])
                files.extend([f for f in Path(input_files).glob("*.zarr")])
            else:
                files = [Path(input_files)]
        elif isinstance(input_files, (Sequence, Iterator)):
            files = [Path(f) for f in input_files]
        else:
            files = input_files
        version_hashes = dict()
        if add_version_hashes:
            for file in files:
                version_hashes[file.name] = find_version_hash(file)

        preprocess_kwargs = dict()
        if preprocess:
            if preprocess == "auto":
                preprocess_kwargs.update(
                    preprocess=partial(preprocessing_corrections, project=project)
                )
            elif isinstance(preprocess, Callable):
                preprocess_kwargs.update(preprocess=preprocess)

        if len(files) == 1:
            ds = xr.open_dataset(files[0], **xr_kwargs)
            for _, process in preprocess_kwargs.items():
                ds = process(ds)
        else:
            ds = xr.open_mfdataset(files, **xr_kwargs, **preprocess_kwargs)
        if version_hashes:
            ds.attrs.update(dict(original_files=str(version_hashes)))

    ds = dataset_corrections(ds, project)

    if domain:
        ds = subset_domain(ds, domain)

    if isinstance(mask, (str, Path)):
        mask = xr.open_dataset(mask)
    if isinstance(mask, (xr.Dataset, xr.DataArray)):
        if regrid:
            mask = conservative_regrid(ds, mask)
        ds = threshold_mask(ds, mask=mask, mask_cutoff=mask_cutoff)

    return ds
