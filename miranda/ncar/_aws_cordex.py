from __future__ import annotations

import ast
import functools
import json
import logging.config
import warnings
from json.decoder import JSONDecodeError
from pathlib import Path

import schema
import xarray as xr
from dask.diagnostics import ProgressBar
from xclim.core import calendar as xcal  # noqa

from miranda.gis import subsetting_domains
from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)

try:
    import intake  # noqa
    import intake_esm  # noqa
    import numcodecs  # noqa
    import s3fs  # noqa
except ImportError:
    warnings.warn(
        f"{__name__} functions require additional dependencies. "
        "Please install them with `pip install miranda[remote]`."
    )


__all__ = [
    "cordex_aws_calendar_correction",
    "cordex_aws_download",
]

_allowed_args = schema.Schema(
    {
        schema.Optional("variable"): schema.Schema(
            schema.Or(
                [
                    "hurs",
                    "huss",
                    "pr",
                    "prec",
                    "ps",
                    "rsds",
                    "sfcWind",
                    "tas",
                    "tasmax",
                    "tasmin",
                    "temp",
                    "tmax",
                    "tmin",
                    "uas",
                    "vas",
                ]
            )
        ),
        schema.Optional("frequency"): "day",
        schema.Optional("scenario"): schema.Schema(
            schema.Or(["eval", "hist", "rcp45", "rcp85", "hist-rcp45", "hist-rcp85"])
        ),
        schema.Optional("grid"): schema.Or(["NAM-22i", "NAM-44i"]),
        schema.Optional("bias_correction"): schema.Or(
            ["raw", "mbcn-Daymet", "mbcn-gridMET"]
        ),
    }
)


def cordex_aws_calendar_correction(ds) -> xr.Dataset | None:
    """AWS-stored CORDEX datasets are all on the same standard calendar, this converts
    the data back to the original calendar, removing added NaNs.

    Credit: Pascal Bourgault (@aulemahal)
    """
    orig_calendar = ds.attrs.get("original_calendar", "standard")

    if orig_calendar in ["365_day", "360_day"]:
        logging.info(f"Converting calendar to {orig_calendar}")
        ds = xcal.convert_calendar(ds, "noleap")  # drops Feb 29th
        if orig_calendar == "360_day":
            time = xcal.date_range_like(ds.time, calendar="360_day")
            ds = ds.where(
                ~ds.time.dt.dayofyear.isin([31, 90, 151, 243, 304]), drop=True
            )
            if ds.time.size != time.size:
                raise ValueError("Conversion of dataset to 360_day calendar failed.")
            ds["time"] = time

    return ds


def cordex_aws_download(
    target_folder: str | Path,
    *,
    search: dict[str, str | list[str]],
    correct_times: bool = False,
    domain: str | None = None,
):
    """Download CORDEX interpolated grid for North America from Amazon S3."""

    def _subset_preprocess(d: xr.Dataset, dom: list[float]) -> xr.Dataset:
        try:
            from clisops.core import subset_bbox

            n, w, s, e = subsetting_domains(dom)
            return subset_bbox(d, lon_bnds=[w, e], lat_bnds=[s, n])
        except ModuleNotFoundError:
            log_msg = (
                "This function requires the `clisops` library which is not installed. "
                "Domain subsetting step will be skipped."
            )
            warnings.warn(log_msg)
            return d

    schema.Schema(_allowed_args).validate(search)

    # Define the catalog description file location.
    catalog_url = (
        "https://ncar-na-cordex.s3-us-west-2.amazonaws.com/catalogs/aws-na-cordex.json"
    )

    # Interpret the "na-cordex-models" column as a list of values.
    col = intake.open_esm_datastore(
        catalog_url,
        csv_kwargs={"converters": {"na-cordex-models": ast.literal_eval}},
    )

    col_subset = col.search(**search)

    additional_kwargs = dict()
    if domain:
        additional_kwargs["preprocess"] = functools.partial(
            _subset_preprocess, dom=domain
        )

    dsets = col_subset.to_dataset_dict(
        zarr_kwargs={"consolidated": True},
        storage_options={"anon": True},
        **additional_kwargs,
    )
    logging.info(f"\nDataset dictionary keys:\n {dsets.keys()}")

    dds = list()
    for key in list(dsets.keys()):
        logging.info(f"Adding {key} to the search criteria.")
        dds.append(dsets[key])

    with ProgressBar():
        for ds in dds:
            scen = ds.attrs["experiment_id"]
            grid = str(ds.attrs["intake_esm_dataset_key"]).split(".")[-2]
            bias_correction = str(ds.attrs["intake_esm_dataset_key"]).split(".")[-1]
            attrs_ref = ds.attrs.copy()

            for i, member in enumerate(ds.member_id):
                for var in ds.variables:
                    if var in search["variable"]:
                        var_out = var
                        break

                new_attrs = dict()
                for key, vals in attrs_ref.items():
                    try:
                        mapped = json.loads(vals)
                        if isinstance(mapped, dict):
                            new_attrs[key] = mapped.get(str(member.values), "")
                    except JSONDecodeError:
                        new_attrs[key] = vals
                    except TypeError:
                        if len(vals) == 1:
                            new_attrs[key] = vals[0]
                        else:
                            raise RuntimeError()

                if correct_times:
                    try:
                        ds = cordex_aws_calendar_correction(ds)
                    except ValueError:
                        logging.warning(
                            f"Calendar failed to convert for {member.values} and variable {var_out}. Skipping..."
                        )
                        continue

                years, datasets = zip(*ds.isel(member_id=i).groupby("time.year"))

                for d in datasets:
                    d.attrs.update(new_attrs)

                out_folder = target_folder.joinpath(f"{member.values}_{scen}")
                out_folder.mkdir(exist_ok=True)

                file_name_pattern = (
                    f"{var_out}_{member.values}_day_{scen}_{grid}_{bias_correction}"
                )

                logging.info(f"Writing out files for {file_name_pattern}.")
                paths = [
                    out_folder.joinpath(f"{file_name_pattern}_{y}.nc")
                    for y in years
                    if not out_folder.joinpath(f"{file_name_pattern}_{y}.nc").exists()
                ]

                datasets = [
                    d
                    for y, d in zip(years, datasets)
                    if not out_folder.joinpath(f"{file_name_pattern}_{y}.nc").exists()
                ]

                if len(datasets) == 0:
                    logging.warning(
                        f"All files currently exist for {scen} and {member.name}. Continuing..."
                    )
                    continue

                logging.info(f"Final count of files: {len(datasets)}")

                xr.save_mfdataset(
                    datasets, paths, engine="h5netcdf", format="NETCDF4_CLASSIC"
                )
