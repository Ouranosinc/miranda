import logging.config
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import dask.config
import numpy as np
import regionmask
import xarray as xr
from clisops.core import subset
from dask.distributed import Client
from xclim.indices import tas

from miranda.scripting import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)


dask.config.set(local_directory=f"{Path(__file__).parent}/dask_workers/")


# Needed pre-processing function
def _drop_those_time_bnds(dataset: xr.Dataset):
    if "time_bnds" in dataset.variables:
        return dataset.drop_vars(["time_bnds"])
    return dataset


outfiles = Path("/path/to/outfiles")
variables = ["tas", "pr"]
start, end = "1950", "2021"
domains = {"QC", "CAN", "AMNO"}


def processing(
    data: Dict[str, List[str]],
    domains: str,
    variables: Sequence[str],
    start: Optional[str] = None,
    end: Optional[str] = None,
    land_sea_mask: Optional[Dict] = None,
    land_sea_percentage: int = 50,
):  # noqa
    for domain in domains:
        output_folder = outfiles.joinpath(domain)
        output_folder.mkdir(exist_ok=True)
        for project, infiles in data.items():
            logging.info(f"Processing {project} data for domain {domain}.")
            for var in variables:
                fname = Path(f"{project}_{var}.nc")
                output = output_folder.joinpath(fname)
                if output.exists():
                    logging.info("File `%s` exists. Continuing..." % fname.name)
                    continue

                # Select only for variable of interest
                multifiles = [x for x in infiles if var in str(x)]
                with Client():
                    logging.info("Resampling variable `%s`." % var)
                    subset_geo = True
                    if domain == "QC":
                        lons = np.array([-79.76, -57.10])
                        lats = np.array([44.99, 62.59])
                    elif domain == "CAN":
                        lons = np.array([-141.02, -52.60])
                        lats = np.array([41.68, 83.14])
                    elif domain == "AMNO":
                        lons = np.array([-180, -13.5])
                        lats = np.array([16, 83.14])
                    else:
                        subset_geo = False
                        if start and end:
                            subset_time = True
                        else:
                            subset_time = False

                    if subset_geo:
                        ds = subset.subset_bbox(
                            xr.open_mfdataset(
                                multifiles,
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
                                multifiles,
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
                            multifiles,
                            combine="by_coords",
                            concat_dim=["time"],
                            chunks={"time": "auto"},
                            engine="netcdf4",
                            preprocess=_drop_those_time_bnds,
                        )

                    if project == "era5-land" and var == "pr":
                        ds["time"] = ds.time - np.timedelta64(1, "h")
                    try:
                        ds = ds.rename({"t2m": "tas"})
                    except ValueError:
                        pass
                    try:
                        ds = ds.rename({"tp": "pr"})
                    except ValueError:
                        pass
                    try:
                        ds = ds.rename({"prbc": "pr"})
                    except ValueError:
                        pass

                    if (
                        project.startswith("era5")
                        and project not in land_sea_mask.keys()
                    ):
                        try:
                            ds = ds.rename({"longitude": "lon", "latitude": "lat"})
                        except ValueError:
                            raise

                    dsOut = xr.Dataset()

                    if project in land_sea_mask.keys():
                        logging.info(
                            f"Masking variable with land-sea mask at {land_sea_percentage} % cutoff."
                        )
                        land_sea_mask, lsm_file = land_sea_mask[project]
                        lsm = subset.subset_bbox(
                            xr.open_dataset(lsm_file),
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
                        dsOut.attrs["land_sea_cutoff"] = f"{land_sea_percentage} %"

                    logging.info("Creating daily upscaled reanalyses.")
                    if var == "tas":
                        if project in ["era5-land", "cfsr", "nrcan"]:
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
                                output = output_folder.joinpath(
                                    f"{output.stem.split('_')[0]}_{v}.nc"
                                )
                                if output.exists():
                                    logging.info(
                                        "File `%s` exists. Continuing..." % output.name
                                    )
                                    continue
                                if v == "tas" and project == "nrcan":
                                    dsOut[v] = tas(tasmax=ds.tasmax, tasmin=ds.tasmin)
                                else:
                                    # Thanks for the help, xclim contributors
                                    r = ds[v_desired].resample(
                                        time="D", keep_attrs=True
                                    )
                                    dsOut[v] = getattr(r, func)(
                                        dim="time", keep_attrs=True
                                    )

                                logging.info(
                                    "Masking individual variable with AR6 regions."
                                )
                                mask = regionmask.defined_regions.ar6.all.mask(
                                    dsOut.lon, dsOut.lat
                                )
                                dsOut = dsOut.assign_coords(region=mask)

                                logging.info("Writing out %s." % output)
                                dsOut.to_netcdf(output)
                                del dsOut[v]
                            continue
                        elif project in ["era5", "merra2", "wfdei"]:
                            dsOut["tasmax"] = ds.tas.resample(time="D").max(
                                dim="time", keep_attrs=True
                            )
                            dsOut["tasmin"] = ds.tas.resample(time="D").min(
                                dim="time", keep_attrs=True
                            )
                            dsOut["tas"] = ds.tas.resample(time="D").mean(
                                dim="time", keep_attrs=True
                            )
                        else:
                            raise ValueError()

                    if var == "pr":
                        logging.info("Converting precipitation units")
                        if project == "era5-land":
                            dsOut["pr"] = (
                                ds.pr.resample(time="D").max(
                                    dim="time", keep_attrs=True
                                )
                                * 1000
                            )
                            dsOut["pr"].attrs["units"] = "mm d-1"
                        elif project == "era5":
                            dsOut["pr"] = (
                                ds.pr.resample(time="D").sum(
                                    dim="time", keep_attrs=True, skipna=False
                                )
                                * 1000
                            )
                            dsOut["pr"].attrs["units"] = "mm d-1"
                        elif project in ["merra2", "cfsr", "nrcan", "wfdei"]:
                            dsOut["pr"] = (
                                ds.pr.resample(time="D").mean(
                                    dim="time", keep_attrs=True
                                )
                                * 86400
                            )
                            dsOut["pr"].attrs["units"] = "mm d-1"

                    logging.info("Masking outputs with AR6 regions.")
                    mask = regionmask.defined_regions.ar6.all.mask(dsOut.lon, dsOut.lat)
                    dsOut = dsOut.assign_coords(region=mask)

                    logging.info("Writing out %s." % output)
                    dsOut.to_netcdf(output)
