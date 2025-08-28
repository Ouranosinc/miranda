"""Module to convert station data to Zarr format."""

from __future__ import annotations
import datetime as dt
import logging
import multiprocessing as mp
import os
import shutil
from collections.abc import Generator
from pathlib import Path
from zoneinfo import ZoneInfo

from miranda.convert._data_corrections import (
    dataset_conversion,
    load_json_data_mappings,
)
from miranda.convert.utils import (
    get_station_meta,
    make_monotonous_time,
    prj_dict,
    q_flag_dict,
    write_zarr,
)
from miranda.eccc._homogenized import create_canhomt_xarray
from miranda.ghcn import create_ghcn_xarray


logger = logging.getLogger("miranda.convert.stationdata")


def convert_statdata_bychunks(
    project: str,
    working_folder: str | os.PathLike[str] | None = None,
    cfvariable_list: list | None = None,
    start_year: int | None = None,
    end_year: int | None = None,
    lon_bnds: list[float] | None = None,
    lat_bnds: list[float] | None = None,
    n_workers: int = 4,
    n_stations: int = 100,
    update_from_raw: bool = False,
) -> None:
    """
    Convert GHCN or CanHomT station data to Zarr format.

    Requires GIS libraries (geopandas).

    Parameters
    ----------
    project : str
        Project name.
    working_folder : str or os.PathLike[str], optional
        The working folder. The default (None) is to use the current working directory.
    cfvariable_list : list, optional
        List of CF variable names. Optional.
    start_year : int, optional
        Start year. Optional.
    end_year : int, optional
        End year. Optional.
    lon_bnds : list of float, optional
        Longitude boundaries.
    lat_bnds : list of float, optional
        Latitude boundaries.
    n_workers : int
        Number of workers to use. Default is 4.
    n_stations : int
        Number of stations to process. Default is 100.
    update_from_raw : bool
        Whether to update from raw data.
    """
    try:
        import geopandas as gpd
        from shapely.geometry import box
    except ModuleNotFoundError as err:
        msg = "GNCN conversion requires the GIS libraries. Install them with `$ pip install miranda[gis]`."
        logger.error(msg)
        raise ModuleNotFoundError(msg) from err

    var_attrs = load_json_data_mappings(project=project)["variables"]
    if cfvariable_list:
        var_attrs = {v: var_attrs[v] for v in var_attrs if var_attrs[v]["_cf_variable_name"] in cfvariable_list}
    readme_url = None
    station_df = get_station_meta(project=project, lon_bnds=lon_bnds, lat_bnds=lat_bnds)
    if project == "ghcnd":
        readme_url = "https://noaa-ghcn-pds.s3.amazonaws.com/readme.txt"
        out_chunks = dict(time=(365 * 4) + 1, station=n_stations)
    # TODO ghcnh not implemented yet
    # elif project == "ghcnh":
    #     readme_url = "https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly/doc/ghcnh_DOCUMENTATION.pdf"
    #     out_chunks = dict(time=(365 * 4) + 1, station=n_stations)
    # logger.info("ghcnh not implemented yet")
    # exit()
    elif project == "canhomt_dly":
        out_chunks = dict(time=(365 * 4) + 1, station=n_stations)
    else:
        msg = f"Unknown project {project}"
        raise ValueError(msg)

    tz_file = Path(__file__).parent.joinpath("data/timezones-with-oceans-now.shapefile.zip")

    tz = gpd.read_file(tz_file).to_crs(epsg=4326)
    # clip to bbox for faster sjoin
    # Create a custom polygon
    polygon = box(lon_bnds[0] - 0.1, lat_bnds[0] - 0.1, lon_bnds[-1] + 0.1, lat_bnds[-1] + 0.1)
    poly_clip = gpd.GeoDataFrame([1], geometry=[polygon], crs=tz.crs)
    tz = tz.clip(poly_clip)
    station_df = gpd.GeoDataFrame(
        station_df,
        geometry=gpd.points_from_xy(station_df.lon, station_df.lat),
        crs=tz.crs,
    ).sjoin(tz, how="left")
    station_df = station_df.rename(columns={"tzid": "timezone"})
    if isinstance(working_folder, str):
        working_folder = Path(working_folder).expanduser()
    working_folder.mkdir(parents=True, exist_ok=True)
    working_folder.joinpath("zarr").mkdir(exist_ok=True)

    if end_year is not None:
        end_date = dt.datetime(end_year, 12, 31, 23, 59, 59, tzinfo=ZoneInfo("UTC"))
    else:
        end_date = dt.datetime.now().astimezone(ZoneInfo("UTC"))
    start_date = dt.datetime(start_year, 1, 1, 0, 0, 0, tzinfo=ZoneInfo("UTC"))

    if update_from_raw:
        for folder in working_folder.joinpath("zarr").iterdir():
            msg = f"Deleting {folder}"
            logger.info(msg)
            shutil.rmtree(folder)

    def _chunk_list(lst: list, n: int) -> Generator[list]:
        """Split list into chunks of size n."""
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    treated = []
    file_list = sorted(list(working_folder.joinpath("raw").rglob(f"*{prj_dict[project]['filetype']}")))
    file_list = [f for f in file_list if f.stem in station_df.station_id.tolist()]
    jobs = []
    for ii, ss in enumerate(_chunk_list(file_list, n_stations)):
        if ii not in treated:
            var_attrs_new = {}
            for vv, meta in var_attrs.items():
                cf_var = var_attrs[vv]["_cf_variable_name"]
                outzarr = working_folder.joinpath("zarr", cf_var, f"{project}_{ii}.zarr")
                if not outzarr.exists() or update_from_raw:
                    var_attrs_new[vv] = meta

            if var_attrs_new:
                dsall_vars = None
                if "ghcn" in project:
                    dsall_vars = create_ghcn_xarray(
                        in_files=ss,
                        variable_meta=var_attrs_new,
                        station_meta=station_df,
                        project=project,
                    )
                elif "canhomt" in project:
                    dsall_vars = create_canhomt_xarray(
                        in_files=ss,
                        variable_meta=var_attrs_new,
                        station_meta=station_df,
                        project=project,
                    )

                if dsall_vars is None:
                    continue
                dsall_vars = dsall_vars.sel(time=slice(str(start_date.year), str(end_date.year)))
                for kk, vv in var_attrs_new.items():
                    cf_var = var_attrs[kk]["_cf_variable_name"]
                    outzarr = working_folder.joinpath("zarr", cf_var, f"{project}_{ii}.zarr")
                    if kk not in dsall_vars.data_vars:
                        continue
                    dsout = dsall_vars.drop_vars([v for v in dsall_vars.data_vars if not v.startswith(kk)])
                    allnull_stat = dsout[kk].isnull().sum(dim="time") == len(dsout.time)
                    dsout = dsout.sel(station=~allnull_stat)
                    dsout = make_monotonous_time(dsout, freq=prj_dict[project]["freq"])

                    ds_corr = dataset_conversion(
                        dsout,
                        project=project,
                        add_version_hashes=False,
                        overwrite=update_from_raw,
                    )
                    ds_corr = ds_corr.rename({f"{kk}_q_flag": f"{cf_var}_q_flag"})
                    for vv in ds_corr.data_vars:
                        if ds_corr[vv].dtype == "float64":
                            ds_corr[vv] = ds_corr[vv].astype("float32")

                    desc_str = "; ".join([f"{k}:{v}" for k, v in q_flag_dict[project].items()])
                    if readme_url is not None:
                        desc_str = f"{desc_str}. See the readme file for information of quality flag (QFLAG1) codes : {readme_url}"
                    attrs = {
                        "flag_values": [c for c in q_flag_dict[project].keys()],
                        "flag_meanings": [c for c in q_flag_dict[project].values()],
                        "standard_name": f"{ds_corr[cf_var].attrs['standard_name']}_quality_flag",
                        "long_name": f"Quality flag for {cf_var}",
                        "description": desc_str,
                    }

                    ds_corr[f"{cf_var}_q_flag"].attrs = attrs

                    jobs.append((ds_corr, outzarr, out_chunks))
                    if len(jobs) >= n_workers:
                        pool = mp.Pool(n_workers)
                        pool.starmap(write_zarr, jobs)
                        pool.close()
                        pool.join()
                        jobs = []

            if len(jobs) > 0:
                pool = mp.Pool(n_workers)
                pool.starmap(write_zarr, jobs)
                pool.close()
                pool.join()
                jobs = []
            treated.append(ii)
