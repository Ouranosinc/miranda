import hashlib
import logging.config
import multiprocessing
import os
import shutil
import sys
from functools import partial
from pathlib import Path
from types import GeneratorType
from typing import List, Mapping, Optional, Union

from miranda import Decoder
from miranda.decode import guess_project
from miranda.scripting import LOGGING_CONFIG
from miranda.utils import discover_data
from miranda.validators import GRIDDED_SCHEMA, SIMULATION_SCHEMA, STATION_OBS_SCHEMA

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "build_path_from_schema",
    "structure_datasets",
]


def generate_version_hashes(in_file: Path, out_file: Path):
    hash_sha256_writer = hashlib.sha256()
    with open(in_file, "rb") as f:
        hash_sha256_writer.update(f.read())
    sha256sum = hash_sha256_writer.hexdigest()

    print(f"Writing sha256sum (ending: {sha256sum[-6:]}) to file: {out_file.name}")
    with open(out_file, "w") as f:
        f.write(sha256sum)

    del hash_sha256_writer
    del sha256sum


def _structure_datasets(
    in_file: Path, out_path: Path, method: str, dry_run: bool = False
):
    if method.lower() in ["move", "copy"]:
        meth = "Moved" if method.lower() == "move" else "Copied"
        output_file = out_path.joinpath(in_file.name)
        try:
            if not dry_run:
                method_mod = ""
                if in_file.is_dir() and method.lower() == "copy":
                    method_mod = "tree"

                if sys.version_info < (3, 9):
                    getattr(shutil, f"{method}{method_mod}")(
                        str(in_file), str(output_file)
                    )
                else:
                    getattr(shutil, f"{method}{method_mod}")(in_file, output_file)
            print(f"{meth} {in_file.name} to {output_file}.")
        except FileExistsError:
            print(f"{in_file.name} already exists at location. Continuing...")


def build_path_from_schema(
    facets: dict, output_folder: Union[str, os.PathLike]
) -> Path:
    """Build a filepath based on a valid data schema.

    Parameters
    ----------
    facets: dict
      Facets for a given dataset.
    output_folder
      Parent folder on which to extend the filetree structure.

    Returns
    -------
    Path
    """
    if facets["type"] == "station-obs":
        STATION_OBS_SCHEMA.validate(facets)
        folder_tree = (
            Path(output_folder)
            / facets["type"]
            / facets["project"]
            / facets["institution"]
            / facets["version"]  # This suggests "date_created"
            / facets["frequency"]
            / facets["variable"]
        )
        if hasattr(facets, "member"):
            return folder_tree / facets["member"]
        return folder_tree

    if facets["type"] in ["forecast", "gridded-obs", "reanalysis"]:
        GRIDDED_SCHEMA.validate(facets)
        return (
            Path(output_folder)
            / facets["type"]
            / facets["institution"]
            / facets["activity"]
            / facets["source"]
            / facets["project"]
            / facets["domain"]
            / facets["frequency"]
            / facets["variable"]
        )

    if facets["type"] == "simulation":
        SIMULATION_SCHEMA.validate(facets)
        if facets["processing_level"] == "raw":
            try:
                if facets["project"] == "CORDEX":
                    return (
                        Path(output_folder)
                        / facets["type"]
                        / facets["processing_level"]
                        / facets["activity"]
                        / facets["mip_era"]
                        / facets["project"]
                        / facets["domain"]
                        / facets["source"]
                        / facets["driving_model"]
                        / facets["experiment"]
                        / facets["member"]
                        / facets["frequency"]
                        / facets["variable"]
                    )
            except KeyError:
                return (
                    Path(output_folder)
                    / facets["type"]
                    / facets["processing_level"]
                    / facets["activity"]
                    / facets["mip_era"]
                    / facets["domain"]
                    / facets["institution"]
                    / facets["source"]
                    / facets["experiment"]
                    / facets["member"]
                    / facets["frequency"]
                    / facets["variable"]
                )
        elif (
            facets["processing_level"] == "biasadjusted"
        ):  # FIXME: remove or keep underline?
            return (
                Path(output_folder)
                / facets["type"]
                / facets["processing_level"]
                / facets["activity"]
                / facets["mip_era"]
                / facets["bias_adjust_institution"]
                / facets["bias_adjust_project"]
                / facets["domain"]
                / facets["institution"]
                / facets["source"]
                / facets["experiment"]
                / facets["member"]
                / facets["frequency"]
                / facets["variable"]
            )

    raise ValueError("No appropriate data schemas found.")


def structure_datasets(
    input_files: Union[str, os.PathLike, List[Union[str, os.PathLike]], GeneratorType],
    output_folder: Union[str, os.PathLike],
    *,
    project: Optional[str] = None,
    guess: bool = True,
    dry_run: bool = False,
    method: str = "copy",
    make_dirs: bool = False,
    set_version_hashes: bool = False,
    filename_pattern: str = "*.nc",
) -> Mapping[Path, Path]:
    """

    Parameters
    ----------
    input_files: str or Path or list of str or Path or GeneratorType
    output_folder: str or Path
    project: {"cordex", "cmip5", "cmip6", "isimip-ft", "reanalysis", "pcic-candcs-u6"}, optional
    guess: bool
      If project not supplied, suggest to decoder that project is the same for all input_files. Default: True.
    dry_run: bool
      Prints changes that would have been made without performing them. Default: False.
    method: {"move", "copy"}
      Method to transfer files to intended location. Default: "move".
    make_dirs:
      Make folder tree if it does not already exist. Default: False.
    set_version_hashes:
      Make an accompanying file with version in filename and sha256sum in contents. Default: False.
    filename_pattern: str
      If pattern ends with "zarr", will 'glob' with provided pattern.
      Otherwise, will perform an 'rglob' (recursive) operation.

    Returns
    -------
    dict
    """
    input_files = discover_data(input_files, filename_pattern)
    if not project and guess:
        # Examine the first file from a list or generator
        for f in input_files:
            project = guess_project(f)
            decoder = Decoder(project)
            decoder.decode(f)
            break
        else:
            raise FileNotFoundError()
        decoder.decode(input_files)
    else:
        decoder = Decoder(project)
        decoder.decode(input_files)

    all_file_paths = dict()
    version_hash_paths = dict()
    for file, facets in decoder.file_facets().items():
        output_filepath = build_path_from_schema(facets, output_folder)
        all_file_paths.update({Path(file): output_filepath})

        if set_version_hashes:
            version_hash_file = f"{Path(file).stem}.{facets['version']}"
            version_hash_paths.update(
                {Path(file): output_filepath.joinpath(version_hash_file)}
            )

    if make_dirs:
        for new_paths in set(all_file_paths.values()):
            Path(new_paths).mkdir(exist_ok=True, parents=True)

    if set_version_hashes:
        with multiprocessing.Pool() as pool:
            pool.starmap(
                generate_version_hashes,
                zip(version_hash_paths.keys(), version_hash_paths.values()),
            )
            pool.close()
            pool.join()

    # multiprocessing copy
    structure_func = partial(_structure_datasets, method=method, dry_run=dry_run)
    with multiprocessing.Pool() as pool:
        pool.starmap(
            structure_func, zip(all_file_paths.keys(), all_file_paths.values())
        )
        pool.close()
        pool.join()

    return all_file_paths
