import hashlib
import logging.config
import multiprocessing
import os
import shutil
import sys
from functools import partial
from pathlib import Path
from types import GeneratorType
from typing import Dict, List, Mapping, Optional, Union

import schema

from miranda.decode import Decoder, guess_project
from miranda.scripting import LOGGING_CONFIG
from miranda.utils import discover_data
from miranda.validators import GRIDDED_SCHEMA, SIMULATION_SCHEMA, STATION_OBS_SCHEMA

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "create_version_hashes",
    "build_path_from_schema",
    "structure_datasets",
]


def _generate_version_hashes(
    in_file: Path, out_file: Path, verify: bool = False
) -> None:
    if not out_file.exists():
        hash_sha256_writer = hashlib.sha256()
        with open(in_file, "rb") as f:
            hash_sha256_writer.update(f.read())
        sha256sum = hash_sha256_writer.hexdigest()

        print(f"Writing sha256sum (ending: {sha256sum[-6:]}) to file: {out_file.name}")
        try:
            with open(out_file, "w") as f:
                f.write(sha256sum)
        except PermissionError:
            logging.error("Unable to write file. Ensure access privileges.")

        del hash_sha256_writer
        del sha256sum
    elif verify:
        hash_sha256_writer = hashlib.sha256()
        with open(in_file, "rb") as f:
            hash_sha256_writer.update(f.read())
        calculated_sha256sum = hash_sha256_writer.hexdigest()

        try:
            with open(out_file) as f:
                found_sha256sum = f.read()

            if calculated_sha256sum != found_sha256sum:
                raise ValueError()
        except ValueError:
            logging.error(
                f"Found sha256sum (ending: {found_sha256sum[-6:]}) "
                f"does not match current value (ending: {calculated_sha256sum[-6:]}) "
                f"for file `{in_file.name}."
            )

    else:
        print(f"Writing sha256sum file `{out_file.name}` exists. Continuing...")


def create_version_hashes(
    input_files: Optional[Union[os.PathLike, List[os.PathLike], GeneratorType]] = None,
    facet_dict: Optional[Dict] = None,
    verify_hash: bool = False,
) -> None:
    if not facet_dict and not input_files:
        raise ValueError()

    if input_files:
        if isinstance(input_files, os.PathLike):
            input_files = [input_files]
        for f in input_files:
            project = guess_project(f)
            decoder = Decoder(project)
            decoder.decode(f)
            break
        else:
            raise FileNotFoundError()
        decoder.decode(input_files)
        facet_dict = decoder.file_facets()

    version_hash_paths = dict()
    for file, facets in facet_dict.items():
        version_hash_file = f"{Path(file).stem}.{facets['version']}"
        version_hash_paths.update(
            {Path(file): Path(file).parent.joinpath(version_hash_file)}
        )

    hash_func = partial(_generate_version_hashes, verify=verify_hash)
    with multiprocessing.Pool() as pool:
        pool.starmap(
            hash_func,
            zip(version_hash_paths.keys(), version_hash_paths.values()),
        )
        pool.close()
        pool.join()


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
) -> Optional[Path]:
    """Build a filepath based on a valid data schema.

    Parameters
    ----------
    facets: dict
      Facets for a given dataset.
    output_folder
      Parent folder on which to extend the filetree structure.

    Returns
    -------
    Path or None
    """
    folder_tree_structure = None
    try:
        if facets["type"] == "station-obs":
            STATION_OBS_SCHEMA.validate(facets)
            folder_tree_structure = (
                "type",
                "institution",
                "source",
                "version",  # This suggests "date_created"
                "frequency",
                "member"
                if hasattr(facets, "member")
                else None,  # This suggests station code
                "variable",
            )

        if facets["type"] in ["forecast", "gridded-obs", "reconstruction"]:
            GRIDDED_SCHEMA.validate(facets)
            folder_tree_structure = (
                "type",
                "institution",
                "source",
                "domain",
                "frequency",
                "variable",
            )

        if facets["type"] == "simulation":
            SIMULATION_SCHEMA.validate(facets)
            folder_tree_structure = (
                "type",
                "processing_level",
                "mip_era",
                "activity"
                if facets["processing_level"] == "raw"
                else "bias_adjust_project",
                "domain",
                "institution",
                "source",
                "driving_model" if facets["activity"] == "CORDEX" else None,
                "experiment",
                "member",
                "frequency",
                "variable",
            )

        if folder_tree_structure:
            facet_tree = list()
            for facet in folder_tree_structure:
                if facet:
                    facet_tree.append(facets[facet])
            return Path(output_folder).joinpath("/".join(facet_tree))
        else:
            raise ValueError()

    except schema.SchemaError as e:
        logging.error(
            f"Validation issues found for file matching schema: {facets}: {e}"
        )
    except ValueError:
        logging.error(
            f"No appropriate data schemas found for file matching schema: {facets}"
        )
    return


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
    verify_hashes: bool = False,
    filename_pattern: str = "*.nc",
) -> Mapping[Path, Path]:
    """

    Parameters
    ----------
    input_files: str or Path or list of str or Path or GeneratorType
    output_folder: str or Path
    project: {"cordex", "cmip5", "cmip6", "isimip-ft", "converted", "pcic-candcs-u6"}, optional
    guess: bool
      If project not supplied, suggest to decoder that activity is the same for all input_files. Default: True.
    dry_run: bool
      Prints changes that would have been made without performing them. Default: False.
    method: {"move", "copy"}
      Method to transfer files to intended location. Default: "move".
    make_dirs: bool
      Make folder tree if it does not already exist. Default: False.
    set_version_hashes: bool
      Make an accompanying file with version in filename and sha256sum in contents. Default: False.
    verify_hashes: bool
      Ensure that any existing she256sum files correspond with companion file. Raise on error. Default: False.
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
    existing_hashes = dict()
    version_hash_paths = dict()
    errored_files = list()
    for file, facets in decoder.file_facets().items():
        output_filepath = build_path_from_schema(facets, output_folder)
        if isinstance(output_filepath, Path):
            all_file_paths.update({Path(file): output_filepath})
        else:
            errored_files.append(Path(file).name)
            continue

        if set_version_hashes:
            version_hash_file = f"{Path(file).stem}.{facets['version']}"
            if Path(file).parent.joinpath(version_hash_file).exists():
                existing_hashes.update(
                    {
                        Path(file).parent.joinpath(
                            version_hash_file
                        ): output_filepath.joinpath(version_hash_file)
                    }
                )
            else:
                version_hash_paths.update(
                    {Path(file): output_filepath.joinpath(version_hash_file)}
                )

    if errored_files:
        logging.warning(
            f"Some files were unable to be structured: [{', '.join(errored_files)}]"
        )

    if make_dirs:
        for new_paths in set(all_file_paths.values()):
            Path(new_paths).mkdir(exist_ok=True, parents=True)

    if set_version_hashes:
        hash_func = partial(_generate_version_hashes, verify=verify_hashes)
        with multiprocessing.Pool() as pool:
            if existing_hashes:
                pool.starmap(
                    getattr(shutil, method),
                    zip(existing_hashes.keys(), existing_hashes.values()),
                )
            pool.starmap(
                hash_func,
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
