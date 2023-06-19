from __future__ import annotations

import hashlib
import logging.config
import multiprocessing
import os
import shutil
from functools import partial
from pathlib import Path
from types import GeneratorType

import yaml
from schema import SchemaError

from miranda.cv import VALIDATION_ENABLED
from miranda.decode import Decoder, DecoderError, guess_project
from miranda.io import discover_data
from miranda.scripting import LOGGING_CONFIG

if VALIDATION_ENABLED:
    from miranda.validators import validation_schemas

logging.config.dictConfig(LOGGING_CONFIG)

__all__ = [
    "create_version_hash_files",
    "build_path_from_schema",
    "structure_datasets",
]


def _verify(hash_value: str, hash_file: os.PathLike) -> None:
    try:
        with open(hash_file) as f:
            found_sha256sum = f.read()

        if hash_value != found_sha256sum:
            raise ValueError()
    except ValueError:
        logging.error(
            f"Found sha256sum (starting: {found_sha256sum[:6]}) "
            f"does not match current value (starting: {hash_value[:6]}) "
            f"for file `{Path(hash_file).name}."
        )


def generate_hash_file(
    in_file: os.PathLike, out_file: os.PathLike, verify: bool = False
) -> None:
    if not Path(out_file).exists():
        hash_sha256_writer = hashlib.sha256()
        with open(in_file, "rb") as f:
            hash_sha256_writer.update(f.read())
        sha256sum = hash_sha256_writer.hexdigest()

        print(
            f"Writing sha256sum (starting: {sha256sum[:6]}) to file: {Path(out_file).name}"
        )
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

        _verify(calculated_sha256sum, out_file)

    else:
        print(f"Writing sha256sum file `{Path(out_file).name}` exists. Continuing...")


def generate_hash_metadata(
    in_file: os.PathLike,
    version: str | None = None,
    hash_file: os.PathLike | None = None,
    verify: bool = False,
) -> dict[str, list[str]]:
    hashversion = dict()

    if version is None:
        version = "vNotFound"

    if not Path(hash_file).exists():
        hash_sha256_writer = hashlib.sha256()
        with open(in_file, "rb") as f:
            hash_sha256_writer.update(f.read())
        sha256sum = hash_sha256_writer.hexdigest()

        print(f"Calculated sha256sum (starting: {sha256sum[:6]})")

        hashversion[Path(in_file).name] = [version, sha256sum]
        del hash_sha256_writer

    else:
        hash_sha256_writer = hashlib.sha256()
        with open(in_file, "rb") as f:
            hash_sha256_writer.update(f.read())
        calculated_sha256sum = hash_sha256_writer.hexdigest()

        if verify:
            _verify(calculated_sha256sum, hash_file)

        hashversion[Path(in_file).name] = [version, calculated_sha256sum]

    return hashversion


def create_version_hash_files(
    input_files: str
    | os.PathLike
    | list[str | os.PathLike]
    | GeneratorType
    | None = None,
    facet_dict: dict | None = None,
    verify_hash: bool = False,
) -> None:
    """Create version hashes based on files or a facets dictionary.

    Parameters
    ----------
    input_files : str, os.PathLike, list of str or os.PathLike, or GeneratorType
    facet_dict : dict, optional
    verify_hash : bool

    Returns
    -------
    None
    """
    if not facet_dict and not input_files:
        raise ValueError("Facets dictionary or sequence of filepaths required.")

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

    hash_func = partial(generate_hash_file, verify=verify_hash)
    with multiprocessing.Pool() as pool:
        pool.starmap(
            hash_func,
            zip(version_hash_paths.keys(), version_hash_paths.values()),
        )
        pool.close()
        pool.join()


def _structure_datasets(
    in_file: str | os.PathLike,
    out_path: str | os.PathLike,
    method: str,
    dry_run: bool = False,
):
    if isinstance(in_file, str):
        in_file = Path(in_file)

    if method.lower() in ["move", "copy"]:
        meth = "Moved" if method.lower() == "move" else "Copied"
        output_file = Path(out_path).joinpath(in_file.name)
        try:
            if not dry_run:
                method_mod = ""
                if in_file.is_dir() and method.lower() == "copy":
                    method_mod = "tree"

                getattr(shutil, f"{method}{method_mod}")(in_file, output_file)
            print(f"{meth} {in_file.name} to {output_file}.")
        except FileExistsError:
            print(f"{in_file.name} already exists at location. Continuing...")


def parse_schema(
    facets: dict, schema: str | os.PathLike | dict, top_folder: str = "datasets"
) -> list:
    """Parse the schema from a YAML schema configuration and construct path using a dictionary of facets.

    Parameters
    ----------
    facets : dict
    schema : str or os.PathLike or dict
    top_folder : str

    Returns
    -------
    list
    """

    def _parse_top_level(schematic: dict, facet_dict: dict, top: str):
        try:
            parent = schematic[top]
        except KeyError:
            logging.error("Schema is not a valid facet-tree reference.")
            raise

        for i, options in enumerate(parent):
            if {"option", "structure", "value"}.issubset(options.keys()):
                option = options["option"]
                value = options["value"]

                if option in facet_dict.keys():
                    if facet_dict[option] == value:
                        return {"branch": value, "structure": options["structure"]}
                    continue
        raise ValueError("Couldn\nt parse top level.")

    def _parse_structure(branch_dict: dict, facet_dict: dict) -> list:
        structure = branch_dict.get("structure")
        folder_tree = list()

        for level in structure:
            if isinstance(level, str):
                folder_tree.append(level)
                continue
            elif isinstance(level, dict):
                if {"option", "is_true"}.issubset(level.keys()):
                    option = level["option"]

                    if option not in facet_dict and "value" in level:
                        raise ValueError(
                            f"Necessary facet not found for schema: `{option}`."
                        )

                    is_true = level.get("is_true")
                    else_value = level.get("else")
                    facet = facet_dict.get(option)

                    if "value" not in level:
                        # The absence of "value" means that "is_true / else" refer to the presence or not of "option" in the facets
                        # We also treat falsy values (empty string, None) as the absence of "option" from the facets
                        if not bool(facet) and else_value:
                            folder_tree.append(else_value)
                        elif bool(facet):
                            folder_tree.append(is_true)
                        else:
                            # "option" absent from the facets and no "else": skip.
                            pass
                    else:
                        value = level["value"]
                        if facet_dict[option] == value:
                            folder_tree.append(is_true)
                        elif else_value:
                            folder_tree.append(else_value)
                        else:
                            # "option" not equal to "value", but no "else" : skip.
                            pass
            else:
                raise ValueError("Supplied schema is invalid.")
        return folder_tree

    if isinstance(schema, (str, os.PathLike)):
        with Path(schema).open() as f:
            schema = yaml.safe_load(f.read())

    branch = _parse_top_level(schema, facets, top_folder)
    tree = list()  # noqa
    tree.extend(_parse_structure(branch, facets))

    return tree


def build_path_from_schema(
    facets: dict,
    output_folder: str | os.PathLike,
    schema: str | os.PathLike | dict | None = None,
    top_folder: str = "datasets",
    validate: bool = True,
) -> Path | None:
    """Build a filepath based on a valid data schema.

    Parameters
    ----------
    facets : dict
        Facets for a given dataset.
    output_folder : str or os.PathLike
        Parent folder on which to extend the filetree structure.
    schema : str or os.PathLike, optional
        Path to YAML schematic of database structure. If None, will use Ouranos schema.
    top_folder : str
        Top-level of supplied schema, used for validation purposes. Default: "datasets".
    validate: bool
        Run facets-validation checks over given file. Default: True.

    Returns
    -------
    Path or None
    """
    if schema is None:
        schema = Path(__file__).parent.joinpath("data").joinpath("ouranos_schema.yml")

    tree = parse_schema(facets, schema, top_folder)
    branch = tree[0]

    if validate and VALIDATION_ENABLED:
        if facets[branch] in validation_schemas.keys():
            try:
                validation_schemas[facets[branch]].validate(facets)
            except SchemaError as e:
                logging.error(
                    f"Validation issues found for file matching schema: {facets}: {e}"
                )
                return
        elif facets[branch] not in validation_schemas.keys():
            logging.error(
                f"No appropriate data schemas found for file matching schema: {facets}",
                DecoderError,
            )
            return
    elif validate and not VALIDATION_ENABLED:
        logging.warning(
            "Facets validation requires pyessv-archive source files. Skipping validation checks."
        )

    file_location = list()
    for facet in tree:
        # Remove spaces in folder paths
        file_location.append(str(facets[facet]).replace(" ", "-"))
    return Path(output_folder).joinpath("/".join(file_location))


def structure_datasets(
    input_files: str | os.PathLike | list[str | os.PathLike] | GeneratorType,
    output_folder: str | os.PathLike,
    *,
    project: str | None = None,
    guess: bool = True,
    dry_run: bool = False,
    method: str = "copy",
    make_dirs: bool = False,
    set_version_hashes: bool = False,
    verify_hashes: bool = False,
    suffix: str = "nc",
) -> dict[Path, Path]:
    """

    Parameters
    ----------
    input_files : str, Path, list of str or Path, or GeneratorType
        Files to be sorted.
    output_folder : str or Path
        The desired location for the folder-tree.
    project : {"cordex", "cmip5", "cmip6", "isimip-ft", "pcic-candcs-u6", "converted"}, optional
        Project used to parse the facets of all supplied datasets.
        If not supplied, will attempt parsing with all available data categories for each file (slow)
        unless `guess` is True.
    guess : bool
        If project not supplied, suggest to decoder that activity is the same for all input_files. Default: True.
    dry_run : bool
        Prints changes that would have been made without performing them. Default: False.
    method : {"move", "copy"}
        Method to transfer files to intended location. Default: "move".
    make_dirs : bool
        Make folder tree if it does not already exist. Default: False.
    set_version_hashes : bool
        Make an accompanying file with version in filename and sha256sum in contents. Default: False.
    verify_hashes : bool
        Ensure that any existing she256sum files correspond with companion file. Raise on error. Default: False.
    suffix : {"nc", "zarr"}
        If "zarr", will perform a 'glob' with provided pattern.
        Otherwise, will perform an 'rglob' (recursive) operation.

    Returns
    -------
    dict[Path, Path]
    """
    input_files = discover_data(input_files, suffix)
    if guess and project is None:
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
        if len(errored_files) < 10:
            logging.warning(
                f"Some files were unable to be structured: [{', '.join(errored_files)}]"
            )
        else:
            logging.warning(
                f"Many files were unable to be structured (n={len(errored_files)})"
            )

    if make_dirs:
        for new_paths in set(all_file_paths.values()):
            Path(new_paths).mkdir(exist_ok=True, parents=True)

    if set_version_hashes:
        hash_func = partial(generate_hash_file, verify=verify_hashes)
        with multiprocessing.Pool() as pool:
            if existing_hashes:
                print(
                    f"Sha256sum signatures exist for {len(existing_hashes)} files. "
                    f"Transferring them via `{method}` method."
                )
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
