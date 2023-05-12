import hashlib
import logging.config
import multiprocessing
import operator as op
import os
import shutil
import sys
from functools import partial, reduce
from pathlib import Path
from types import GeneratorType
from typing import Dict, List, Optional, Tuple, Union

import yaml
from pandas import isna
from schema import SchemaError

from miranda.convert.utils import date_parser
from miranda.cv import VALIDATION_ENABLED
from miranda.decode import Decoder, DecoderError, freq_to_timedelta, guess_project
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
    version: Optional[str] = None,
    hash_file: Optional[os.PathLike] = None,
    verify: bool = False,
) -> Dict[str, List[str]]:
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
    input_files: Optional[
        Union[str, os.PathLike, List[Union[str, os.PathLike]], GeneratorType]
    ] = None,
    facet_dict: Optional[Dict] = None,
    verify_hash: bool = False,
) -> None:
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
    in_file: Union[str, os.PathLike],
    out_file: Union[str, os.PathLike],
    method: str,
    dry_run: bool = False,
):
    if isinstance(in_file, str):
        in_file = Path(in_file)

    if method.lower() in ["move", "copy"]:
        meth = "Moved" if method.lower() == "move" else "Copied"
        output_file = Path(out_file)
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


def _parse_option(option: dict, facets: dict):
    """Parse an option element of the facet schema tree."""
    facet_value = facets[option["option"]]
    if "value" in option:
        if isinstance(option["value"], str):
            answer = facet_value == option["value"]
        else:  # A list
            answer = facet_value in option["value"]
    else:
        answer = not isna(facet_value)

    if "is_true" in option and answer:
        return option["is_true"]
    if "else" in option and not answer:
        return option["else"]
    return answer


def _parse_level(schema: Union[dict, str], facets: dict):
    if isinstance(schema, str):
        # A single facet:
        if isna(facets[schema]):
            raise ValueError(f"Schema requires a value for facet {schema}.")
        return facets[schema]
    if "option" in schema:
        answer = _parse_option(schema, facets)
        if isinstance(answer, bool) and not answer:
            # Test failed with no "else" value, we skip this level.
            return None
        return _parse_level(answer, facets)
    if "concat" in schema:
        parts = []
        for element in schema["concat"]:
            part = _parse_level(element, facets)
            if not isna(part):
                parts.append(part)
        return "_".join(parts)
    if "text" in schema:
        return schema["text"]
    raise ValueError(f"Invalid schema : {schema}")


def _parse_dates(facets):
    if facets["frequency"] == "fx":
        return "fx"

    start = date_parser(facets["date_start"], output_type="datetime")
    end = date_parser(facets["date_end"], output_type="datetime")
    freq = freq_to_timedelta(facets["frequency"])

    # Full years : Starts on Jan 1st and is either annual or ends on Dec 31st (accepting Dec 30 for 360 cals)
    if (
        start.month == 1
        and start.day == 1
        and (freq >= freq_to_timedelta("yr") or (end.month == 12 and end.day > 29))
    ):
        if start.year == end.year:
            return f"{start:%Y}"
        return f"{start:%Y}-{end:%Y}"
    # Full months : Starts on the 1st and is either montly or ends on the last day
    if start.day == 1 and (freq >= freq_to_timedelta("mon") or end.day > 27):
        # Full months
        if (start.year, start.month) == (end.year, end.month):
            return f"{start:%Y%m}"
        return f"{start:%Y%m}-{end:%Y%m}"
    # The full range
    return f"{start:%Y%m%d}-{end:%Y%m%d}"


def _parse_filename(schema: list, facets: dict) -> str:
    return "_".join(
        [
            facets[element] if element != "dates" else _parse_dates(facets)
            for element in schema
            if element == "dates" or not isna(facets[element])
        ]
    )


def _parse_structure(schema: list, facets: dict) -> list:
    folder_tree = list()
    for level in schema:
        part = _parse_level(level, facets)
        if not isna(part):
            folder_tree.append(part)
    return folder_tree


def parse_schema(
    facets: dict, schema: Union[str, os.PathLike, dict], top_folder: str = "datasets"
) -> Tuple[List[str], str]:
    """Parse the schema from a YAML schema configuration and construct path using a dictionary of facets.

    Parameters
    ----------
    facets : dict
    schema : str or os.PathLike or dict
    top_folder : str

    Returns
    -------
    list of folders, filename without suffix
    """
    if isinstance(schema, (str, os.PathLike)):
        with Path(schema).open() as f:
            schema = yaml.safe_load(f.read())

    try:
        parent = schema[top_folder]
    except KeyError:
        logging.error("Schema is not a valid facet-tree reference.")
        raise

    for i, structure in enumerate(parent):
        if {"with", "structure", "filename"} != set(structure.keys()):
            raise ValueError("Invalid schema specification.")

        match = reduce(
            op.and_, map(partial(_parse_option, facets=facets), structure["with"])
        )
        if match:
            return _parse_structure(structure["structure"], facets), _parse_filename(
                structure["filename"], facets
            )
    raise ValueError(
        f"This file doesn't match any registered structure. Facets:\n{facets}"
    )


def build_path_from_schema(
    facets: dict,
    output_folder: Union[str, os.PathLike],
    schema: Optional[Union[str, os.PathLike, dict]] = None,
    top_folder: str = "datasets",
    validate: bool = True,
) -> Optional[Path]:
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
      The Path includes the filename, without suffix.
    """
    if validate and VALIDATION_ENABLED:
        if facets["type"] in validation_schemas.keys():
            try:
                validation_schemas[facets["type"]].validate(facets)
            except SchemaError as e:
                logging.error(
                    f"Validation issues found for file matching schema: {facets}: {e}"
                )
                return
        elif facets["type"] not in validation_schemas.keys():
            logging.error(
                f"No appropriate data schemas found for file matching schema: {facets}",
                DecoderError,
            )
            return
    elif validate and not VALIDATION_ENABLED:
        logging.warning(
            "Facets validation requires pyessv-archive source files. Skipping validation checks."
        )

    if schema is None:
        schema = Path(__file__).parent.joinpath("data").joinpath("ouranos_schema.yml")

    tree, filename = parse_schema(facets, schema, top_folder)
    return Path(output_folder).joinpath("/".join(tree)) / filename


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
    suffix: str = "nc",
) -> Dict[Path, Path]:
    """

    Parameters
    ----------
    input_files : str or Path or list of str or Path or GeneratorType
    output_folder : str or Path

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
    Dict[Path, Path]
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
            all_file_paths.update(
                {Path(file): output_filepath.with_suffix(Path(file).suffix)}
            )
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
