"""Testing utilities module."""

from __future__ import annotations
import importlib.metadata as ilm
import importlib.resources as ilr
import logging
import os
import re
import time
import warnings
from collections.abc import Callable
from datetime import datetime as dt
from functools import wraps
from io import StringIO
from pathlib import Path
from shutil import copytree
from typing import IO, Any, TextIO
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlparse
from urllib.request import urlretrieve

from filelock import FileLock
from packaging.version import Version
from xarray import Dataset
from xarray import open_dataset as _open_dataset
from xclim.testing.utils import show_versions as _show_versions

import miranda


try:
    import pooch
except ImportError:
    warnings.warn("The `pooch` library is not installed. The default cache directory for testing data will not be set.", stacklevel=2)
    pooch = None

logger = logging.getLogger("miranda")


__all__ = [
    "TESTDATA_BRANCH",
    "TESTDATA_CACHE_DIR",
    "TESTDATA_REPO_URL",
    "audit_url",
    "cassini",
    "default_testdata_cache",
    "default_testdata_repo_url",
    "default_testdata_version",
    "gather_testing_data",
    "open_dataset",
    "populate_testing_data",
    "publish_release_notes",
    "show_versions",
    "testing_setup_warnings",
]

default_testdata_version = "v2025.5.16"
"""Default version of the testing data to use when fetching datasets."""

default_testdata_repo_url = "https://raw.githubusercontent.com/Ouranosinc/miranda-testdata/"
"""Default URL of the testing data repository to use when fetching datasets."""

try:
    default_testdata_cache = Path(pooch.os_cache("miranda-testdata"))
    """Default location for the testing data cache."""
except (AttributeError, TypeError):
    default_testdata_cache = None

TESTDATA_REPO_URL = str(os.getenv("MIRANDA_TESTDATA_REPO_URL", default_testdata_repo_url))
"""
Sets the URL of the testing data repository to use when fetching datasets.

Notes
-----
When running tests locally, this can be set for both `pytest` and `tox` by exporting the variable:

.. code-block:: console

    $ export MIRANDA_TESTDATA_REPO_URL="https://github.com/my_username/miranda-testdata"

or setting the variable at runtime:

.. code-block:: console

    $ env MIRANDA_TESTDATA_REPO_URL="https://github.com/my_username/miranda-testdata" pytest
"""

TESTDATA_BRANCH = str(os.getenv("MIRANDA_TESTDATA_BRANCH", default_testdata_version))
"""
Sets the branch of the testing data repository to use when fetching datasets.

Notes
-----
When running tests locally, this can be set for both `pytest` and `tox` by exporting the variable:

.. code-block:: console

    $ export MIRANDA_TESTDATA_BRANCH="my_testing_branch"

or setting the variable at runtime:

.. code-block:: console

    $ env MIRANDA_TESTDATA_BRANCH="my_testing_branch" pytest
"""

TESTDATA_CACHE_DIR = os.getenv("MIRANDA_TESTDATA_CACHE_DIR", default_testdata_cache)
"""
Sets the directory to store the testing datasets.

If not set, the default location will be used (based on ``platformdirs``, see :func:`pooch.os_cache`).

Notes
-----
When running tests locally, this can be set for both `pytest` and `tox` by exporting the variable:

.. code-block:: console

    $ export MIRANDA_TESTDATA_CACHE_DIR="/path/to/my/data"

or setting the variable at runtime:

.. code-block:: console

    $ env MIRANDA_TESTDATA_CACHE_DIR="/path/to/my/data" pytest
"""


# Publishing Tools ###


def publish_release_notes(
    style: str = "md",
    file: os.PathLike[str] | StringIO | TextIO | None = None,
    changes: str | os.PathLike[str] | None = None,
) -> str | None:
    """
    Format release notes in Markdown or ReStructuredText.

    Parameters
    ----------
    style : {"rst", "md"}
        Use ReStructuredText formatting or Markdown. Default: Markdown.
    file : {os.PathLike, StringIO, TextIO}, optional
        If provided, prints to the given file-like object. Otherwise, returns a string.
    changes : str or os.PathLike[str], optional
        If provided, manually points to the file where the changelog can be found.
        Assumes a relative path otherwise.

    Returns
    -------
    str, optional
        If `file` not provided, the formatted release notes.

    Notes
    -----
    This function is used solely for development and packaging purposes.
    """
    if isinstance(changes, str | Path):
        changes_file = Path(changes).absolute()
    else:
        changes_file = Path(__file__).absolute().parents[3].joinpath("CHANGELOG.rst")

    if not changes_file.exists():
        raise FileNotFoundError("Changelog file not found in miranda folder tree.")

    with Path(changes_file).open(encoding="utf-8") as hf:
        changes = hf.read()

    if style == "rst":
        hyperlink_replacements = {
            r":issue:`([0-9]+)`": r"`GH/\1 <https://github.com/Ouranosinc/miranda/issues/\1>`_",
            r":pull:`([0-9]+)`": r"`PR/\1 <https://github.com/Ouranosinc/miranda/pull/\>`_",
            r":user:`([a-zA-Z0-9_.-]+)`": r"`@\1 <https://github.com/\1>`_",
        }
    elif style == "md":
        hyperlink_replacements = {
            r":issue:`([0-9]+)`": r"[GH/\1](https://github.com/Ouranosinc/miranda/issues/\1)",
            r":pull:`([0-9]+)`": r"[PR/\1](https://github.com/Ouranosinc/miranda/pull/\1)",
            r":user:`([a-zA-Z0-9_.-]+)`": r"[@\1](https://github.com/\1)",
        }
    else:
        msg = f"Formatting style not supported: {style}"
        raise NotImplementedError(msg)

    for search, replacement in hyperlink_replacements.items():
        changes = re.sub(search, replacement, changes)

    if style == "md":
        changes = changes.replace("=========\nChangelog\n=========", "# Changelog")

        titles = {r"\n(.*?)\n([\-]{1,})": "-", r"\n(.*?)\n([\^]{1,})": "^"}
        for title_expression, level in titles.items():
            found = re.findall(title_expression, changes)
            for grouping in found:
                fixed_grouping = str(grouping[0]).replace("(", r"\(").replace(")", r"\)")
                search = rf"({fixed_grouping})\n([\{level}]{'{' + str(len(grouping[1])) + '}'})"
                replacement = f"{'##' if level == '-' else '###'} {grouping[0]}"
                changes = re.sub(search, replacement, changes)

        link_expressions = r"[\`]{1}([\w\s]+)\s<(.+)>`\_"
        found = re.findall(link_expressions, changes)
        for grouping in found:
            search = rf"`{grouping[0]} <.+>`\_"
            replacement = f"[{str(grouping[0]).strip()}]({grouping[1]})"
            changes = re.sub(search, replacement, changes)

    if not file:
        return changes
    if isinstance(file, Path | os.PathLike):
        with Path(file).open(mode="w", encoding="utf-8") as f:
            print(changes, file=f)
    else:
        print(changes, file=file)
    return None


def show_versions(
    file: os.PathLike | StringIO | TextIO | None = None,
    deps: list | None = None,
) -> str | None:
    """
    Print the versions of miranda and its dependencies.

    Parameters
    ----------
    file : {os.PathLike, StringIO, TextIO}, optional
        If provided, prints to the given file-like object. Otherwise, returns a string.
    deps : list, optional
        A list of dependencies to gather and print version information from. Otherwise, prints `miranda` dependencies.

    Returns
    -------
    str or None
        The formatted version information if `file` is not provided, otherwise None.
    """

    def _get_miranda_dependencies():
        miranda_metadata = ilm.metadata("miranda")
        requires = miranda_metadata.get_all("Requires-Dist")
        requires = [req.split("[")[0].split(";")[0].split(">")[0].split("<")[0].split("=")[0].split("!")[0].strip() for req in requires]
        sorted_deps = sorted(list(set(requires) - {"miranda"}))

        return ["miranda"] + sorted_deps

    if deps is None:
        deps = _get_miranda_dependencies()

    return _show_versions(file=file, deps=deps)


# Test Data Utilities ###


def testing_setup_warnings():
    """Warn users about potential incompatibilities between miranda and miranda-testdata versions."""
    if re.match(r"^\d+\.\d+\.\d+$", miranda.__version__) and TESTDATA_BRANCH != default_testdata_version:
        # This does not need to be emitted on GitHub Workflows and ReadTheDocs
        if not os.getenv("CI") and not os.getenv("READTHEDOCS"):
            warnings.warn(
                f"`miranda` stable ({miranda.__version__}) is running tests against a non-default "
                f"branch of the testing data. It is possible that changes to the testing data may "
                f"be incompatible with some assertions in this version. "
                f"Please be sure to check {TESTDATA_REPO_URL} for more information.",
                stacklevel=2,
            )

    if re.match(r"^v\d+\.\d+\.\d+", TESTDATA_BRANCH):
        # Find the date of last modification of miranda source files to generate a calendar version
        install_date = dt.strptime(
            time.ctime(Path(miranda.__file__).stat().st_mtime),
            "%a %b %d %H:%M:%S %Y",
        )
        install_calendar_version = f"{install_date.year}.{install_date.month}.{install_date.day}"

        if Version(TESTDATA_BRANCH) > Version(install_calendar_version):
            warnings.warn(
                f"The installation date of `miranda` ({install_date.ctime()}) "
                f"predates the last release of testing data ({TESTDATA_BRANCH}). "
                "It is very likely that the testing data is incompatible with this build of `miranda`.",
                stacklevel=2,
            )


def load_registry(branch: str = TESTDATA_BRANCH, repo: str = TESTDATA_REPO_URL) -> dict[str, str]:
    """
    Load the registry file for the test data.

    Parameters
    ----------
    branch : str
        Branch of the repository to use when fetching testing datasets.
    repo : str
        URL of the repository to use when fetching testing datasets.

    Returns
    -------
    dict
        Dictionary of filenames and hashes.
    """
    if not repo.endswith("/"):
        repo = f"{repo}/"
    remote_registry = audit_url(
        urljoin(
            urljoin(repo, branch if branch.endswith("/") else f"{branch}/"),
            "data/registry.txt",
        )
    )

    if repo != default_testdata_repo_url:
        external_repo_name = urlparse(repo).path.split("/")[-2]
        external_branch_name = branch.split("/")[-1]
        registry_file = Path(str(ilr.files("miranda").joinpath(f"testing/registry.{external_repo_name}.{external_branch_name}.txt")))
        urlretrieve(remote_registry, registry_file)  # noqa: S310

    elif branch != default_testdata_version:
        custom_registry_folder = Path(str(ilr.files("miranda").joinpath(f"testing/{branch}")))
        custom_registry_folder.mkdir(parents=True, exist_ok=True)
        registry_file = custom_registry_folder.joinpath("registry.txt")
        urlretrieve(remote_registry, registry_file)  # noqa: S310

    else:
        registry_file = Path(str(ilr.files("miranda").joinpath("testing/registry.txt")))

    if not registry_file.exists():
        raise FileNotFoundError(f"Registry file not found: {registry_file}")

    # Load the registry file
    with registry_file.open(encoding="utf-8") as f:
        registry = {line.split()[0]: line.split()[1] for line in f}
    return registry


def cassini(
    repo: str = TESTDATA_REPO_URL,
    branch: str = TESTDATA_BRANCH,
    cache_dir: str | Path = TESTDATA_CACHE_DIR,
    allow_updates: bool = True,
):
    """
    Pooch registry instance for miranda test data.

    Parameters
    ----------
    repo : str
        URL of the repository to use when fetching testing datasets.
    branch : str
        Branch of repository to use when fetching testing datasets.
    cache_dir : str or Path
        The path to the directory where the data files are stored.
    allow_updates : bool
        If True, allow updates to the data files. Default is True.

    Returns
    -------
    pooch.Pooch
        The Pooch instance for accessing the miranda testing data.

    Notes
    -----
    There are three environment variables that can be used to control the behaviour of this registry:
        - ``MIRANDA_TESTDATA_CACHE_DIR``: If this environment variable is set, it will be used as the
          base directory to store the data files.
          The directory should be an absolute path (i.e., it should start with ``/``).
          Otherwise, the default location will be used (based on ``platformdirs``, see :py:func:`pooch.os_cache`).
        - ``MIRANDA_TESTDATA_REPO_URL``: If this environment variable is set, it will be used as the URL of
          the repository to use when fetching datasets. Otherwise, the default repository will be used.
        - ``MIRANDA_TESTDATA_BRANCH``: If this environment variable is set, it will be used as the branch of
          the repository to use when fetching datasets. Otherwise, the default branch will be used.

    Examples
    --------
    Using the registry to download a file:

    .. code-block:: python

        import xarray as xr
        from miranda.testing import cassini

        example_file = cassini().fetch("example.nc")
        data = xr.open_dataset(example_file)
    """
    if pooch is None:
        raise ImportError(
            "The `pooch` package is required to fetch the miranda testing data. "
            "You can install it with `pip install pooch` or `pip install miranda[dev]`."
        )
    if not repo.endswith("/"):
        repo = f"{repo}/"
    remote = audit_url(urljoin(urljoin(repo, branch if branch.endswith("/") else f"{branch}/"), "data"))

    _cassini = pooch.create(
        path=cache_dir,
        base_url=remote,
        version=default_testdata_version,
        version_dev=branch,
        allow_updates=allow_updates,
        registry=load_registry(branch=branch, repo=repo),
    )

    # Add a custom fetch method to the Pooch instance
    # Needed to address: https://github.com/readthedocs/readthedocs.org/issues/11763
    # Fix inspired by @bjlittle (https://github.com/bjlittle/geovista/pull/1202)
    _cassini.fetch_diversion = _cassini.fetch

    # Overload the fetch method to add user-agent headers
    @wraps(_cassini.fetch_diversion)
    def _fetch(*args, **kwargs: bool | Callable) -> str:  # numpydoc ignore=GL08  # *args: str
        def _downloader(
            url: str,
            output_file: str | IO,
            poocher: pooch.Pooch,
            check_only: bool | None = False,
        ) -> None:
            """Download the file from the URL and save it to the save_path."""
            headers = {"User-Agent": f"miranda ({miranda.__version__})"}
            downloader = pooch.HTTPDownloader(headers=headers)
            return downloader(url, output_file, poocher, check_only=check_only)

        # default to our http/s downloader with user-agent headers
        kwargs.setdefault("downloader", _downloader)
        return _cassini.fetch_diversion(*args, **kwargs)

    # Replace the fetch method with the custom fetch method
    _cassini.fetch = _fetch

    return _cassini


def open_dataset(
    name: str,
    cassini_kwargs: dict[str, Path | str | bool] | None = None,
    **xr_kwargs: Any,
) -> Dataset:
    r"""
    Convenience function to open a dataset from the miranda testing data using the `cassini` class.

    This is a thin wrapper around the `cassini` class to make it easier to open miranda testing datasets.

    Parameters
    ----------
    name : str
        Name of the file containing the dataset.
    cassini_kwargs : dict
        Keyword arguments passed to the cassini function.
    **xr_kwargs : Any
        Keyword arguments passed to xarray.open_dataset.

    Returns
    -------
    xarray.Dataset
        The dataset.

    See Also
    --------
    xarray.open_dataset : Open and read a dataset from a file or file-like object.
    cassini : Pooch wrapper for accessing the miranda testing data.
    """
    if cassini_kwargs is None:
        cassini_kwargs = {}
    return _open_dataset(cassini(**cassini_kwargs).fetch(name), **xr_kwargs)


def populate_testing_data(
    temp_folder: Path | None = None,
    repo: str = TESTDATA_REPO_URL,
    branch: str = TESTDATA_BRANCH,
    local_cache: Path = TESTDATA_CACHE_DIR,
) -> None:
    """
    Populate the local cache with the testing data.

    Parameters
    ----------
    temp_folder : Path, optional
        Path to a temporary folder to use as the local cache. If not provided, the default location will be used.
    repo : str, optional
        URL of the repository to use when fetching testing datasets.
    branch : str, optional
        Branch of miranda-testdata to use when fetching testing datasets.
    local_cache : Path
        The path to the local cache. Defaults to the location set by the platformdirs library.
        The testing data will be downloaded to this local cache.
    """
    # Create the Pooch instance
    n = cassini(repo=repo, branch=branch, cache_dir=temp_folder or local_cache)

    # Download the files
    errored_files = []
    for file in load_registry():
        try:
            n.fetch(file)
        except HTTPError:  # noqa: PERF203
            msg = f"File `{file}` not accessible in remote repository."
            logger.error(msg)
            errored_files.append(file)
        else:
            logger.info("Files were downloaded successfully.")

    if errored_files:
        logger.error(
            "The following files were unable to be downloaded: %s",
            errored_files,
        )


def gather_testing_data(
    worker_cache_dir: str | os.PathLike[str] | Path,
    worker_id: str,
    _cache_dir: str | os.PathLike[str] | None = TESTDATA_CACHE_DIR,
) -> None:
    """
    Gather testing data across workers.

    Parameters
    ----------
    worker_cache_dir : str or Path
        The directory to store the testing data.
    worker_id : str
        The worker ID.
    _cache_dir : str or Path, optional
        The directory to store the testing data. Default is None.

    Raises
    ------
    ValueError
        If the cache directory is not set.
    FileNotFoundError
        If the testing data is not found.
    """
    if _cache_dir is None:
        raise ValueError("The cache directory must be set. Please set the `cache_dir` parameter or the `MIRANDA_DATA_DIR` environment variable.")
    cache_dir = Path(_cache_dir)

    if worker_id == "master":
        populate_testing_data(branch=TESTDATA_BRANCH)
    else:
        cache_dir.mkdir(exist_ok=True, parents=True)
        lockfile = cache_dir.joinpath(".lock")
        test_data_being_written = FileLock(lockfile)
        with test_data_being_written:
            # This flag prevents multiple calls from re-attempting to download testing data in the same pytest run
            populate_testing_data(branch=TESTDATA_BRANCH)
            cache_dir.joinpath(".data_written").touch()
        with test_data_being_written.acquire():
            if lockfile.exists():
                lockfile.unlink()
        copytree(cache_dir.joinpath(default_testdata_version), worker_cache_dir)


# Testing Utilities ###


def audit_url(url: str, context: str | None = None) -> str:
    """
    Check if the URL is well-formed.

    Parameters
    ----------
    url : str
        The URL to check.
    context : str, optional
        Additional context to include in the error message. Default is None.

    Returns
    -------
    str
        The URL if it is well-formed.

    Raises
    ------
    URLError
        If the URL is not well-formed.
    """
    msg = ""
    result = urlparse(url)
    if result.scheme == "http":
        msg = f"{context if context else ''} URL is not using secure HTTP: '{url}'".strip()
    if not all([result.scheme, result.netloc]):
        msg = f"{context if context else ''} URL is not well-formed: '{url}'".strip()

    if msg:
        logger.error(msg)
        raise URLError(msg)
    return url
