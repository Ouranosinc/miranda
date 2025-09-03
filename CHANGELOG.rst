.. :changelog:

=========
Changelog
=========

..
    `Unreleased <https://github.com/Ouranosinc/miranda>`_ (latest)
    --------------------------------------------------------------

    Contributors:

    Changes
    ^^^^^^^
    * No change.

    Fixes
    ^^^^^
    * No change.

.. _changes_0.6.0:

`v0.6.0 <https://github.com/Ouranosinc/miranda/tree/v0.6.0>`_ (2025-09-03)
--------------------------------------------------------------------------
Contributors to this version: Travis Logan (:user:`tlogan2000`), Trevor James Smith (:user:`Zeitsperre`), Aslı Beşe (:user:`aslibese`).

Announcements
^^^^^^^^^^^^^
* `miranda` boilerplate code is now versioned with `cruft <https://cruft.github.io/cruft>`_ and the `Ouranosinc/cookiecutter-pypackage <https://github.com/Ouranosinc/cookiecutter-pypackage>`_ template.
* The `miranda` library now requires Python 3.10 or higher.
* The `miranda` library has undergone a significant refactoring to remove unused code and modules that were out of scope for the library.
    * The modules `miranda.archive` and `miranda.remote` have been removed and will be moved to a separate project.
    * Data treatment functions have been moved to `miranda.treatments`.
    * Pre-processing logic and functions have been moved to `miranda.preprocess`.
    * Several obsolete/nonexistent configuration files have been removed from `miranda.convert.corrections`.
* The `miranda.validators` submodule has become `miranda.validate` and now contains all controlled vocabulary and validation functions.
* `miranda` now leverages `pooch` to fetch and cache testing datasets from `miranda-testdata <https://github.com/Ouranosinc/miranda-testdata>`_:
    * ``miranda.testing.cassini`` is used to create an instance of `Cassini` for fetching testing data.
    * ``miranda.testing.registry.txt`` is a text file containing the list of datasets available in `miranda-testdata`.
* `miranda` has dropped `black`, `isort` and `blackdocs`, and now relies solely on `ruff` and `flake8`/`flake8-rst-docstrings` for linting and formatting.

New features
^^^^^^^^^^^^
* Conversion for `CanHomTv4 daily` ECCC 4th generation of the adjusted and homogenized station.
* Conversions for variables in the `GHCN-D` weather station network dataset are now supported.
* Conversion support has been added for the `ORRC`, `CaSR v3.1`, and `RDRS v2.1` datasets.
* ECMWF: Added support for the `"era5-single-levels-monthly-means"` project.
* ECMWF: Added support for ocean variables (`sst`, `siconc`), convective precipitation variables (`'cp'`, `'cape'`), and wind speed (`'u'`, `'v'`).
* Aggregation operations now support more variables (`'hur'`, `'hurs'`, `'huss'`, `'rlds'`, `'ta'`, `'tdp'`, `'ua'`, `'uas'`, `'va'`, `'vas'`)
* Minimum values of ``"0 kg m2 s-1"`` has been set for both `'tp'` and `'sf'` variables in ERA5 and ERA5-Land projects.
* Project user and developer documentation has been greatly expanded. All public functions and modules now have `numpy`-based docstrings.
* The `miranda` library now uses a `src` layout for better packaging and distribution.
* `ruff` checks and formatting standards have been adopted for the entire codebase.
* Added a new configuration for converting the NRCAN gridded climate dataset (`NRCANmet`).
* Conversion configuration JSON files are now validated against `schema <https://github.com/keleshev/schema>`_ schemas.

Breaking changes
^^^^^^^^^^^^^^^^
* Removed modules `miranda.archive` and `miranda.remote` (split into a separate project yet to be published).
* ``miranda.utils.show_versions`` has been moved to ``miranda.testing.show_versions``. It now uses ``xclim.testing.show_versions`` to display the versions of all dependencies.
* Python 3.8 and Python 3.9 are no longer supported.
* The `dev` recipe now requires `pooch` (>=1.8.0).
* many dependencies have been updated to more modern versions, including:
    * `numpy` (>=1.25.0)
    * `xarray` (>=2023.11.0)
    * `xclim` (>=0.57.0)
* Logging has been significantly improved and standardized across the library.
    * Logging within modules has been standardized to use the ``miranda`` logger and never ``root``.
    * Submodules no longer configure message logging to standard output and instead use the `miranda` logger.

Bug fixes
^^^^^^^^^
* Transformation docstrings are now only updated when the transformation is actually applied.
* Added a missing helper function to `miranda.units` (``group_by_length``) that was mistakenly removed in a previous change.

Internal changes
^^^^^^^^^^^^^^^^
* `miranda` now has a security policy (``SECURITY.md``) for disclosing sensitive issues using secure communication channels. This has also been added to the documentation.
* `miranda` now applies the `numpydoc` documentation style to all publicly-exposed docstrings.
* GitHub Workflows now use commit hashes for both running GitHub Actions and installing Python dependencies from PyPI.
* `miranda` now has a ``CODE_OF_CONDUCT.md`` file for setting community standards and expectations.
* Now using the GitHub Ouranos bot for automatic version bumping via ``bump-version.yml`` GitHub Workflow.
* Adjusted calls using `os.path` to use `pathlib` for better cross-platform compatibility.
* Added new `pytest` fixtures for the new `miranda-testdata` repository:
    * ``cassini``: `pytest` fixture for fetching local filepaths of cached testing data.
    * ``open_dataset``: `pytest` fixture for one-off fetching and opening of a registered test data set.
    * ``era5_precip``: `pytest` fixture fetching and opening a zip file containing a subset of the ERA5 precipitation dataset.
    * ``timeseries``: `pytest` fixture for generating an artificial CF-compliant time series dataset using `xclim` and `xarray`.
    * ``multivariable_dataset``: `pytest` fixture for generating an artificial `xarray` multivariable dataset.
* The ``tox.ini`` and ``pyproject.toml`` dependency pins have been synchronized.
* `schema` schemas have been defined for all conversion JSON files, and are now used to validate the JSON files as part of the testing suite.
* The code formatting now follows `ruff` standards, and `black`, `isort`, and `blackdocs` have been removed from the project. The `pre-commit` configuration has been updated accordingly. Line lengths have been increased from 88 to 150.
* `pre-commit` hook versions have been updated and new hooks have been added for checking variable spelling and security issues. Hooks for `mypy` and `vulture` have been staged for eventual inclusion in the CI testing suite.
* Allow some variables that are lacking a ``standard_name`` attribute to be converted if ``_standard_name`` is explicitly set as ``False``.

.. _changes_0.5.0:

`v0.5.0 <https://github.com/Ouranosinc/miranda/tree/v0.5.0>`_ (2023-06-19)
--------------------------------------------------------------------------
Contributors to this version: Juliette Lavoie (:user:`juliettelavoie`), Trevor James Smith (:user:`Zeitsperre`).

New features
^^^^^^^^^^^^
* Added support for collecting and converting `ptype` ECMWF ERA5 variable.
* A new ``"_frequency": true`` toggle for returning the output frequency of converted data.
* Added a new JSON template for NEX-GDDP-CMIP6 datasets.
* `miranda` is now `PEP 517 <https://peps.python.org/pep-0517/>`_ and `PEP 621 <https://peps.python.org/pep-0621/>`_ compliant, using the `flit <https://flit.pypa.io/en/stable/>`_ backend.

Internal changes
^^^^^^^^^^^^^^^^
* Various fixes to existing docstrings.
* Time frequency checks are more resilient when converting Monthly time-step data.
* Masking and regridding of datasets when running ``convert_dataset`` is now optional or automatic.
* Updated templates to newest API.
* Created a `gis` recipe for exclusively installing GIS libraries.
* Removed many unneeded dependencies, cleaned up Makefile.
* All public-facing functions now contain at least a minimal docstring for documentation generation.

.. _changes_0.4.0:

`v0.4.0 <https://github.com/Ouranosinc/miranda/tree/v0.4.0>`_ (2023-03-30)
--------------------------------------------------------------------------
Contributors to this version: Trevor James Smith (:user:`Zeitsperre`), Pascal Bourgault (:user:`aulemahal`), Travis Logan (:user:`tlogan2000`).

New features
^^^^^^^^^^^^
* Improvements have been made to the development documentation; Project URLs, ReadTheDocs theming, and other quality of life changes.
* Conversion JSON definitions now support pre-processing to render dimensions and variable names consistent before running corrections/conversions.
* New datasets with CF-like attributes conversion supported:
    - RDRS (ECCC)
    - GRNCH (ETS)
* Preliminary ``miranda.io`` module for organizing output-writing functionality.
* New ``miranda.io.fetch_chunk_config`` function for "rechunking" datasets according to project presets.
* New ``mirands.io.utils.name_output_file`` for generating names from Dataset facets or from a dictionary.
* New ``mirands.gis.subset_domain`` for clipping dataset to a preconfigured region.

Bug fixes
^^^^^^^^^
* Many data-related utilities now have more accurate static typing.
* Converted dataset global attributes are now synchronized for consistency.
* ECMWF-based datasets now implement more consistent conversion factors and metadata.
* ``miranda.storage.file_size`` now handles dictionaries of Pathlib objects.

Internal changes
^^^^^^^^^^^^^^^^
* Pre-commit version updates.
* Improvements have been made to the development documentation; Project URLs, ReadTheDocs theming, installation methods, and other quality of life changes.
* Schema and folder structure updates:
    - `gridded-obs` -> `reconstruction`
    - `bias-adjust-project` is used when present and not just when `level=="biasadjusted"`
* CI now using `tox>=4.0` and `ubuntu-latest` virtual machine images.

.. _changes_0.3.0:

`v0.3.0 <https://github.com/Ouranosinc/miranda/tree/v0.3.0>`_ (2022-11-24)
--------------------------------------------------------------------------
Contributors to this version: Trevor James Smith (:user:`Zeitsperre`), Pascal Bourgault (:user:`aulemahal`), David Huard (:user:`huard`), Travis Logan (:user:`tlogan2000`), Gabriel Rondeau-Genesse (:user:`RondeauG`), and Sébastien Biner (:user:`sbiner`).

Announcements
^^^^^^^^^^^^^
* First public release on PyPI.

New features
^^^^^^^^^^^^
* Dataset conversion tools (``miranda.convert``) use a JSON-definition file to dynamically populate metadata, run data quality checks, and convert units to CF-compliant standard. Supported datasets are:
    - ERA5/ERA5-Land (complete)
    - MELCC (stations) (beta)
    - ECCC (stations) (alpha)
    - NASA DayMet (WIP)
    - NASA AgMerra/AgCFSR (WIP)
    - Hydro Québec (stations) (WIP)
    - DEH (stations) (WIP)
    - WFDEI-GEM-CAPA (WIP)
* Module (``miranda.eccc``) for ECCC station data and ECCC Adjusted and Homogenized Canadian Climate Data (AHCCD) conversion (WIP).
* Module (``miranda.ncar``) for fetching interpolated CORDEX-NAM (22i/44i) from NCAR AWS data storage.
* Module (``miranda.ecmwf``) for fetching ECMWF ERA5/-Land (single-levels, pressure-levels, monthly-means) datasets via CDSAPI.
* Module (``miranda.gis``) for setting specific subsetting domains used when converting gridded datasets.
* Modules (``miranda.archive`` and ``miranda.remote``) for performing data archiving actions locally and remotely (powered by `fabric <https://github.com/fabric/fabric>`_ and `paramiko <https://github.com/paramiko/paramiko>`_) (WIP).
* Module (``miranda.decode``) for ingesting and parsing dataset metadata based on filename and dataset attributes. Supported datasets are:
    - `miranda` converted datasets
    - CMIP6
    - CMIP5
    - CMIP5-CORDEX
    - ISIMIP-FT
    - CanDCS-U6 (PCIC)
* Module (``miranda.structure``) for create constructing file-tree databases based on YAML-defined metadata schemas (WIP).
* Modules (``miranda.cv`` and ``miranda.validators``) for validating metadata using ESGF controlled vocabularies (taken from `pyessv-archive <https://github.com/ES-DOC/pyessv-archive>`_) and schema definitions (powered by `schema <https://github.com/keleshev/schema>`_), respectively (WIP).
