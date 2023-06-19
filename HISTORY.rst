.. :changelog:

=======
History
=======

v0.5.0 (2023-06-19)
-------------------
Contributors to this version: Juliette Lavoie (:user: `juliettelavoie`), Trevor James Smith (:user:`Zeitsperre`).

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

v0.4.0 (2023-03-30)
-------------------
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

v0.3.0 (2022-11-24)
-------------------
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
