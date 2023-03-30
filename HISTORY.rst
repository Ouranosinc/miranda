.. :changelog:

=======
History
=======

v0.4.0 (2023-03-30)
-------------------
Contributors to this version: Trevor James Smith (:user:`Zeitsperre`), Pascal Bourgault (:user:`aulemahal`),  Travis Logan (:user:`tlogan2000`).

New features
^^^^^^^^^^^^
* Improvements have been made to the development documentation; Project URLs, ReadTheDocs theming, and other quality of life changes.
* Conversion JSON definitions now support pre-processing to render dimensions and variable names consistent before running corrections/conversions.
* Dataset conversions now support ETS-GRNCH and ECCC RDRS projects.
* Preliminary ``miranda.io`` module for organizing output-writing functionality.
* New ``fetch_chunk_config`` function for "rechunking" datasets according to project presets.
* The ``name_output_file`` function now accepts a dictionary of facets for constructing filenames.

Internal changes
^^^^^^^^^^^^^^^^
* Pre-commit version updates.
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
