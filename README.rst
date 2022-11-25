==============
Miranda |logo|
==============

|build| |coveralls| |black|

Python utilities for climate data collection, conversion, and management

* Documentation: |docs|
* Free Software: |license|

Features
--------

Data collection functions for climate and forecast data hosted at:
    * ECMWF (ERA5, ERA5-Land, TIGGE)
    * ECCC (Canada) (Monthly Climate Summaries, ECCC GEOAPI - In development)
    * NCAR (CORDEX-NA on AWS)

Data conversion for `Climate and Forecasting (CF) <https://cfconventions.org/>`_ Variable and Metadata compliance:
    * ECMWF (ERA5, ERA5-Land, TIGGE - In Development)
    * ECCC (Canada) (Flat File Observations, Monthly Climate Summaries, Adjusted and Homogenized Climate Data, ECCC GEOAPI - In Development)
    * MELCC (Québec)
    * Hydro-Québec (In Development)

Database structuring and facets validation:
    * Simulations:
       * WCRP (CMIP5, CMIP6, CMIP5-CORDEX, CORDEX-ADJUST, ISIMIP, etc.)
    * Station-Observations:
       * MELCC (Québec) (Needs `mdbtools <https://github.com/mdbtools/mdbtools>`_ installed)
       * ECCC (Canada) (In Development)
       * Hydro-Québec (In Development)
    * Gridded-Observations:
       * NRCAN (Canada) (Future)
       * MELCC (Future)
    * Reanalyses:
       * ECMWF (ERA5, ERA5-Land, TIGGE)
       * NASA (DayMET, AgMerra/AgCFSR, MERRA2) - In Development
       * NCEP (CFSR/CFSv2) - In Development
       * WFDEI-GEM-CaPa (University of Saskatchewan) - In Development

Installation
------------
`miranda` can be installed from PyPI::

    $ pip install miranda

Some functionalities require complex-to-install dependencies.
In order to gain access to them, we strongly suggest using `Anaconda <https://www.anaconda.com/products/distribution>`_ to manage your environment::

    $ conda env create -f environment.yml
    $ conda activate miranda
    $ pip install miranda[full]

For more information about Anaconda/Miniconda/conda-forge:
 * Miniconda: https://docs.conda.io/en/latest/miniconda.html
 * conda-forge: https://conda-forge.org/#about

`miranda` also relies on `PyESSV <https://github.com/ES-DOC/pyessv>`_ for its climate data controlled vocabulary.
This library is optional for users who do not require validation checks,
but enabling this feature requires additional installation steps::


    $ mkdir -p ~/.esdoc
    $ git clone git@github.com:ES-DOC/pyessv-archive.git ~/.esdoc/pyessv-archive

Contributing
------------
See the contributing documentation: https://miranda.readthedocs.io/en/latest/contributing.html

.. |build| image:: https://github.com/Ouranosinc/miranda/actions/workflows/main.yml/badge.svg
        :target: https://github.com/Ouranosinc/miranda/actions/workflows/main.yml
        :alt: Build Status

.. |coveralls| image:: https://coveralls.io/repos/github/Ouranosinc/miranda/badge.svg
        :target: https://coveralls.io/github/Ouranosinc/miranda
        :alt: Coveralls

.. |docs| image:: https://readthedocs.org/projects/miranda/badge
        :target: https://miranda.readthedocs.io/en/latest
        :alt: Documentation Status

.. |license| image:: https://img.shields.io/github/license/Ouranosinc/miranda.svg
        :target: https://github.com/Ouranosinc/miranda/blob/master/LICENSE
        :alt: License

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/psf/black
        :alt: Python Black

.. |logo| image:: https://raw.githubusercontent.com/Ouranosinc/miranda/main/docs/_static/images/miranda-logo-small.png
        :target: https://github.com/Ouranosinc/miranda
        :alt: Miranda
