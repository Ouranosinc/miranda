==============
Miranda |logo|
==============

|build| |coveralls| |codefactor| |black|

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
    * MELCC (Québec) (In Development)
    * Hydro-Québec (In Development)

Database structuring and facets validation:
    * Simulations:
       * WCRP (CORDEX, CORDEX-ADJUST, CMIP5, CMIP6, ISIMIP, etc.)
    * Station-Observations:
       * ECCC (Canada) (In Development)
       * MELCC (Québec) (In Development)
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
`miranda` is not yet available on PyPI, so the suggested method to install is as follows::

    $ git clone git@github.com:Ouranosinc/miranda.git
    $ cd miranda

    # If using Anaconda:
    $ conda create -n miranda -f environment.yml
    $ conda activate miranda

    $ pip install miranda[full]

`miranda` also relies on `PyESSV <https://github.com/ES-DOC/pyessv>`_ for its climate data controlled vocabulary. This library requires additional installation steps::

    $ mkdir -p ~/.esdoc
    $ git clone git@github.com:ES-DOC/pyessv-archive.git ~/.esdoc/pyessv-archive

*We strongly suggest using Anaconda3/miniconda3 (with the conda-forge repository enabled) to manage your environment and dependencies*
 * Anaconda: https://www.anaconda.com/products/distribution
 * Miniconda: https://docs.conda.io/en/latest/miniconda.html
 * conda-forge: https://conda-forge.org/#about

Contributing
------------
See the contributing documentation: https://miranda.readthedocs.io/en/latest/contributing.html

.. |build| image:: https://github.com/Ouranosinc/miranda/actions/workflows/main.yml/badge.svg
        :target: https://github.com/Ouranosinc/miranda/actions/workflows/main.yml
        :alt: Build Status

.. |coveralls| image:: https://coveralls.io/repos/github/Ouranosinc/miranda/badge.svg
        :target: https://coveralls.io/github/Ouranosinc/miranda
        :alt: Coveralls

.. |codefactor| image:: https://www.codefactor.io/repository/github/Ouranosinc/miranda/badge
        :target: https://www.codefactor.io/repository/github/Ouranosinc/miranda
        :alt: CodeFactor

.. |docs| image:: https://readthedocs.org/projects/miranda/badge
        :target: https://miranda.readthedocs.io/en/latest
        :alt: Documentation Status

.. |license| image:: https://img.shields.io/github/license/Ouranosinc/miranda.svg
        :target: https://github.com/Ouranosinc/miranda/blob/master/LICENSE
        :alt: License

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/psf/black
        :alt: Python Black

.. |logo| image:: https://raw.githubusercontent.com/Ouranosinc/miranda/add_logo/docs/_static/images/miranda-logo-small.png
        :target: https://github.com/Ouranosinc/miranda
        :alt: Miranda
