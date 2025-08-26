==============
Miranda |logo|
==============

+----------------------------+-----------------------------------------------------+
| Versions                   | |pypi|                                              |
+----------------------------+-----------------------------------------------------+
| Documentation and Support  | |docs| |versions|                                   |
+----------------------------+-----------------------------------------------------+
| Open Source                | |license| |ossf-score|                              |
+----------------------------+-----------------------------------------------------+
| Coding Standards           | |ruff| |pre-commit|                                 |
+----------------------------+-----------------------------------------------------+
| Development Status         | |status| |build| |coveralls|                        |
+----------------------------+-----------------------------------------------------+

Python utilities for climate data collection, conversion, and management.

* Free software: Apache Software License 2.0
* Documentation: https://miranda.readthedocs.io.

Features
--------

Data collection functions for climate and forecast data hosted at:
    * ECMWF (Europe) (ERA5, ERA5-Land, TIGGE)
    * ECCC (Canada) (Monthly Climate Summaries, ECCC GEOAPI) (In development)
    * NCAR (CORDEX-NA on AWS)

Data conversion for `Climate and Forecasting (CF) <https://cfconventions.org/>`_ Variable and Metadata compliance:
    * ECMWF (Europe) (ERA5, ERA5-Land, TIGGE) (In Development)
    * ECCC (Canada) (Flat File Observations, Monthly Climate Summaries, Adjusted and Homogenized Climate Data, ECCC GEOAPI) (In Development)
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
       * NRCAN (Canada) (In Development)
       * MELCC (Québec) (In Development)
    * Reanalyses:
       * ECCC (Canada) (CanSWE, CaSR, RDRS)
       * ECMWF (Europe) (ERA5, ERA5-Land, TIGGE)
       * NASA (DayMET, NEX-GDDP) (In Development)
       * USask (Canada) (EMDNA)

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

Credits
-------

This package was created with Cookiecutter_ and the `Ouranosinc/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/cookiecutter/cookiecutter
.. _`Ouranosinc/cookiecutter-pypackage`: https://github.com/Ouranosinc/cookiecutter-pypackage

.. |build| image:: https://github.com/Ouranosinc/miranda/actions/workflows/main.yml/badge.svg
        :target: https://github.com/Ouranosinc/miranda/actions
        :alt: Build Status

..
    .. |conda| image:: https://img.shields.io/conda/vn/conda-forge/miranda.svg
            :target: https://anaconda.org/conda-forge/miranda
            :alt: Conda-forge Build Version

.. |coveralls| image:: https://coveralls.io/repos/github/Ouranosinc/miranda/badge.svg
        :target: https://coveralls.io/github/Ouranosinc/miranda
        :alt: Coveralls

.. |docs| image:: https://readthedocs.org/projects/miranda/badge/?version=latest
        :target: https://miranda.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

.. |license| image:: https://img.shields.io/github/license/Ouranosinc/miranda.svg
        :target: https://github.com/Ouranosinc/miranda/blob/main/LICENSE
        :alt: License

..
    .. |ossf-bp| image:: https://bestpractices.coreinfrastructure.org/projects/9945/badge
            :target: https://bestpractices.coreinfrastructure.org/projects/9945
            :alt: Open Source Security Foundation Best Practices

.. |ossf-score| image:: https://api.securityscorecards.dev/projects/github.com/Ouranosinc/miranda/badge
        :target: https://securityscorecards.dev/viewer/?uri=github.com/Ouranosinc/miranda
        :alt: OpenSSF Scorecard

.. |logo| image:: https://raw.githubusercontent.com/Ouranosinc/miranda/main/docs/_static/images/miranda-logo-small.png
        :target: https://github.com/Ouranosinc/miranda
        :alt: Miranda

.. |pre-commit| image:: https://results.pre-commit.ci/badge/github/Ouranosinc/miranda/main.svg
        :target: https://results.pre-commit.ci/latest/github/Ouranosinc/miranda/main
        :alt: pre-commit.ci status

.. |pypi| image:: https://img.shields.io/pypi/v/miranda.svg
        :target: https://pypi.python.org/pypi/miranda
        :alt: PyPI

.. |ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
        :target: https://github.com/astral-sh/ruff
        :alt: Ruff

.. |status| image:: https://www.repostatus.org/badges/latest/active.svg
        :target: https://www.repostatus.org/#active
        :alt: Project Status: Active – The project has reached a stable, usable state and is being actively developed.

.. |versions| image:: https://img.shields.io/pypi/pyversions/miranda.svg
        :target: https://pypi.python.org/pypi/miranda
        :alt: Supported Python Versions
