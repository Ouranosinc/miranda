============
Installation
============


Stable release
--------------

To install miranda, run this command in your terminal:

.. code-block:: console

    $ python -m pip install miranda

This is the preferred method to install miranda, as it will always install the most recent stable release.

To make use of remote operations (`miranda.remote`) and some dataset downloading functions (`miranda.ncar` `miranda.ecmwf`), additional libraries are needed.
They can can be installed with the following::

    $ pip install miranda[remote]

For better RAM usage when converting datasets, some additional/optional GIS libraries can be installed as well::

    $ pip install miranda[gis]

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: https://docs.python-guide.org/starting/installation/

From sources
------------

The sources for miranda can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git@github.com:Zeitsperre/miranda

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/Zeitsperre/miranda/tarball/main

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python -m pip install .


.. _Github repo: https://github.com/Zeitsperre/miranda
.. _tarball: https://github.com/Zeitsperre/miranda/tarball/main


Creating a Conda environment
----------------------------

To create a conda development environment including all miranda dependencies, enter the following command from within your cloned repo::

    $ conda create -n my_miranda_env python=3.8 --file=environment.yml
    $ conda activate my_miranda_env
    $ pip install -e .[dev]
