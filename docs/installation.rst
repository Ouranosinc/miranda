============
Installation
============

At the command line either via pip::

    $ pip install miranda

This is the preferred method to install miranda, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

From sources
------------
The sources for miranda can be downloaded from the `Github repo`_.

You can either clone the public repository::

    $ git clone git://github.com/Zeitsperre/miranda

Or download the `tarball`_::

    $ curl -OL https://github.com/Zeitsperre/miranda/tarball/master

Once you have a copy of the source, you can install it with::

    $ python setup.py install

Alternatively, you can also install a local copy via pip::

    $ pip install .

.. _Github repo: https://github.com/Ouranosinc/xclim
.. _tarball: https://github.com/Ouranosinc/xclim/tarball/master

Creating a Conda environment
----------------------------

To create a conda development environment including all xclim dependencies, enter the following command from within your cloned repo::

    $ conda create -n my_miranda_env python=3.6 --file=requirements_dev.txt
    $ conda activate my_miranda_env
    $ python setup.py install
