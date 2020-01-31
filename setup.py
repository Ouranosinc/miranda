#!/usr/bin/env python
import os
import sys

from setuptools import find_packages
from setuptools import setup

NAME = "miranda"
VERSION = "0.1.0-beta"
DESCRIPTION = "Python utilities for climate data collection and management"
KEYWORDS = "climate meteorology archiving collection NetCDF"
URL = "https://github.com/Ouranosinc/miranda"
AUTHOR = "Trevor James Smith"
AUTHOR_EMAIL = "smith.trevorj@ouranos.ca"
REQUIRES_PYTHON = ">=3.5.0"
LICENSE = "Apache Software License 2.0"

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist bdist_wheel upload")
    sys.exit()

requirements = list()
with open("requirements.txt") as req:
    for dependency in req.readlines():
        requirements.append(dependency)

dev_requirements = list()
with open("requirements_dev.txt") as dev:
    for dependency in dev.readlines():
        dev_requirements.append(dependency)

docs_requirements = ["sphinx", "guzzle-sphinx-theme", "nbsphinx", "pandoc", "ipython"]

readme = open("README.rst").read()
doclink = """
Documentation
-------------

The full documentation is at http://miranda.readthedocs.io/en/latest."""
history = open("HISTORY.rst").read().replace(".. :changelog:", "")

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=readme + "\n\n" + doclink + "\n\n" + history,
    long_description_content_type="text/x-rst",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    packages=find_packages(),
    package_dir={"miranda": "miranda"},
    include_package_data=True,
    install_requires=requirements,
    extras_require={"docs": docs_requirements, "dev": dev_requirements},
    license=LICENSE,
    zip_safe=False,
    keywords="climate meteorology archiving collection NetCDF",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
)
