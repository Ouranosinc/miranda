.PHONY: help clean clean-pyc clean-build list test test-all coverage docs release sdist

help:
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "compliant - make changes to code to satisfy lint rules"
	@echo "lint - check style with modified pep8 and black"
	@echo "test - run tests quickly with the default Python"
	@echo "test-all - run tests on every Python version with tox"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	@echo "release - package and upload a release"
	@echo "sdist - package"
	@echo "install - install the library in developer mode"

clean: clean-build clean-pyc  ## remove all build, test, coverage and Python artifacts

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-docs:
	rm -f docs/miranda*.rst
	rm -f docs/modules.rst
	$(MAKE) -C docs clean

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

compliant:
	auotopep8 miranda tests
	black miranda tests

lint:
	flake8 miranda tests
	black  --check miranda tests

test:
	pytest tests

test-all:
	tox

coverage:
	coverage run --source miranda setup.py test
	coverage report -m
	coverage html
	open htmlcov/index.html

docs: clean-docs
	sphinx-apidoc -o docs/ --module-first miranda
	$(MAKE) -C docs linkcheck html
ifndef READTHEDOCS
	xdg-open docs/_build/html/index.html
endif

release: clean
	python setup.py sdist upload
	python setup.py bdist_wheel upload

sdist: clean
	python setup.py sdist
	python setup.py bdist_wheel upload
	ls -l dist

install: clean ## install the package to the active Python's site-packages
	pip install -e .
