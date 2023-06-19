============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/Ouranosinc/miranda/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Miranda could always use more documentation, whether as part of the
official Miranda docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/Ouranosinc/miranda/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `miranda` for
local development.

1. Fork_ the `miranda` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/miranda.git

3. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

Now you can make your changes locally.

4. Begin by installing a development build of your branch::

    # To install miranda with its development environment dependencies
    $ pip install -e .[dev]
    # To install miranda with GIS libraries
    $ pip install -e .[gis]
    # To install miranda with its documentation dependencies
    $ pip install -e .[docs]
    # To install miranda with its remote API dependencies
    $ pip install -e .[remote]

5. When you're done making changes, check that your changes pass style and unit
   tests, including testing other Python versions with tox::

    $ tox

To get tox, just pip install it.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

.. _Fork: https://github.com/Ouranosinc/miranda/fork

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

#. The pull request should include tests.
#. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
#. The pull request should work for Python 3.7, 3.8, 3.9, 3.10, and 3.11.
   Check https://github.com/Ouranosinc/miranda/actions
   for active pull request builds or run the ``tox`` command and
   make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

	 $ pytest test/test_miranda.py
