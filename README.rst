========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |requires|
        | |coveralls| |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/python-snakewrap/badge/?style=flat
    :target: https://readthedocs.org/projects/python-snakewrap
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/dohlee/python-snakewrap.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/dohlee/python-snakewrap

.. |requires| image:: https://requires.io/github/dohlee/python-snakewrap/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/dohlee/python-snakewrap/requirements/?branch=master

.. |coveralls| image:: https://coveralls.io/repos/dohlee/python-snakewrap/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/dohlee/python-snakewrap

.. |codecov| image:: https://codecov.io/github/dohlee/python-snakewrap/coverage.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/dohlee/python-snakewrap

.. |version| image:: https://img.shields.io/pypi/v/snakewrap.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/snakewrap

.. |commits-since| image:: https://img.shields.io/github/commits-since/dohlee/python-snakewrap/v0.0.0.svg
    :alt: Commits since latest release
    :target: https://github.com/dohlee/python-snakewrap/compare/v0.0.0...master

.. |wheel| image:: https://img.shields.io/pypi/wheel/snakewrap.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/snakewrap

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/snakewrap.svg
    :alt: Supported versions
    :target: https://pypi.org/project/snakewrap

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/snakewrap.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/snakewrap


.. end-badges

Utility package for snakemake wrappers.

* Free software: MIT license

Installation
============

::

    pip install snakewrap

Documentation
=============


https://python-snakewrap.readthedocs.io/


Development
===========

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
