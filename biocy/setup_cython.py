#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
from Cython.Build import cythonize

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

setup(
    ext_modules = cythonize("biocy/Graph.pyx", annotate=True)
)

