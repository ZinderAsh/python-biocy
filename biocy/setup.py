#!/usr/bin/env python

"""The setup script."""

# distutils: language = c++

from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['cython', ]

test_requirements = ['pytest>=3', ]

extensions = [
    Extension("Graph",
              ["biocy/Graph.pyx", "biocy/cpp/Graph.cpp", "biocy/cpp/KmerFinder.cpp", "biocy/cpp/GFA.cpp", "biocy/cpp/utils.cpp"],
              include_dirs=[numpy.get_include()]),
]

extensions = cythonize(extensions, annotate=True)

setup(
    ext_modules=extensions,
    author="Sindre Ask Vestaberg",
    author_email='sindre.ask.vestaberg@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Temporary Description.",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='biocy',
    name='biocy',
    packages=find_packages(include=['biocy', 'biocy.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/zinderash/biocy',
    version='0.1.0',
    zip_safe=False,
)
