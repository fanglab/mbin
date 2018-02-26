#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

globals = {}
execfile("mbin/__init__.py", globals)
__version__ = globals["__version__"]
__author__ = globals["__author__"]
__email__ = globals["__email__"]

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'pysam >= 0.10.0',
    'numpy >= 1.7.1, < 1.14',
    'h5py >= 2.0.1',
    'pbcore >= 0.9.4',
    'scipy >= 0.12.0',
    'biopython >= 1.6.1',
    'matplotlib >= 1.5.0'
]

setup_requirements = [
    # TODO(jbeaulaurier): put setup requirements (distutils extensions, etc.) here
    "nose",
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='mbin',
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="mBin: a methylation-based binning framework for metagenomic SMRT sequencing reads",
    long_description=readme + '\n\n' + history,
    url='https://github.com/fanglab/mbin',
    packages=find_packages(include=['mbin']),
    include_package_data=True,
    # data_files=[('', ['bin/bh_tsne'])],
    install_requires=requirements,
    license="CC BY-NC-SA 4.0 license",
    zip_safe=False,
    keywords='mbin methylation sequencing binning metagenomic',
    # scripts=["bin/mbin"],
    entry_points={"console_scripts" : ["buildcontrols = mbin.controls:launch", \
                                       "filtermotifs = mbin.motifs:launch", \
                                       "methylprofiles = mbin.profiles:launch", \
                                       "mapfeatures = mbin.visualize:launch", \
                                       ]},
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: Free for non-commercial use',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='nose.collector',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
