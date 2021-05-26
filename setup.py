#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from setuptools import setup, find_packages
from wfcommons.version import __version__

with open('README.md', 'r') as fh:
    long_description = fh.read()

# # Fetch the version
# exec(open('wfcommons/version.py').read())

setup(
    name='wfcommons',
    version=str(__version__),
    license='GPLv3',
    author='WfCommons team',
    author_email='support@wfcommons.org',
    description='Community Framework for Enabling Scientific Workflow Research and Education',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/wfcommons/wfcommons',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'jsonschema',
        'matplotlib',
        'networkx',
        'numpy',
        'pygraphviz',
        'python-dateutil',
        'requests',
        'scipy',
        'setuptools',
        'pyyaml'
    ],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Documentation :: Sphinx',
        'Topic :: System :: Distributed Computing'
    ],
    python_requires='>=3.6',
    entry_points = {
        'console_scripts': [
            'wfchef=wfcommons.wfchef.chef:main',
        ],
        'workflow_recipes': [
            # 'epigenomics_recipe = wfcommons.wfchef.wfcommons/.epigenomics:EpigenomicsRecipe',
            'epigenomics_recipe = wfcommons.wfchef.recipes/.epigenomics:EpigenomicsRecipe',
            'montage_recipe = wfcommons.wfchef.recipes/.montage:MontageRecipe',
            'cycles_recipe = wfcommons.wfchef.recipes/.cycles:CyclesRecipe',
            'seismology_recipe = wfcommons.wfchef.recipes/.seismology:SeismologyRecipe',
            'soykb_recipe = wfcommons.wfchef.recipes/.soykb:SoykbRecipe',
            'srasearch_recipe = wfcommons.wfchef.recipes/.srasearch:SrasearchRecipe',
            'genome_recipe = wfcommons.wfchef.recipes/.genome:GenomeRecipe',
            'blast_recipe = wfcommons.wfchef.recipes/.blast:BlastRecipe',
            'bwa_recipe = wfcommons.wfchef.recipes/.bwa:BwaRecipe',
        ]
    },
)
