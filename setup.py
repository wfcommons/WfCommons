#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import sys
import subprocess

from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext


class Build(build_ext):
    """Customized setuptools build command - builds protos on build."""

    def run(self):
        protoc_command = ["make"]
        if subprocess.call(protoc_command) != 0:
            sys.exit(-1)
        super().run()


with open('README.md', 'r') as fh:
    long_description = fh.read()

# Fetch the version
exec(open('wfcommons/version.py').read())

setup(
    name='wfcommons',
    version=str(__version__),
    license='GPLv3',
    author='WfCommons team',
    author_email='support@wfcommons.org',
    description='A Framework for Enabling Scientific Workflow Research and Education',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/wfcommons/wfcommons',
    packages=find_packages(),
    include_package_data=True,
    has_ext_modules=lambda: True,
    cmdclass={
        'build_ext': Build,
    },
    install_requires=[
        'jsonschema',
        'matplotlib',
        'networkx',
        'numpy',
        'python-dateutil',
        'requests',
        'scipy',
        'setuptools',
        'pyyaml',
        'pandas',
        'stringcase'
    ],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Documentation :: Sphinx',
        'Topic :: System :: Distributed Computing'
    ],
    python_requires='>=3.8',
    data_files=[
        ('bin', ['bin/cpu-benchmark'])
    ],
    entry_points={
        'console_scripts': [
            'wfchef=wfcommons.wfchef.chef:main'
        ],
        'workflow_recipes': [
            'epigenomics_recipe = wfcommons.wfchef.recipes:EpigenomicsRecipe',
            'montage_recipe = wfcommons.wfchef.recipes:MontageRecipe',
            'cycles_recipe = wfcommons.wfchef.recipes:CyclesRecipe',
            'seismology_recipe = wfcommons.wfchef.recipes:SeismologyRecipe',
            'soykb_recipe = wfcommons.wfchef.recipes:SoykbRecipe',
            'srasearch_recipe = wfcommons.wfchef.recipes:SrasearchRecipe',
            'genome_recipe = wfcommons.wfchef.recipes:GenomeRecipe',
            'blast_recipe = wfcommons.wfchef.recipes:BlastRecipe',
            'bwa_recipe = wfcommons.wfchef.recipes:BwaRecipe',
        ]
    },
)
