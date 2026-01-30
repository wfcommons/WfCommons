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
import os

from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext

class Build(build_ext):
    """Customized setuptools build command - builds cpu-benchmark on build."""

    def run(self):
        # Try to build the cpu-benchmark, but make it optional
        # This allows installation on Windows where make/g++ may not be available
        try:
            result = subprocess.call(["make"], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
            if result != 0:
                sys.stderr.write("Warning: 'make' build failed. cpu-benchmark will not be available.\n")
                sys.stderr.write("This is expected on Windows. To build cpu-benchmark, install make and g++.\n")
        except (FileNotFoundError, OSError):
            sys.stderr.write("Warning: 'make' is not installed. cpu-benchmark will not be available.\n")
            sys.stderr.write("This is expected on Windows. To build cpu-benchmark, install make and g++.\n")
        super().run()

# Conditionally include cpu-benchmark if it exists
data_files = []
cpu_benchmark_path = 'bin/cpu-benchmark'
if True:
#if os.path.exists(cpu_benchmark_path):
    data_files.append(('bin', [cpu_benchmark_path]))

setup(
    packages=find_packages(),
    include_package_data=True,
    has_ext_modules=lambda: True,
    cmdclass={
        'build_ext': Build,
    },
    data_files=data_files,
    scripts=[
        'bin/wfbench'
    ],
    entry_points={
        'console_scripts': [
            'wfchef=wfcommons.wfchef.chef:main',
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
            'rnaseq_recipe = wfcommons.wfchef.recipes:RnaseqRecipe',
        ]
    },
)
