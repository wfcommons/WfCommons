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

setup(
    packages=find_packages(),
    include_package_data=True,
    has_ext_modules=lambda: True,
    cmdclass={
        'build_ext': Build,
    },
    data_files=[
        ('bin', ['bin/cpu-benchmark', 'bin/wfbench'])
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
