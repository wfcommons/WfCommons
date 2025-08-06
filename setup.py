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
import os
import stat
import shutil
import site
import pathlib
import subprocess

from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext


class Build(build_ext):
    """Customized setuptools build command - builds protos on build."""

    def run(self):
        protoc_command = ["make"]
        if subprocess.call(protoc_command) != 0:
            print("Error: 'make' is not installed. Please install 'make' and try again.")
            sys.exit(-1) 
        super().run()

        # Do a by-hand copy of cpu-benchmark to the bin directory, as this is so 
        # hard to do using data_file in the setup() declaration
        scripts_dir = os.path.join(site.USER_BASE, "bin")
        if not pathlib.Path(scripts_dir).is_dir():
            scripts_dir = "/opt/local/bin"
        source_path = os.path.join("bin", "cpu-benchmark")
        target_path = os.path.join(scripts_dir, "cpu-benchmark")
        # Ensure it's executable (just in case)
        st = os.stat(source_path)
        os.chmod(source_path, st.st_mode | stat.S_IEXEC)
        # Copy to scripts directory
        shutil.copy2(source_path, target_path)

setup(
    packages=find_packages(),
    include_package_data=True,
    has_ext_modules=lambda: True,
    cmdclass={
        'build_ext': Build,
    },
    #data_files=[
    #    (user_bin_dir, ['bin/cpu-benchmark'])
    #],
    scripts=['bin/wfbench'],
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
            'rnaseq_recipe = wfcommons.wfchef.recipes:RnaseqRecipe',
        ]
    },
)
