#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from setuptools import setup, find_packages
import pathlib

thisdir = pathlib.Path(__file__).resolve().parent

workflow_recipes = [
    line.strip()
    for line in thisdir.joinpath("workflow_recipes.txt").read_text().splitlines()
    if line.strip()
]

setup(
    name='wfcommons.wfchef.recipe.PACKAGE_NAME',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'networkx',
        'wfcommons'
    ],
    entry_points={
        'workflow_recipes': workflow_recipes,
    }
)
