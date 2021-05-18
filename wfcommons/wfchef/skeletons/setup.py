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

setup(
    name='wfchef.recipe.skeleton',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'networkx',
        'wfcommons'
    ],
    entry_points = {
        'worfklow_recipes': [
            "skeleton_recipe = skeleton_recipe:SkeletonRecipe"
        ],
    }
)