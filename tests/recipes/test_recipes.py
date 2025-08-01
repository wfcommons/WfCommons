#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import pytest
import wfcommons.utils

from typing import Dict, List, Tuple
from wfcommons.wfchef.recipes import GenomeRecipe
from wfcommons.wfchef.recipes import SeismologyRecipe
from wfcommons.wfchef.recipes import MontageRecipe
from wfcommons.wfchef.recipes import RnaseqRecipe
from wfcommons.wfchef.recipes import BwaRecipe
from wfcommons.wfchef.recipes import SoykbRecipe
from wfcommons.wfchef.recipes import CyclesRecipe
from wfcommons.wfchef.recipes import BlastRecipe
from wfcommons.wfchef.recipes import EpigenomicsRecipe
from wfcommons.wfchef.recipes import SrasearchRecipe

from wfcommons import WorkflowGenerator

class TestRecipes:

    recipe_class_list = [SeismologyRecipe,
            MontageRecipe,
            RnaseqRecipe,
            BwaRecipe,
            SoykbRecipe,
            CyclesRecipe,
            BlastRecipe,
            EpigenomicsRecipe,
            SrasearchRecipe,
        ]
    
    # @pytest.mark.unit
    # @pytest.mark.parametrize(
    #     "recipe_class",
    #     recipe_class_list
    # )
    # def test_recipes(self, recipe_class) -> None:
    #
    #     recipe = recipe_class.from_num_tasks(num_tasks=200, runtime_factor=1.1, input_file_size_factor=1.5,
    #                                              output_file_size_factor=0.8)
    #     workflow = WorkflowGenerator(recipe).build_workflow()
    #
    #
    # @pytest.mark.unit
    # @pytest.mark.parametrize(
    #     "recipe_class",
    #     recipe_class_list
    # )
    # def test_recipes_errors(self, recipe_class) -> None:
    #     # Not enough tasks
    #     recipe = recipe_class.from_num_tasks(num_tasks=2, runtime_factor=1.1, input_file_size_factor=1.5,
    #                                          output_file_size_factor=0.8)
    #     with pytest.raises(ValueError):
    #         WorkflowGenerator(recipe).build_workflow()
    #
    #     # Bogus parameters
    #     with pytest.raises(ValueError):
    #         recipe_class.from_num_tasks(num_tasks=2, runtime_factor=-1.1, input_file_size_factor=1.5,
    #                                     output_file_size_factor=0.8)
    #     with pytest.raises(ValueError):
    #         recipe_class.from_num_tasks(num_tasks=2, runtime_factor=1.1, input_file_size_factor=-1.5,
    #                                     output_file_size_factor=0.8)
    #     with pytest.raises(ValueError):
    #         recipe_class.from_num_tasks(num_tasks=2, runtime_factor=1.1, input_file_size_factor=-1.5,
    #                                     output_file_size_factor=-0.8)
