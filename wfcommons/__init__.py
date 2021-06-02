#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from .version import __version__

__author__ = 'WfCommons Team - https://wfcommons.org'
__credits__ = 'University of Southern California, University of Hawaii at Manoa'

import logging

from .wfchef.recipes import BlastRecipe, BwaRecipe, CyclesRecipe, EpigenomicsRecipe, GenomeRecipe, MontageRecipe, \
    SeismologyRecipe, SoykbRecipe, SrasearchRecipe
from .wfgen import WorkflowGenerator
from .wfinstances import Instance, InstanceAnalyzer, InstanceElement

logging.getLogger('wfcommons').addHandler(logging.NullHandler())
