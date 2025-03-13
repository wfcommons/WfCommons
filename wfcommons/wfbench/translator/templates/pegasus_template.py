#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import os
from Pegasus.api import *


def which(file):
    for path in os.environ['PATH'].split(os.pathsep):
        if os.path.exists(os.path.join(path, file)):
            return os.path.join(path, file)
    return None

# catalogs
tc = TransformationCatalog()
rc = ReplicaCatalog()

t_cpu_benchmark = Transformation('cpu-benchmark', site='local',
pfn = os.getcwd() + '/bin/cpu-benchmark', is_stageable=True)
tc.add_transformations(t_cpu_benchmark)
transformation_path = os.getcwd() + '/bin/wfbench'

task_output_files = {}

# Generated code goes here