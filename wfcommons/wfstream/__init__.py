#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Streaming-workflow support: dispel4py traces <-> WfFormat, and recipes built
from them.

`pipeline.on_new_workflow` / `pipeline.on_new_size_run` are the entry points a
registry calls; the other modules are the individual steps.
"""

from .dispel_fwd_converter import build_workflow, convert
from .pipeline import on_new_size_run, on_new_workflow

# `simulate` and the other steps stay module-level: exporting the function here
# would shadow the module of the same name for `from wfcommons.wfstream import
# simulate`.
