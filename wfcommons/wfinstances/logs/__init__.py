#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from .makeflow import MakeflowLogsParser
from .taskvine import TaskVineLogsParser
from .nextflow import NextflowLogsParser
from .pegasus import PegasusLogsParser
from .pegasusrec import HierarchicalPegasusLogsParser
