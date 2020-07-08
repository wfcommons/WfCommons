#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

__version__ = "1.0"
__author__ = 'WorkflowHub Team'

from .logger import DEFAULT_LOGGER, configure_logger

# Configure the default logger 'workflowhub'
configure_logger(DEFAULT_LOGGER)

from .errors import TraceNotValid
from .types import JsonDict
from .utils import read_json
from .file import File
from .machine import Machine
from .job import Job
from .trace import Trace
