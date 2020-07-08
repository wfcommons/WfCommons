#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging
import json
from logging import Logger
from typing import Optional

from .types import JsonDict

"""
I/O from/to JSON traces
"""


def read_json(trace: str, logger: Optional[Logger] = None) -> JsonDict:
	"""
		Reads the json from the file path.
		returns the json object loaded with json data from the file
		:param trace:
		:param logger: the logger uses to output debug information
		:return: json object
	"""
	if logger is None:
		logger = logging.getLogger("workflowhub")

	with open(trace) as data:
		data_json: JsonDict = json.load(data)
		logger.info("parsed JSON trace: " + trace)
		return data_json
