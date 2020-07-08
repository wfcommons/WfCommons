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
from logging import Logger

DEFAULT_LOGGER = logging.getLogger("workflowhub")


def configure_logger(logger) -> None:
	"""
	Logger configuration
	"""
	logger.setLevel(logging.DEBUG)

	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	formatter = logging.Formatter(
		'[%(asctime)s] %(levelname)-8s %(message)s',
		datefmt='%Y-%m-%d %H:%M:%S',
	)
	ch.setFormatter(formatter)
	logger.addHandler(ch)
