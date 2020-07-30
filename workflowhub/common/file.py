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
from typing import Optional


class File:
	"""
		Representation of a file
	"""

	def __init__(self, name: str, size: int, link: str, logger: Optional[Logger] = None) -> None:
		"""
			A file uses by jobs

			:param name: the name of the file
			:type name: str
			:param name: File size in KB
			:type name: int
			:param name: Whether it is an input or output data (possible value: input or output)
			:type name: str
			:param logger: the logger where to log information/warning or errors
			:type logger: Logger
		"""

		if logger is None:
			self.logger: Logger = logging.getLogger(__name__)
		else:
			self.logger = logger

		self.name: str = name
		self.size: int = size
		self.link: str = link
