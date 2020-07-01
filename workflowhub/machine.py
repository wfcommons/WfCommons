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
from typing import Dict, Any, Union, Optional
from logging import Logger

class Machine():
	"""
		Representation of one compute machine
	"""
	def __init__(self,
			name: str,
			cpu: Dict[str, Union[int, str]],
			system: Optional[str],
			architecture: Optional[str],
			memory: Optional[int],
			release: Optional[str],
			hashcode: Optional[str],
			logger: Optional[Logger] = None
		) -> None:
		"""
			A workflow trace.

			:param name: Machine node name
			:type name: Optional[str]
			:param cpu: a dictionnary containing information about the CPU specification. 
				- Must at least contains two fields: 
					- count: int (Number of CPU cores)
					- speed: int (CPU speed of each core in MHz)
			:type cpu: Dict[str, Union[int, str]]
			:param system: Machine system (linux, macos, windows)
			:type system: Optional[str]
			:param architecture: Machine architecture (e.g., x86_64, ppc)
			:type architecture: Optional[str]
			:param memory: Total machine's RAM memory in KB
			:type memory: Optional[int]
			:param release: Machine release
			:type release: Optional[str]
			:param hashcode: MD5 Hashcode for the Machine
			:type hashcode: Optional[str]
			:param logger: the logger where to log information/warning or errors
			:type logger: Optional[Logger]
		"""
		if logger is None:
			self.logger: Logger = logging.getLogger(__name__)
		else:
			self.logger = logger

		self.name: str = name
		self.cpu: Dict[str, Union[int, str]] = cpu
		self.system: str = system
		self.architecture: str = architecture
		self.memory: int = memory
		self.release: str = release
		self.hashcode = hashcode

		self.cores: int = cpu['count']
		self.speed: int = cpu['speed']
		self.flops: int = cpu['count']*cpu['speed']*10^6

		self.logger.info("created machine: {0} with {1} cores and {2} FLOPS.".format(
			self.name, self.cores, self.flops)
		)
