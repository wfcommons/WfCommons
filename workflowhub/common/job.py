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
from typing import List, Dict, Any, Union, Optional
from logging import Logger

from .machine import Machine
from .file import File


class Job:
	"""
		Representation of a job
	"""

	def __init__(self, name: str, jtype: str, runtime: int, cores: int, machine: Optional[Machine], args: List[str],
				 avg_cpu: Optional[float], bytes_read: Optional[int], bytes_written: Optional[int],
				 memory: Optional[int], energy: Optional[int], avg_power: Optional[float],
				 priority: Optional[int],
				 files: List[File] = [],
				 logger: Optional[Logger] = None
				 ) -> None:
		"""
			A job in a workflow.

			:param name: the name of the job
			:type name: str
			:param jtype: the type of the job (whether it is a compute or an auxiliary job)
			:type jtype: str
			:param runtime: Job runtime in seconds
			:type runtime: int
			:param cores: Number of cores required by the job
			:type cores: int
			:param machine: Machine on which is the job has been executed
			:type machine: Machine
			:param args: List of job arguments
			:type args: List[str]
			:param avg_cpu: Average CPU utilization in %
			:type avg_cpu: Optional[float]
			:param bytes_read: Total bytes read in KB
			:type bytes_read: Optional[int]
			:param bytes_written: Total bytes written in KB
			:type bytes_written: Optional[int]
			:param memory: Memory (resident set) size of the process in KB
			:type memory: Optional[int]
			:param energy: Total energy consumption in kWh
			:type energy: Optional[int]
			:param avg_energy: Average power consumption in W
			:type avg_energy: Optional[float]
			:param priority: 
			:type priority: Optional[int]
			:param files: List of input/output files used by the job
			:type files: List[File]
			:param logger: the logger where to log information/warning or errors
			:type logger: Logger
		"""

		if logger is None:
			self.logger: Logger = logging.getLogger(__name__)
		else:
			self.logger = logger

		self.name: str = name
		self.type: str = jtype
		self.runtime: int = runtime
		self.cores: Optional[int] = cores
		self.args: List[str] = args
		self.avg_cpu: Optional[float] = avg_cpu
		self.bytes_read: Optional[int] = bytes_read
		self.bytes_written: Optional[int] = bytes_written
		self.memory: Optional[int] = memory
		self.energy: Optional[int] = energy
		self.avg_power: Optional[float] = avg_power
		self.files: List[File] = files
		self.machine: Machine = machine

		self.logger.info("created {0} job {1}: runtime => {2} secondes.".format(
			self.type, self.name, self.runtime)
		)
