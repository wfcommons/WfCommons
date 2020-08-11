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

from typing import List, Optional
from logging import Logger

from .machine import Machine
from .file import File
from ..utils import NoValue


class JobType(NoValue):
    """Job type."""
    COMPUTE = 'compute'


class Job:
    """Representation of a job.

    :param name: The name of the job.
    :type name: str
    :param job_type: The type of the job.
    :type job_type: JobType
    :param runtime: Job runtime in seconds.
    :type runtime: float
    :param cores: Number of cores required by the job.
    :type cores: int
    :param machine: Machine on which is the job has been executed.
    :type machine: Machine
    :param args: List of job arguments.
    :type args: List[str]
    :param avg_cpu: Average CPU utilization in %.
    :type avg_cpu: float
    :param bytes_read: Total bytes read in KB.
    :type bytes_read: int
    :param bytes_written: Total bytes written in KB.
    :type bytes_written: int
    :param memory: Memory (resident set) size of the process in KB.
    :type memory: int
    :param energy: Total energy consumption in kWh.
    :type energy: int
    :param avg_power: Average power consumption in W.
    :type avg_power: float
    :param priority: Job priority.
    :type priority: int
    :param files: List of input/output files used by the job.
    :type files: List[File]
    :param logger: The logger where to log information/warning or errors.
    :type logger: Logger
    """

    def __init__(self,
                 name: str,
                 job_type: JobType,
                 runtime: float,
                 cores: int,
                 machine: Optional[Machine],
                 args: List[str],
                 avg_cpu: Optional[float],
                 bytes_read: Optional[int],
                 bytes_written: Optional[int],
                 memory: Optional[int],
                 energy: Optional[int],
                 avg_power: Optional[float],
                 priority: Optional[int],
                 files: List[File] = [],
                 logger: Optional[Logger] = None
                 ) -> None:
        """A job in a workflow."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        self.name: str = name
        self.type: JobType = job_type
        self.runtime: float = runtime
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

        self.logger.debug("created {0} job {1}: runtime => {2} secondes.".format(
            self.type, self.name, self.runtime)
        )
