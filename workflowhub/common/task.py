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


class TaskType(NoValue):
    """Task type."""
    COMPUTE = 'compute'


class Task:
    """Representation of a task.

    :param name: The name of the task.
    :type name: str
    :param task_type: The type of the task.
    :type task_type: TaskType
    :param runtime: Task runtime in seconds.
    :type runtime: float
    :param cores: Number of cores required by the task.
    :type cores: int
    :param machine: Machine on which is the task has been executed.
    :type machine: Machine
    :param args: List of task arguments.
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
    :param priority: Task priority.
    :type priority: int
    :param files: List of input/output files used by the task.
    :type files: List[File]
    :param logger: The logger where to log information/warning or errors.
    :type logger: Logger
    """

    def __init__(self,
                 name: str,
                 task_type: TaskType,
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
        """A task in a workflow."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        self.name: str = name
        self.type: TaskType = task_type
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

        self.logger.debug("created {0} task {1}: runtime => {2} secondes.".format(
            self.type, self.name, self.runtime)
        )
