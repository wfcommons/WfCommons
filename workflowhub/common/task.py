#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging

from typing import Dict, List, Optional
from logging import Logger

from .machine import Machine
from .file import File
from ..utils import NoValue


class TaskType(NoValue):
    """Task type."""
    COMPUTE = 'compute'
    AUXILIARY = 'auxiliary'
    TRANSFER = 'transfer'


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
                 machine: Optional[Machine] = None,
                 args: List[str] = [],
                 avg_cpu: Optional[float] = None,
                 bytes_read: Optional[int] = None,
                 bytes_written: Optional[int] = None,
                 memory: Optional[int] = None,
                 energy: Optional[int] = None,
                 avg_power: Optional[float] = None,
                 priority: Optional[int] = None,
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
        self.priority = priority

        self.logger.debug("created {0} task {1}: runtime => {2} seconds.".format(
            self.type, self.name, self.runtime)
        )

    def as_dict(self) -> Dict:
        """A JSON representation of the task.

        :return: A JSON object representation of the task.
        :rtype: Dict
        """
        task_files = []
        for f in self.files:
            task_files.append(f.as_dict())

        task_obj = {
            'name': self.name,
            'type': self.type.value,
            'runtime': self.runtime,
            'parents': [],
            'children': [],
            'files': task_files,
        }
        if self.cores:
            task_obj['cores'] = self.cores
        if self.avg_cpu:
            task_obj['avgCPU'] = self.avg_cpu
        if self.bytes_read:
            task_obj['bytesRead'] = self.bytes_read
        if self.bytes_written:
            task_obj['bytesWritten'] = self.bytes_written
        if self.memory:
            task_obj['memory'] = self.memory
        if self.energy:
            task_obj['energy'] = self.energy
        if self.avg_power:
            task_obj['avgPower'] = self.avg_power
        if self.priority:
            task_obj['priority'] = self.priority
        if self.args:
            task_obj['arguments'] = self.args
        if self.machine:
            task_obj['machine'] = self.machine.name

        return task_obj
