#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2023 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging

from datetime import datetime
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
    SUBWORKFLOW = 'subworkflow'

class Task:
    """Representation of a task.

    :param name: The name of the task.
    :type name: str
    :param task_id: Task unique ID (e.g., ID0000001).
    :type task_id: str
    :param runtime: Task runtime in seconds.
    :type runtime: float
    :param input_files: List of input files used by the task.
    :type input_files: Optional[List[File]]
    :param output_files: List of output files used by the task.
    :type output_files: Optional[List[File]]
    :param cores: Number of cores required by the task.
    :type cores: float
    :param category: Task category (can be used, for example, to define tasks that use the same program).
    :type category: Optional[str]
    :param machines: Machines on which is the task has been executed.
    :type machines: Optional[List[Machine]]
    :param program: Program name.
    :type program: Optional[str]
    :param args: List of task arguments.
    :type args: Optional[List[str]]
    :param avg_cpu: Average CPU utilization in %.
    :type avg_cpu: Optional[float]
    :param bytes_read: Total bytes read in KB.
    :type bytes_read: Optional[int]
    :param bytes_written: Total bytes written in KB.
    :type bytes_written: Optional[int]
    :param memory: Memory (resident set) size of the process in bytes.
    :type memory: Optional[int]
    :param energy: Total energy consumption in kWh.
    :type energy: Optional[int]
    :param avg_power: Average power consumption in W.
    :type avg_power: Optional[float]
    :param priority: Task priority.
    :type priority: Optional[int]
    :param task_type: The type of the task.
    :type task_type: TaskType
    :param logger: The logger where to log information/warning or errors.
    :type logger: Optional[Logger]
    """

    def __init__(self,
                 name: str,
                 task_id: str,
                 runtime: float,
                 input_files: Optional[List[File]] = None,
                 output_files: Optional[List[File]] = None,
                 cores: float = 1.0,
                 category: Optional[str] = None,
                 machines: Optional[List[Machine]] = None,
                 program: Optional[str] = None,
                 args: Optional[List[str]] = None,
                 avg_cpu: Optional[float] = None,
                 bytes_read: Optional[int] = None,
                 bytes_written: Optional[int] = None,
                 memory: Optional[int] = None,
                 energy: Optional[int] = None,
                 avg_power: Optional[float] = None,
                 priority: Optional[int] = None,
                 executedAt: Optional[str] = None,
                 task_type: Optional[TaskType] = None,
                 launch_dir: Optional[str] = None,
                 logger: Optional[Logger] = None,
                 ) -> None:
        """A task in a workflow."""
        self.logger: Logger = logging.getLogger(
            __name__) if logger is None else logger
        self.name: str = name
        self.task_id: str = task_id
        self.runtime: float = runtime
        self.cores: Optional[float] = cores
        self.category: Optional[str] = category
        self.program: Optional[str] = program
        self.args: List[str] = args if args else []
        self.avg_cpu: Optional[float] = avg_cpu
        self.bytes_read: Optional[int] = bytes_read
        self.bytes_written: Optional[int] = bytes_written
        self.memory: Optional[int] = memory
        self.energy: Optional[int] = energy
        self.avg_power: Optional[float] = avg_power
        self.input_files: List[File] = input_files if input_files else []
        self.output_files: List[File] = output_files if output_files else []
        self.machines: Optional[List[Machine]] = machines
        self.priority: Optional[int] = priority
        self.type: Optional[TaskType] = task_type
        self.launch_dir: Optional[str] = launch_dir
        self.start_time: Optional[str] = str(datetime.now().astimezone().isoformat()) if not executedAt else executedAt
        self.logger.debug(
            f"created task {self.task_id}: runtime => {self.runtime} seconds.")

    def specification_as_dict(self) -> Dict:
        """A JSON representation of the task specification.

        :return: A JSON object representation of the task.
        :rtype: Dict
        """
        return {
            'name': self.name,
            'id': self.task_id,
            'parents': [],
            'children': [],
            'inputFiles': [f.file_id for f in self.input_files],
            'outputFiles': [f.file_id for f in self.output_files]
        }

    def execution_as_dict(self) -> Dict:
        """A JSON representation of the task execution.

        :return: A JSON object representation of the task.
        :rtype: Dict
        """
        task_obj = {
            'id': self.task_id,
            'runtimeInSeconds': self.runtime,
            'command': {}
        }
        if self.cores is not None:
            task_obj['coreCount'] = self.cores
        if self.avg_cpu is not None:
            task_obj['avgCPU'] = self.avg_cpu
        if self.bytes_read is not None:
            task_obj['readBytes'] = self.bytes_read
        if self.bytes_written is not None:
            task_obj['writtenBytes'] = self.bytes_written
        if self.memory is not None:
            task_obj['memoryInBytes'] = self.memory
        if self.energy is not None:
            task_obj['energyInKWh'] = self.energy
        if self.avg_power is not None:
            task_obj['avgPowerInKWh'] = self.avg_power
        if self.priority is not None:
            task_obj['priority'] = self.priority
        if self.program is not None:
            task_obj['command']['program'] = self.program
        if self.args is not None:
            task_obj['command']['arguments'] = self.args
        if self.machines:
            task_obj['machines'] = [m.name for m in self.machines]               
        if self.start_time:
            task_obj['executedAt'] = self.start_time
        if self.category is not None:
            task_obj['category'] = self.category
        if self.launch_dir:
            task_obj['launchDir'] = self.launch_dir
        return task_obj

    def __str__(self) -> str:
        return self.task_id
