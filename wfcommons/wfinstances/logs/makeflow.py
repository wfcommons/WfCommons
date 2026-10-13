#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import itertools
import math
import pathlib
import re

from datetime import datetime, timezone
from logging import Logger
from typing import Any, Dict, List, Optional

from .abstract_logs_parser import LogsParser
from ...common.file import File
from ...common.machine import Machine
from ...common.task import Task, TaskType
from ...common.workflow import Workflow
from ...utils import sanitize_wfformat_id

# Makeflow variable assignment (e.g., "OUTPUT=data/out.txt") and reference (e.g., "$(OUTPUT)" or "${OUTPUT}")
_VARIABLE_ASSIGNMENT = re.compile(r'^\s*(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)\s*=(.*)$')
_VARIABLE_REFERENCE = re.compile(r'\$(?:\(([A-Za-z_][A-Za-z0-9_]*)\)|\{([A-Za-z_][A-Za-z0-9_]*)\})')

# Resource Monitor units (cctools uses binary multiples)
_BYTES_UNITS = {'B': 1, 'kB': 1024, 'KB': 1024, 'KiB': 1024, 'MB': 1024 ** 2, 'MiB': 1024 ** 2,
                'GB': 1024 ** 3, 'GiB': 1024 ** 3, 'TB': 1024 ** 4, 'TiB': 1024 ** 4}
_SECONDS_UNITS = {'us': 1e-6, 'ms': 1e-3, 's': 1, 'm': 60, 'h': 3600}


class MakeflowLogsParser(LogsParser):
    """
    Parse Makeflow submit directory to generate workflow instance.

    :param execution_dir: Makeflow workflow execution directory (contains .mf/.makeflow and .makeflowlog files).
    :type execution_dir: pathlib.Path
    :param resource_monitor_logs_dir: Resource Monitor log files directory (created with `makeflow ----monitor=... ...`)
    :type resource_monitor_logs_dir: pathlib.Path
    :param description: Workflow instance description.
    :type description: Optional[str]
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Optional[Logger]
    """

    def __init__(self,
                 execution_dir: pathlib.Path,
                 resource_monitor_logs_dir: pathlib.Path,
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the makeflow log parser."""
        super().__init__('Makeflow', wms_url='https://ccl.cse.nd.edu/software/makeflow/', description=description,
                         logger=logger)

        # Sanity checks
        if not execution_dir.is_dir():
            raise OSError(f'The provided path does not exist or is not a folder: {execution_dir}')
        if not resource_monitor_logs_dir.is_dir():
            raise OSError(f'The provided path does not exist or is not a folder: {resource_monitor_logs_dir}')

        # Makeflow file
        files: List[pathlib.Path] = list(execution_dir.glob('*.mf'))
        if len(files) > 1:
            raise OSError(f'Multiple .mf files in: {execution_dir}')
        if len(files) == 0:
            files: List[pathlib.Path] = list(execution_dir.glob('*.makeflow'))
            if len(files) > 1:
                raise OSError(f'Multiple .makeflow files in: {execution_dir}')
            if len(files) == 0:
                raise OSError(f'Unable to find a .mf or .makeflow file in: {execution_dir}')
        self.mf_file: pathlib.Path = files[0]

        # Log file
        files = list(execution_dir.glob('*.makeflowlog'))
        if len(files) == 0:
            raise OSError(f'Unable to find .makeflowlog file in: {execution_dir}')
        if len(files) > 1:
            raise OSError(f'Multiple .makeflowlog files in: {execution_dir}')
        self.mf_log_file: pathlib.Path = files[0]
        if self.mf_log_file.read_text().count("# NODE") == 0:
            raise OSError(f'Not sufficiently verbose log file {self.mf_log_file}. Re-run the workflow with `makeflow --log-verbose ...`')

        self._execution_dir: pathlib.Path = execution_dir
        self._resource_monitor_logs_dir: pathlib.Path = resource_monitor_logs_dir
        self._files_map = {}
        self._args_map = {}
        self._variables: Dict[str, str] = {}
        self._machines: Dict[str, Machine] = {}

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow instance based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: Optional[str]

        :return: A workflow instance object.
        :rtype: Workflow
        """
        self.workflow_name = workflow_name

        # create base workflow instance object
        self.workflow = Workflow(name=self.workflow_name,
                                 description=self.description,
                                 runtime_system_name=self.wms_name,
                                 runtime_system_version=self.wms_version,
                                 runtime_system_url=self.wms_url)

        # parse workflow file
        self._parse_workflow_file()

        # parse makeflow log file
        self._parse_makeflow_log_file()

        # parse resource monitor files
        self._parse_resource_monitor_logs()

        return self.workflow

    def _parse_workflow_file(self) -> None:
        """Parse the makeflow workflow file and build the workflow structure."""
        task_id_counter = 1

        with open(self.mf_file) as f:
            outputs = []
            inputs = []
            for line in f:
                # print(f"Processing line: {line}")
                if line.lstrip().startswith('#'):
                    continue
                if '\t' not in line and _VARIABLE_ASSIGNMENT.match(line):
                    # variable assignment
                    match = _VARIABLE_ASSIGNMENT.match(line)
                    self._variables[match.group(1)] = self._expand_variables(match.group(2).strip().strip('"'))
                    continue
                if ':' in line and '\t' not in line:
                    line = self._expand_variables(line)
                    outputs = [self._get_file_name(file) for file in line.split(':')[0].split()]
                    inputs = [self._get_file_name(file) for file in line.split(':')[1].split()]

                    for file in itertools.chain(outputs, inputs):
                        if file not in self._files_map:
                            self._files_map[file] = {'task_name': None, 'children': [], 'file': []}

                elif '\t' in line:
                    # task execution command (likely olf here)
                    line = self._expand_variables(line)
                    prefix = line.replace('./', '').strip().split()[1 if 'LOCAL' in line else 0]
                    task_name = "{}_ID{:07d}".format(sanitize_wfformat_id(prefix), task_id_counter)

                    # create list of input and output files
                    output_files = self._create_files(outputs, "output", task_name)
                    input_files = self._create_files(inputs, "input", task_name)

                    # create task
                    args = ' '.join(line.split())
                    task = Task(name=task_name,
                                task_id=task_name,
                                category=prefix,
                                runtime=0,
                                program=prefix,
                                args=args.split(),
                                cores=1,
                                input_files=input_files,
                                output_files=output_files,
                                logger=self.logger)
                    self.workflow.add_task(task)
                    args = args.replace('\\\\', '\\')
                    self._args_map[args] = task
                    task_id_counter += 1

        # adding edges
        for file in self._files_map:
            for child in self._files_map[file]['children']:
                parent = self._files_map[file]['task_name']
                if parent and parent != child:
                    self.workflow.add_dependency(parent, child)

    def _expand_variables(self, line: str) -> str:
        """
        Expand the references to the Makeflow variables defined so far (e.g., $(VAR) or ${VAR}).

        :param line: A line of the makeflow file.
        :type line: str

        :return: The line with the known variables expanded.
        :rtype: str
        """
        def _replace(match: re.Match) -> str:
            name = match.group(1) or match.group(2)
            return self._variables.get(name, match.group(0))

        return _VARIABLE_REFERENCE.sub(_replace, line)

    @staticmethod
    def _get_file_name(file: str) -> str:
        """
        Get the workflow-level name of a file from a makeflow rule, in which a file can be
        renamed in the task sandbox (e.g., "local_name->remote_name").

        :param file: The file as written in a makeflow rule.
        :type file: str

        :return: The workflow-level file name.
        :rtype: str
        """
        return file.split('->')[0] if '->' in file and not file.startswith('->') else file

    def _create_files(self, files_list: List[str], input_or_output: str, task_name: str) -> List[File]:
        """
        Create a list of files objects.

        :param files_list: list of file names.
        :rtype files_list: List[str]
        :param input_or_output: Whether the files in the list are input or output.
        :rtype link: stc
        :param task_name: Task name.
        :rtype task_name: str

        :return: List of file objects.
        :rtype: List[File]
        """
        list_files = []
        for file in files_list:
            if self._files_map[file]['file']:
                list_files.append(
                    self._files_map[file]['file'][0] if input_or_output == "input" else self._files_map[file]['file'][1])
            else:
                size = 0
                file_path = self._execution_dir.joinpath(file)
                if file_path.is_dir():
                    size = sum(f.stat().st_size for f in file_path.glob("*") if f.is_file())
                elif file_path.is_file():
                    size = int(file_path.stat().st_size)

                file_obj_in = File(file_id=self._make_file_id(file),
                                   size=size,
                                   logger=self.logger)
                file_obj_out = File(file_id=self._make_file_id(file),
                                    size=size,
                                    logger=self.logger)
                list_files.append(file_obj_in if input_or_output == "input" else file_obj_out)
                self._files_map[file]['file'].extend([file_obj_in, file_obj_out])

            # files dependencies
            if input_or_output == "input":
                self._files_map[file]['children'].append(task_name)
            else:
                self._files_map[file]['task_name'] = task_name

        return list_files

    def _parse_makeflow_log_file(self):
        """Parse the makeflow log file and update workflow task information."""
        with open(self.mf_log_file) as f:
            start_time = 0

            for line in f:
                if 'STARTED' in line:
                    start_time = int(line.split()[2])
                    self.workflow.executed_at = datetime.fromtimestamp(start_time / 1000000, tz=timezone.utc).strftime(
                        '%Y-%m-%dT%H:%M:%S+00:00')


                elif 'COMPLETED' in line:
                    self.workflow.makespan = float('%.2f' % ((int(line.split()[2]) - start_time) / 1000000))

                elif line.startswith('# FILE') and 'condorlog' not in line:
                    file_name = self._get_file_name(line.split()[3])
                    if file_name in self._files_map:
                        size = int(line.split()[5])
                        for file_obj in self._files_map[file_name]['file']:
                            file_obj.size = size

    def _parse_resource_monitor_logs(self):
        """Parse the log files produced by resource monitor"""
        for file in self._resource_monitor_logs_dir.glob("*.summary"):
            with open(file) as f:
                data = json.load(f)

                # task
                command = data['command'].replace('perl', '').strip()
                if command not in self._args_map:
                    self.logger.warning(f'Unable to find the task of the resource monitor summary {file.name}: '
                                        f'{command}')
                    continue
                task = self._args_map[command]
                task.runtime = self._to_seconds(data['wall_time'])
                if 'start' in data:
                    start_time = datetime.fromtimestamp(self._to_seconds(data['start']), tz=timezone.utc)
                    task.start_time = start_time.isoformat()
                if 'cores' in data:
                    # a task runs on at least one core (the peak number of cores may be rounded down to 0)
                    task.cores = max(1, data['cores'][0] if isinstance(data['cores'], list) else data['cores'])
                # resource monitor reports memory and I/O in MB
                task.memory = self._to_bytes(data['memory'])
                task.bytes_read = self._to_bytes(data['bytes_read'])
                task.bytes_written = self._to_bytes(data['bytes_written'])
                if task.runtime > 0:
                    task.avg_cpu = float('%.4f' % (self._to_seconds(data['cpu_time']) / task.runtime * 100))
                if data.get('host'):
                    if data['host'] not in self._machines:
                        cpu = {'coreCount': int(data['machine_cpus'][0])} if 'machine_cpus' in data else None
                        self._machines[data['host']] = Machine(name=data['host'], cpu=cpu, logger=self.logger)
                    task.machines = [self._machines[data['host']]]

                # workflow (resource monitor and makeflow are both part of cctools, e.g., "7.8.2 FINAL")
                if data.get('monitor_version'):
                    self.wms_version = str(data['monitor_version']).split()[0]
                    self.workflow.runtime_system_version = self.wms_version

    @staticmethod
    def _to_bytes(measurement: Any) -> int:
        """
        Convert a resource monitor measurement (e.g., [12, "MB"]) into bytes.

        :param measurement: A [value, unit] measurement (a value with no unit is assumed in MB).
        :type measurement: Any

        :return: The measurement in bytes.
        :rtype: int
        """
        value, unit = (measurement[0], measurement[1]) if isinstance(measurement, list) else (measurement, 'MB')
        return int(round(float(value) * _BYTES_UNITS.get(unit, _BYTES_UNITS['MB'])))

    @staticmethod
    def _to_seconds(measurement: Any) -> float:
        """
        Convert a resource monitor measurement (e.g., [1.25, "s"] or [1700000001000000, "us"]) into seconds.

        :param measurement: A [value, unit] measurement (a value with no unit is assumed in seconds).
        :type measurement: Any

        :return: The measurement in seconds.
        :rtype: float
        """
        value, unit = (measurement[0], measurement[1]) if isinstance(measurement, list) else (measurement, 's')
        return float(value) * _SECONDS_UNITS.get(unit, 1)
