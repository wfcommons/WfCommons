#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import glob
import json
import itertools
import math
import os

from datetime import datetime
from logging import Logger
from typing import List, Optional

from .abstract_logs_parser import LogsParser
from ...common.file import File, FileLink
from ...common.machine import Machine
from ...common.task import Task, TaskType
from ...common.workflow import Workflow


class MakeflowLogsParser(LogsParser):
    """
    Parse Makeflow submit directory to generate workflow instance.

    :param execution_dir: Makeflow workflow execution directory (contains .mf and .makeflowlog files).
    :type execution_dir: str
    :param resource_monitor_logs_dir: Resource Monitor log files directory.
    :type resource_monitor_logs_dir: str
    :param description: Workflow instance description.
    :type description: str
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 execution_dir: str,
                 resource_monitor_logs_dir: str,
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the makeflow log parser."""
        super().__init__('Makeflow', 'http://ccl.cse.nd.edu/software/makeflow/', description, logger)

        # Sanity check
        if not os.path.isdir(execution_dir):
            raise OSError('The provided path does not exist or is not a folder: {}'.format(execution_dir))

        files = glob.glob('{}/*.mf'.format(execution_dir))
        if len(files) == 0:
            raise OSError('Unable to find .mf file in: {}'.format(execution_dir))
        self.mf_file = files[0]

        files = glob.glob('{}/*.makeflowlog'.format(execution_dir))
        if len(files) == 0:
            raise OSError('Unable to find .makeflowlog file in: {}'.format(execution_dir))
        self.mf_log_file = files[0]

        if not os.path.isdir(resource_monitor_logs_dir):
            raise OSError('The provided path does not exist or is not a folder: {}'.format(resource_monitor_logs_dir))

        self.execution_dir = execution_dir

        # self.mf_log_file = mf_log_file
        self.resource_monitor_logs_dir = resource_monitor_logs_dir
        self.files_map = {}
        self.args_map = {}

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow instance based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: str

        :return: A workflow instance object.
        :rtype: Workflow
        """
        self.workflow_name = workflow_name

        # create base workflow instance object
        self.workflow = Workflow(name=self.workflow_name,
                                 description=self.description,
                                 wms_name=self.wms_name,
                                 wms_url=self.wms_url)

        # parse workflow file
        self._parse_workflow_file()

        # parse makeflow log file
        self._parse_makeflow_log_file()

        # parse resource monitor files
        self._parse_resource_monitor_logs()

        return self.workflow

    def _parse_workflow_file(self):
        """Parse the makeflow workflow file and build the workflow structure."""
        task_id_counter = 1

        with open(self.mf_file) as f:
            outputs = []
            inputs = []
            for line in f:
                if ':' in line:
                    outputs = line.split(':')[0].split()
                    inputs = line.split(':')[1].split()

                    for file in itertools.chain(outputs, inputs):
                        if file not in self.files_map:
                            self.files_map[file] = {'task_name': None, 'children': [], 'file': []}

                elif len(line.strip()) > 0:
                    # task execution command
                    prefix = line.replace('./', '').replace('perl', '').strip().split()[1 if 'LOCAL' in line else 0]
                    task_name = "{}_ID{:07d}".format(prefix, task_id_counter)

                    # create list of task files
                    list_files = []
                    list_files.extend(self._create_files(outputs, FileLink.OUTPUT, task_name))
                    list_files.extend(self._create_files(inputs, FileLink.INPUT, task_name))

                    # create task
                    args = ' '.join(line.replace('LOCAL', '').replace('perl', '').strip().split())
                    task = Task(name=task_name,
                                task_id="ID{:07d}".format(task_id_counter),
                                category=prefix,
                                task_type=TaskType.COMPUTE,
                                runtime=0,
                                program=prefix,
                                args=args.split(),
                                cores=1,
                                files=list_files,
                                logger=self.logger)
                    self.workflow.add_node(task_name, task=task)
                    self.args_map[args] = task
                    task_id_counter += 1

        # adding edges
        for file in self.files_map:
            for child in self.files_map[file]['children']:
                if self.files_map[file]['task_name']:
                    self.workflow.add_edge(self.files_map[file]['task_name'], child)

    def _create_files(self, files_list: List[str], link: FileLink, task_name: str):
        """ Create a list of files objects.

        :param files_list: list of file names.
        :rtype files_list: List[str]
        :param link: Link type for the files in the list.
        :rtype link: FileLink
        :param task_name: Task name.
        :rtype task_name: str

        :return: List of file objects.
        :rtype: List[File]
        """
        list_files = []
        for file in files_list:
            if self.files_map[file]['file']:
                list_files.append(
                    self.files_map[file]['file'][0] if link == FileLink.INPUT else self.files_map[file]['file'][1])
            else:
                size = 0
                if os.path.isdir('{}/{}'.format(self.execution_dir, file)):
                    size = sum(os.path.getsize(f) for f in os.listdir('.') if os.path.isfile(f))
                elif os.path.isfile('{}/{}'.format(self.execution_dir, file)):
                    size = int(math.ceil(os.stat('{}/{}'.format(self.execution_dir, file)).st_size / 1000))  # B to KB

                file_obj_in = File(name=file,
                                   size=size,
                                   link=FileLink.INPUT,
                                   logger=self.logger)
                file_obj_out = File(name=file,
                                    size=size,
                                    link=FileLink.OUTPUT,
                                    logger=self.logger)
                list_files.append(file_obj_in if link == FileLink.INPUT else file_obj_out)
                self.files_map[file]['file'].extend([file_obj_in, file_obj_out])

            # files dependencies
            if link == FileLink.INPUT:
                self.files_map[file]['children'].append(task_name)
            else:
                self.files_map[file]['task_name'] = task_name

        return list_files

    def _parse_makeflow_log_file(self):
        """Parse the makeflow log file and update workflow task information."""
        with open(self.mf_log_file) as f:
            start_time = 0

            for line in f:
                if 'STARTED' in line:
                    start_time = int(line.split()[2])
                    self.workflow.executed_at = datetime.utcfromtimestamp(start_time / 1000000).strftime(
                        '%Y-%m-%dT%H:%M:%S+00:00')

                elif 'COMPLETED' in line:
                    self.workflow.makespan = float('%.2f' % ((int(line.split()[2]) - start_time) / 1000000))

                elif line.startswith('# FILE') and not 'condorlog' in line:
                    file_name = line.split()[3]
                    if file_name in self.files_map:
                        size = int(math.ceil(int(line.split()[5]) / 1000))  # B to KB
                        for file_obj in self.files_map[file_name]['file']:
                            file_obj.size = size

    def _parse_resource_monitor_logs(self):
        """Parse the log files produced by resource monitor"""
        for file in glob.glob('{}/*.summary'.format(self.resource_monitor_logs_dir)):
            with open(file) as f:
                data = json.load(f)

                # task
                task = self.args_map[data['command'].replace('perl', '').strip()]
                task.runtime = float(data['wall_time'][0])
                task.cores = float(data['cores'][0])
                task.memory = int(data['memory'][0]) * 1000  # MB to KB
                task.bytes_read = int(data['bytes_read'][0] * 1000)  # MB to KB
                task.bytes_written = int(data['bytes_written'][0] * 1000)  # MB to KB
                task.avg_cpu = float('%.4f' % (float(data['cpu_time'][0]) / float(data['wall_time'][0]) * 100))
                task.machine = Machine(name=data['host'],
                                       cpu={'count': int(data['machine_cpus'][0]), 'speed': 0, 'vendor': ''},
                                       logger=self.logger)

                # workflow
                self.workflow.wms_version = data['monitor_version']
