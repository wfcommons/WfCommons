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
import os

from logging import Logger
from typing import Optional

from .abstract_logs_parser import LogsParser
from ...common.file import File, FileLink
from ...common.machine import Machine, MachineSystem
from ...common.task import Task, TaskType
from ...common.workflow import Workflow


def _storage_unit_conversion(value: str) -> int:
    if '-' in value:
        return 0
    if 'KB' in value:
        return round(float(value.replace('KB', '')))
    if 'MB' in value:
        return round(float(value.replace('MB', '')) * 1000)


class NextflowLogsParser(LogsParser):
    """
    Parse Nextflow submit directory to generate workflow trace.

    :param execution_dir: Nextflow's execution directory.
    :type execution_dir: str
    :param description: Workflow instance description.
    :type description: str
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 execution_dir: str,
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the nextflow log parser."""
        super().__init__('Nextflow', 'https://www.nextflow.io', description, logger)

        # Sanity check
        if not os.path.isdir(execution_dir):
            raise OSError('The provided path does not exist or is not a directory: {}'.format(execution_dir))

        self.execution_dir = execution_dir
        self.files_map = {}
        self.text_files = None
        self.line_count = None

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow trace based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: str

        :return: A workflow trace object.
        :rtype: Workflow
        """
        self.workflow_name = workflow_name
        self.workflow = Workflow(name=self.workflow_name,
                                 description=self.description,
                                 executed_at=self.executed_at,
                                 makespan=self.makespan)

        self._parse_workflow_file()

        return self.workflow

    def _parse_workflow_file(self):
        """Parse the Nextflow workflow file and build the workflow structure."""
        # find execution trace file
        files = glob.glob('{}/execution_trace_*.txt'.format(self.execution_dir))
        if len(files) == 0:
            raise OSError('Unable to find execution_trace_*.txt file in: {}'.format(self.execution_dir))
        execution_trace_file = files[0]

        with open(execution_trace_file) as f:
            for line in f:
                if line.startswith('task_id'):
                    continue

                contents = line.strip().split('\t')
                task_id = "ID{:06d}".format(int(contents[0]))
                category = contents[3].lower().split(' ')[0]
                category = category[category.rfind(':') + 1:]
                task_name = '{}_{}'.format(category, task_id)
                duration = contents[7].split(' ')
                runtime = 0
                for d in duration:
                    if 'ms' in d:
                        runtime += float(d.replace('ms', '')) / 100
                    elif 's' in d:
                        runtime += float(d.replace('s', ''))
                    elif 'm' in d:
                        runtime += float(d.replace('m', '')) * 60

                task = Task(name=task_name,
                            task_id=task_id,
                            category=category,
                            task_type=TaskType.COMPUTE,
                            runtime=runtime,
                            program=category,
                            args=[],
                            cores=1,
                            files=[],
                            avg_cpu=float(contents[9].replace('-', '0').replace('%', '')),
                            bytes_read=_storage_unit_conversion(contents[12]),
                            bytes_written=_storage_unit_conversion(contents[13]),
                            memory=_storage_unit_conversion(contents[10]),
                            logger=self.logger)
                self.workflow.add_node(task_name, task=task)
                print(task.as_dict())
