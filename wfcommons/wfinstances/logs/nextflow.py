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
import os

from logging import Logger
from typing import Optional

from .abstract_logs_parser import LogsParser
from ...common.file import File, FileLink
from ...common.machine import Machine, MachineSystem
from ...common.task import Task, TaskType
from ...common.workflow import Workflow


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
        # find execution report file
        files = glob.glob('{}/execution_report_*.html'.format(self.execution_dir))
        if len(files) == 0:
            raise OSError('Unable to find execution_report_*.html file in: {}'.format(self.execution_dir))
        execution_report_file = files[0]
        trace_data = None

        # parsing execution report file
        with open(execution_report_file) as f:
            read_trace_data = False
            for line in f:
                if 'Nextflow report data' in line:
                    read_trace_data = True
                    continue

                if not read_trace_data:
                    continue

                if 'window.data =' in line:
                    trace_data = line.replace('window.data = ', '').strip()
                else:
                    trace_data += line.replace('\\\\', '').replace('\\/', '/').replace('\\\'', '').replace(';', '')
                    read_trace_data = False

        trace_data = json.loads(trace_data)
        for t in trace_data['trace']:
            task_id = "ID{:06d}".format(int(t['task_id']))
            category = t['process'].lower().split(' ')[0]
            category = category[category.rfind(':') + 1:]
            task_name = '{}_{}'.format(category, task_id)
            task = Task(name=task_name,
                        task_id=task_id,
                        category=category,
                        task_type=TaskType.COMPUTE,
                        runtime=float(t['duration']) / 1000,
                        program=category,
                        args=list(filter(None, t['script'].replace('\n', '').split(' '))),
                        cores=float(t['cpus']),
                        files=[],
                        avg_cpu=float(t['%cpu'].replace('-', '0')),
                        bytes_read=round((int(t['rchar']) + int(t['read_bytes'])) / 1024),
                        bytes_written=round((int(t['wchar']) + int(t['write_bytes'])) / 1024),
                        memory=round(int(t['rss']) / 1024),
                        logger=self.logger)
            self.workflow.add_node(task_name, task=task)
