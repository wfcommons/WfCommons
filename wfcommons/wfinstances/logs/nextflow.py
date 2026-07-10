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
import pathlib
from datetime import datetime, timedelta
import csv


from logging import Logger
from typing import Dict, Optional

from .abstract_logs_parser import LogsParser
from ...common.task import Task, TaskType
from ...common.workflow import Workflow


class NextflowLogsParser(LogsParser):
    """"
    Parse Nextflow execution directory to generate a workflow trace. Note that the workflow
    reconstruction is not perfect, and will likely only capture file-based data dependencies.
    The workflow must have been executed with two features enabled: 1) the nf-prov plugin, and
    2) execution tracing. This can be achieved by invoking Nextflow with a config
    file, e.g.::

        nextflow run -c nextflow-wfcommons.config ...

    where ``nextflow-wfcommons.config`` contains:

    .. code-block:: groovy

        plugins {
            id 'nf-prov'
        }

        prov {
            enabled = true
            formats {
                wrroc {
                    file = 'ro-crate-metadata.json'
                    overwrite = true
                }
            }
        }

        trace {
            enabled = true
            file = "results/pipeline_info/execution_trace_${new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')}.txt"
            overwrite = true
        }

    :param execution_dir: Nextflow's execution directory.
    :type execution_dir: pathlib.Path
    :param nextflow_version: The Nextflow version used to execute the workflow
    :type nextflow_version: str
    :param trace_file_name_pattern: The trace file name pattern to find the trace file (default: "*trace*.txt")
    :type trace_file_name_pattern: str
    :param description: Workflow instance description.
    :type description: Optional[str]
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Optional[Logger]
    """

    def __init__(self,
                 execution_dir: pathlib.Path,
                 nextflow_version: str,
                 trace_file_name_pattern: Optional[str] = "*trace*.txt",
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the nextflow log parser."""
        super().__init__('Nextflow', wms_version=nextflow_version, wms_url='https://www.nextflow.io',
                         description=description, logger=logger)

        self.execution_dir = execution_dir

        # Load the Nextflow execution trace, and create a task runtime dictionary
        self.nextflow_execution_trace_files = list(self.execution_dir.rglob(trace_file_name_pattern))
        if len(self.nextflow_execution_trace_files) == 0:
            raise FileNotFoundError("No execution_trace_*.txt file found in Nextflow execution directory.")


    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow instance based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: Optional[str]

        :return: A workflow instance object.
        :rtype: Workflow
        """

        # Parse the Nextflow (most recent) execution trace file to create a dict of task runtimes
        nextflow_execution_trace_file: pathlib.Path = max(self.nextflow_execution_trace_files,
                                                          key=lambda p: p.stat().st_mtime)
        nextflow_task_runtimes = self._load_nextflow_trace(nextflow_execution_trace_file)


        # Create an RO-Create parser
        from wfcommons.wfinstances import ROCrateLogsParser
        ro_crate_parser = ROCrateLogsParser(self.execution_dir,
                                            wms_name="Nextflow",
                                            wms_version=self.wms_version,
                                            wms_url=self.wms_url,
                                            description=self.description, logger=self.logger,
                                            steps_to_ignore=None,
                                            file_extensions_to_ignore=None,
                                            instruments_to_ignore=None,
                                            task_runtimes=nextflow_task_runtimes)

        return ro_crate_parser.build_workflow(workflow_name)


    def _load_nextflow_trace(self, trace_file: pathlib.Path) -> dict[str, tuple[str, str]]:
        """
        Parse a Nextflow execution trace file.
        Returns dict of task_name -> (start_iso, end_iso).
        """
        times = {}
        with open(trace_file, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                name = row['name'].strip()
                submit_str = row['submit'].strip()
                duration_str = row['duration'].strip()
                if not submit_str or not duration_str:
                    continue
                try:
                    start = datetime.strptime(submit_str, '%Y-%m-%d %H:%M:%S.%f')
                    duration_sec = self._parse_duration(duration_str)
                    end = start + timedelta(seconds=duration_sec)
                    # Store as ISO strings to match _time_diff expectations
                    times[name] = (start.isoformat(), end.isoformat())
                except Exception:
                    continue
        return times

    def _parse_duration(self, duration_str: str) -> float:
        """Parse Nextflow duration strings like '2m 40s', '519ms', '1h 3m 2s'."""
        import re
        total = 0.0
        for value, unit in re.findall(r'(\d+(?:\.\d+)?)\s*(h|m|s|ms)', duration_str):
            value = float(value)
            if unit == 'h':
                total += value * 3600
            elif unit == 'm':
                total += value * 60
            elif unit == 's':
                total += value
            elif unit == 'ms':
                total += value / 1000
        return total

#     def __init__(self,
#                  execution_dir: pathlib.Path,
#                  description: Optional[str] = None,
#                  logger: Optional[Logger] = None) -> None:
#         """Create an object of the nextflow log parser."""
#         super().__init__('Nextflow', 'https://www.nextflow.io', description, logger)
#
#         # Sanity check
#         if not execution_dir.is_dir():
#             raise OSError(f'The provided path does not exist or is not a directory: {execution_dir}')
#
#         self.execution_dir = execution_dir
#         self.files_map = {}
#         self.tasknames_to_taskids = {}
#         self.text_files = None
#         self.line_count = None
#
#     def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
#         """
#         Create workflow trace based on the workflow execution logs.
#
#         :param workflow_name: The workflow name.
#         :type workflow_name: Optional[str]
#
#         :return: A workflow trace object.
#         :rtype: Workflow
#         """
#         self.workflow_name = workflow_name
#         self.workflow = Workflow(name=self.workflow_name,
#                                  description=self.description,
#                                  executed_at=self.executed_at,
#                                  makespan=self.makespan,
#                                  # wms_name=self.wms_name,
#                                  # wms_url=self.wms_url
#                                  )
#
#         self._parse_execution_report_file()
#         self._parse_execution_timeline_file()
#
#         return self.workflow
#
#     def _parse_execution_report_file(self) -> None:
#         """Parse the Nextflow execution report file and gather the tasks information."""
#         trace_data = self._read_data('execution_report*.html')
#
#         for t in trace_data['trace']:
#             task_id = "ID{:06d}".format(int(t['task_id']))
#             category = t['process'].lower().split(' ')[0]
#             category = _parse_task_name(category)
#             task_name = _parse_task_name(t['name'])
#             task = Task(name=task_name,
#                         task_id=task_id,
#                         category=category,
#                         task_type=TaskType.COMPUTE,
#                         runtime=float(_parse_number(t['duration'])) / 1000,
#                         program=category,
#                         args=list(filter(None, t['script'].replace('\n', '').split(' '))),
#                         cores=float(t['cpus']),
#                         #files=[],
#                         avg_cpu=float(_parse_number(t['%cpu'])),
#                         bytes_read=round((int(_parse_number(t['rchar'])) + int(_parse_number(t['read_bytes']))) / 1024),
#                         bytes_written=round(
#                             (int(_parse_number(t['wchar'])) + int(_parse_number(t['write_bytes']))) / 1024),
#                         memory=round(int(_parse_number(t['rss'])) / 1024),
#                         logger=self.logger)
#             self.workflow.add_task(task)
#             self.tasknames_to_taskids[task_name] = task_id
#
#     def _parse_execution_timeline_file(self) -> None:
#         """Parse the Nextflow execution timeline file and build the workflow structure."""
#         timeline_data = self._read_data('execution_timeline*.html')
#         tasks_map = {}
#         max_index = 0
#
#         for e in timeline_data['processes']:
#             task_name = _parse_task_name(e['label'])
#             index = int(e['index'])
#             max_index = max(index, max_index)
#             if index not in tasks_map:
#                 tasks_map[index] = []
#             tasks_map[index].append(task_name)
#
#         for index in range(max_index + 1):
#             if index > 0:
#                 for c in tasks_map[index]:
#                     for p in tasks_map[index - 1]:
#                         self.workflow.add_edge(self.tasknames_to_taskids[p], self.tasknames_to_taskids[c])
#
#         self.workflow.makespan = float(
#             (int(timeline_data['endingMillis']) - int(timeline_data['beginningMillis'])) / 1000)
#
#     def _read_data(self, file_format: str) -> Dict:
#         """
#         Read data into a JSON from a file that matches the format.
#
#         :param file_format: File format to be searched
#         :type file_format: str
#
#         :return: Data in JSON format
#         :rtype: Dict
#         """
#         files = glob.glob(f'{self.execution_dir}/{file_format}')
#         if len(files) == 0:
#             raise OSError(f'Unable to find {file_format} in {self.execution_dir}')
#
#         data = None
#
#         # parsing execution report file
#         with open(files[0]) as f:
#             self.logger.debug(f'Reading data from: {files[0]}')
#             read_trace_data = False
#             read_nextflow_version = False
#
#             for line in f:
#                 if 'Nextflow report data' in line:
#                     read_trace_data = True
#                     continue
#
#                 if 'Nextflow version' in line:
#                     read_nextflow_version = True
#                     continue
#
#                 if not read_trace_data and not read_nextflow_version:
#                     continue
#
#                 if read_nextflow_version:
#                     version = line.strip().split(' ')[2].replace(',', '')
#                     self.workflow.wms_version = version
#                     read_nextflow_version = False
#                     continue
#
#                 if 'window.data =' in line:
#                     data = line.replace('window.data = ', '').strip()
#                 elif line.startswith(';') or '</script>' in line:
#                     read_trace_data = False
#                 elif len(line) > 0:
#                     data += line.replace('\\\\', '').replace('\\/', '/').replace('\\\'', '').replace(';', '')
#
#         return json.loads(data)
#
#
# def _parse_task_name(task_name: str):
#     """
#     Format the task name.
#
#     :param task_name: Raw task name
#     :type task_name: str
#
#     :return: Formatted task name
#     :rtype: str
#     """
#     return task_name[task_name.rfind(':') + 1:].replace('(', '').replace(')', '').replace(' ', '_').lower()
#
#
# def _parse_number(number: str):
#     """
#     Format a number.
#
#     :param number: Raw number
#     :type number: str
#
#     :return: Formatted number
#     :rtype: str
#     """
#     return number.replace('-', '0')
