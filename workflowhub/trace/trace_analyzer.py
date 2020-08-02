#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from os import path
from typing import Any, Dict, List, Optional, Tuple
from .trace import Trace
from ..common.job import Job
from ..common.file import FileLink
from ..utils import best_fit_distribution


class TraceAnalyzer:
    def __init__(self):
        self.traces: List[Trace] = []
        self.jobs_summary: Dict[str, List:Job] = {}
        self.traces_summary: Dict[str, Dict[str, Any]] = {}

    def append_trace(self, trace: Trace):
        if trace not in self.traces:
            self.traces.append(trace)

    def build_summary(self, jobs_list: List[str], include_raw_data: Optional[bool] = True):
        """
        :param jobs_list:
        :param include_raw_data:
        """
        # build jobs summary
        for trace in self.traces:
            for node in trace.workflow.nodes.data():
                job: Job = node[1]['job']
                job_name: str = [j for j in jobs_list if j in job.name][0]

                if job_name not in self.jobs_summary:
                    self.jobs_summary[job_name] = []
                self.jobs_summary[job_name].append(job)

        # build traces summary
        for job_name in self.jobs_summary:
            runtime_list: List[float] = []
            inputs_dict: Dict[str, Any] = {}
            outputs_dict: Dict[str, Any] = {}

            for job in self.jobs_summary[job_name]:
                runtime_list.append(job.runtime)

                for file in job.files:
                    extension: str = path.splitext(file.name)[1] if '.' in file.name else file.name
                    if file.link == FileLink.INPUT:
                        self._append_file_to_dict(extension, inputs_dict, file.size)
                    elif file.link == FileLink.OUTPUT:
                        self._append_file_to_dict(extension, outputs_dict, file.size)

            self._best_fit_distribution_for_file(inputs_dict, include_raw_data)
            self._best_fit_distribution_for_file(outputs_dict, include_raw_data)

            self.traces_summary[job_name] = {
                'runtime': {
                    'min': min(runtime_list),
                    'max': max(runtime_list),
                    'distribution': self._json_format_distribution_fit(best_fit_distribution(runtime_list))
                },
                'input': inputs_dict,
                'output': outputs_dict
            }
            if include_raw_data:
                self.traces_summary[job_name]['runtime']['data'] = runtime_list

        return self.traces_summary

    def _append_file_to_dict(self, extension, dict_obj, file_size):
        """
        :param extension:
        :param dict_obj:
        :param file_size:
        """
        if extension not in dict_obj:
            dict_obj[extension] = {'data': [], 'distribution': None}
        dict_obj[extension]['data'].append(file_size)

    def _best_fit_distribution_for_file(self, dict_obj, include_raw_data):
        """
        :param dict_obj:
        :param include_raw_data:
        """
        for ext in dict_obj:
            dict_obj[ext]['min'] = min(dict_obj[ext]['data'])
            dict_obj[ext]['max'] = max(dict_obj[ext]['data'])
            if dict_obj[ext]['min'] != dict_obj[ext]['max']:
                dict_obj[ext]['distribution'] = self._json_format_distribution_fit(
                    best_fit_distribution(dict_obj[ext]['data']))
            if not include_raw_data:
                del dict_obj[ext]['data']

    def _json_format_distribution_fit(self, dist_tuple: Tuple):
        """
        :param dist_tuple:
        """
        formatted_entry = {'name': dist_tuple[0], 'params': []}
        for p in dist_tuple[1]:
            formatted_entry['params'].append(p)
        return formatted_entry
