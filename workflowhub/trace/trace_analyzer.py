#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import re

from typing import List, Tuple
from .trace import Trace
from ..common.job import Job
from ..types import OutputFormat, JsonDict
from ..utils import best_fit_distribution


class TraceAnalyzer:
	def __init__(self):
		self.traces: List[Trace] = []
		self.jobs_summary: JsonDict[str, List[Tuple[float, int]]] = {}

	def append_trace(self, trace: Trace):
		if trace not in self.traces:
			self.traces.append(trace)

	def build_summary(self, output_format: OutputFormat = OutputFormat.CSV):
		# build jobs summary
		for trace in self.traces:
			for node in trace.workflow.nodes.data():
				job: Job = node[1]['job']
				job_name: str = re.sub(r'_ID[0-9]*', '', job.name)
				input_size: int = 0
				for file in job.files:
					if file.link == "input":
						input_size += file.size
				if job_name not in self.jobs_summary:
					self.jobs_summary[job_name] = [(job.runtime, input_size)]
				else:
					self.jobs_summary[job_name].append((job.runtime, input_size))

		# build traces summary
		traces_summary: JsonDict[str, Tuple[float, int]]

		for job in self.jobs_summary:
			runtime_list: List[float] = []
			inputs_list: List[int] = []
			for entry in self.jobs_summary[job]:
				runtime_list.append(entry[0])
				inputs_list.append(entry[1])

			runtime_best_distribution = best_fit_distribution(runtime_list)
			inputs_best_distribution = best_fit_distribution(inputs_list)
			print(job)
			print(runtime_best_distribution)
			print(inputs_best_distribution)
