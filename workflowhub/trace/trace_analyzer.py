#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from typing import List
from .trace import Trace
from ..types import OutputFormat


class TraceAnalyzer:
	def __init__(self):
		self.traces: List[Trace] = []

	def append_trace(self, trace: Trace):
		if trace not in self.traces:
			self.traces.append(trace)

	def build_summary(self, output_format: OutputFormat = OutputFormat.CSV, normalized: bool = True):
		print(output_format)
