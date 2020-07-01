#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from workflowhub.trace import Trace

# TODO: If the schema is not available, download it from the git
SCHEMA_TRACES = "data/schemes/traces.json"

def main(trace_path: str) -> None:
	trace = Trace(input_trace = trace_path, schema = SCHEMA_TRACES)
	# print(trace.schema.keys())
	print(trace.makespan)
	# for e in trace.workflow.edges:
	# 	print(e[0], trace.workflow.get_edge_data(*e), list(trace.workflow.adj[e[0]]))

	# Print job following a topological order !
	for job in trace:
		print(job)
	trace.write_dot()

if __name__ == '__main__':
	main("test.json")
