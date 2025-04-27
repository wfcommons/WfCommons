#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import pathlib
import pytest

from wfcommons.wfinstances import Instance, InstanceAnalyzer


class TestInstanceAnalyzer:

    @pytest.fixture
    def instance(self) -> pathlib.Path:
        instance = {
            "name": "workflow_test",
            "description": "Instance generate for WfCommons Test",
            "createdAt": "2020-12-30T02:19:01.238077",
            "schemaVersion": "1.5",
            "author": {
                "name": "wfcommons",
                "email": "support@wfcommons.org"
            },
            "workflow": {
                "specification": {
                    "tasks": [
                        {
                            "name": "task_01",
                            "id": "task_01",
                            "children": ["task_02", "task_03"],
                            "inputFiles": ["input_01.txt", "input_02.txt"],
                            "outputFiles": ["output_01.txt"],
                            "parents": []
                        }
                    ],
                    "files": [
                        {"id": "input_01.txt", "sizeInBytes": 10}, {"id": "input_02.txt", "sizeInBytes": 20},
                        {"id": "output_01.txt", "sizeInBytes": 15}
                    ]
                },
                "execution": {
                    "makespanInSeconds": 111.5,
                    "executedAt": "2024-11-30T03:25:55+00:00",
                    "tasks": [
                        {
                            "id": "task_01",
                            "runtimeInSeconds": 0.052203,
                            "command": {
                                "program": "program_1",
                                "arguments": ["a", "b"]
                            },
                            "coreCount": 1,
                            "avgCPU": 90.123,
                            "readBytes": 150,
                            "writtenBytes": 200,
                            "memoryInBytes": 300,
                            "machines": ["worker-1"]
                        }
                    ],
                    "machines": [
                        {"nodeName": "worker-1", "cpu": {"coreCount": 24}}
                    ]
                }
            },
            "runtimeSystem": {"name": "wfcommons", "version": "1.0", "url": "https://wfcommons.org"}
        }
        instance_path = pathlib.Path("/tmp/wfcommons_test_workflow.json")
        instance_path.write_text(json.dumps(instance))
        return instance_path

    @pytest.mark.unit
    def test_instance_analyzer_summary(self, instance : pathlib.Path) -> None:
        instance = Instance(instance)
        analyzer = InstanceAnalyzer()
        analyzer.append_instance(instance)
        workflow_tasks = ['task_01']
        summary = analyzer.build_summary(workflow_tasks, include_raw_data=False)
        assert(summary == {'task_01': {'runtime': {'min': 0.052203, 'max': 0.052203, 'distribution': {'name': 'dweibull', 'params': [0.19451622991431783, 4.155260197271249e-34, 1.161132255170155]}}, 'input': {'.txt': {'distribution': {'name': 'rdist', 'params': [1.5504806356651624, 0.0013236200991527764, 0.0013236200991527767]}, 'min': 10, 'max': 20}}, 'output': {'.txt': {'distribution': None, 'min': 15, 'max': 15}}}})
