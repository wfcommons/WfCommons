#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pytest

from wfcommons.common import Task, TaskType, File, FileLink, Machine, MachineSystem


task_name = "Task Test"
task_id = "task_test_1"
input_files = [File(file_id="file_in_1", size=10, link=FileLink.INPUT), File(file_id="file_in_2", size=20, link=FileLink.INPUT)]
output_files = [File(file_id="file_out_1", size=30, link=FileLink.OUTPUT), File(file_id="file_out_2", size=40, link=FileLink.OUTPUT)]


class TestTask:
   
    @pytest.fixture
    def task(self) -> Task:
        return Task(
            name=task_name,
            task_id=task_id,
            runtime=123.45,
            input_files=input_files,
            output_files=output_files,
            category="task_test",
            machines=[Machine(name="machine_1", cpu = {"coreCount": 48, "speedInMHz": 1200, "vendor": "Vendor Name"})],
            program="program",
            args=["arg_1", "arg_2"],
            avg_cpu=0.5,
            bytes_read=12345,
            bytes_written=54321,
            memory=10,
            energy=100,
            avg_power=1.0,
            priority=100,
            executedAt="2024-09-15T08:59:33.699321-04:00",
            task_type=TaskType.COMPUTE,
            launch_dir="/tmp",
        )

    @pytest.mark.unit
    def test_task_specification(self, task: Task) -> None:

        task_specification = {
            "name": task_name,
            "id": task_id,
            "parents": [],
            "children": [],
            "inputFiles": [str(f) for f in input_files],
            "outputFiles": [str(f) for f in output_files]
        }

        assert(task_specification == task.specification_as_dict())
