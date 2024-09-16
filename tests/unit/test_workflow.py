#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import pytest

from wfcommons.common import Task, Workflow
from wfcommons.version import __version__, __schema_version__


tasks_list = [
    ([Task(name="task_1", task_id="task_1", runtime=15.0)]),
    ([Task(name="task_1", task_id="task_1", runtime=15.0), Task(name="task_2", task_id="task_2", runtime=30.0)]),
    ([Task(name="task_1", task_id="task_1", runtime=15.0), Task(name="task_2", task_id="task_2", runtime=30.0), Task(name="task_3", task_id="task_3", runtime=60.0)]),
]


class TestWorkflow:
   
    @pytest.fixture
    def workflow(self) -> Workflow:
        return Workflow(
            name="Workflow Test",
            makespan=100.0
        )

    @pytest.mark.unit
    def test_workflow_creation(self, workflow: Workflow) -> None:
        workflow.write_json(pathlib.Path("/tmp/workflow_test.json"))

        workflow_json = {
            "name": "Workflow Test",
            "description": "Instance generated with WfCommons - https://wfcommons.org",
            "createdAt": workflow.created_at,
            "schemaVersion": f"{__schema_version__}",
            "author": {
                "name": f"{workflow.author_name}",
                "email": "support@wfcommons.org"
            },
            "workflow": {
                "specification": {
                    "tasks": [],
                    "files": []
                },
                "execution": {
                    "makespanInSeconds": 100.0,
                    "executedAt": workflow.executed_at,
                    "tasks": []
                }
            },
            "runtimeSystem": {
                "name": "WfCommons",
                "version": f"{__version__}",
                "url": "https://docs.wfcommons.org/en/v1.1.dev/"
            }
        }

        assert(workflow.workflow_json == workflow_json)

    @pytest.mark.unit
    def test_workflow_empty_tasks(self, workflow: Workflow) -> None:
        assert(not workflow.tasks)

    @pytest.mark.parametrize(("tasks"), tasks_list)
    @pytest.mark.unit
    def test_workflow_add_tasks(self, workflow: Workflow, tasks: list[Task]) -> None:
        tasks_list = []
        for task in tasks:
            workflow.add_task(task)
            tasks_list.append(task.task_id)
        assert(workflow.leaves() == tasks_list)

    @pytest.mark.parametrize(("tasks"), tasks_list)
    @pytest.mark.unit
    def test_workflow_add_dependencies_roots(self, workflow: Workflow, tasks: list[Task]) -> None:
        previous_task = None
        for task in tasks:
            workflow.add_task(task)
            if previous_task:
                workflow.add_dependency(previous_task.task_id, task.task_id)
            previous_task = task
        assert(workflow.roots() == [tasks[0].task_id])

    @pytest.mark.parametrize(("tasks"), tasks_list)
    @pytest.mark.unit
    def test_workflow_add_dependencies_leaves(self, workflow: Workflow, tasks: list[Task]) -> None:
        previous_task = None
        for task in tasks:
            workflow.add_task(task)
            if previous_task:
                workflow.add_dependency(previous_task.task_id, task.task_id)
            previous_task = task
        assert(workflow.leaves() == [previous_task.task_id])
