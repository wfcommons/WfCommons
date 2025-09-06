#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import pytest
import requests
import json

from datetime import datetime
from wfcommons.common import Task, Workflow
from wfcommons.wfinstances import Instance
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
            makespan=100.0,
            executed_at=str(datetime.now().astimezone().isoformat())
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
                "url": f"https://docs.wfcommons.org/en/v{__version__}/"
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

    @pytest.mark.unit
    def test_workflow_json_generation(self):

        # Put a JSON file in /tmp
        url = "https://raw.githubusercontent.com/wfcommons/WfInstances/refs/heads/main/makeflow/blast/blast-chameleon-small-001.json"
        response = requests.get(url)
        local_file_name = url.split("/")[-1]
        with open("/tmp/" + local_file_name, 'wb') as f:
            f.write(response.content)

        # Create an instance from the JSON File and write it back to a JSON
        instance = Instance(pathlib.Path("/tmp") / local_file_name)

        # Testing the "iterator" capability
        assert(len(list(instance)) == len(instance.instance["workflow"]["specification"]["tasks"]))

        # Test writing to JSON
        instance.workflow.write_json(pathlib.Path("/tmp/written_workflow.json"))

        # Get the two jsons as objects
        with open("/tmp/" + local_file_name) as f1, open("/tmp/written_workflow.json") as f2:
            original_json = json.load(f1)
            written_json = json.load(f2)

        # Fix things that will be always (rightly) different in the written_json
        written_json["description"] = original_json["description"]
        written_json["createdAt"] = original_json["createdAt"]
        written_json["author"] = original_json["author"]
        written_json["workflow"]["execution"]["executedAt"] = original_json["workflow"]["execution"]["executedAt"]

        original_json["workflow"]["specification"]["files"] = sorted(original_json["workflow"]["specification"]["files"], key=lambda x: x['id'])
        written_json["workflow"]["specification"]["files"] = sorted(written_json["workflow"]["specification"]["files"], key=lambda x: x['id'])

        # Compare the two jsons!
        assert(original_json == written_json)

    @pytest.mark.unit
    def test_workflow_dot_file(self):

        # Put a JSON file in /tmp
        url = "https://raw.githubusercontent.com/wfcommons/WfInstances/refs/heads/main/makeflow/blast/blast-chameleon-small-001.json"
        response = requests.get(url)
        local_file_name = url.split("/")[-1]
        with open("/tmp/" + local_file_name, 'wb') as f:
            f.write(response.content)

        # Create an instance from the JSON File and write it back to a JSON
        instance = Instance(pathlib.Path("/tmp") / local_file_name)

        # Capture some metrics
        num_tasks = len(instance.workflow.tasks)
        num_dependencies = len(instance.workflow.edges)

        # # Create a dot file
        dot_path = pathlib.Path("/tmp/written_workflow.dot")
        instance.workflow.write_dot(dot_path)
        assert dot_path.exists()
        with open(str(dot_path), "r", encoding="utf-8") as f:
            content = f.read()
            assert(num_tasks == content.count("label") - 1)  # Extra "label" in file for \N
            assert(num_dependencies == content.count("->"))  # Extra "label" in file for \N

        # Read it back
        instance.workflow.read_dot(dot_path)
        assert(num_tasks == len(instance.workflow.tasks))
        assert(num_tasks == len(instance.workflow.nodes))
        assert(num_dependencies == len(instance.workflow.edges))




