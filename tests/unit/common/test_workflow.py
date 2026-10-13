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
import math
from numbers import Real
from typing import Any

from datetime import datetime
from wfcommons.common import File, Task, Workflow
from wfcommons.wfinstances import Instance, SchemaValidator
from wfcommons.version import __version__, __schema_version__


tasks_list = [
    ([Task(name="task_1", task_id="task_1", runtime=15.0)]),
    ([Task(name="task_1", task_id="task_1", runtime=15.0), Task(name="task_2", task_id="task_2", runtime=30.0)]),
    ([Task(name="task_1", task_id="task_1", runtime=15.0), Task(name="task_2", task_id="task_2", runtime=30.0), Task(name="task_3", task_id="task_3", runtime=60.0)]),
]


def assert_json_close(
    actual: Any,
    expected: Any,
    *,
    rel_tol: float = 1e-9,
    abs_tol: float = 1e-12,
    path: str = "$",
) -> None:
    """Assert that two JSON-like objects are equal, with close numerics."""

    actual_is_number = (
        isinstance(actual, Real)
        and not isinstance(actual, bool)
    )
    expected_is_number = (
        isinstance(expected, Real)
        and not isinstance(expected, bool)
    )

    if actual_is_number or expected_is_number:
        assert actual_is_number and expected_is_number, (
            f"Numeric type mismatch at {path}: "
            f"actual={actual!r}, expected={expected!r}"
        )

        assert math.isclose(
            actual,
            expected,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
        ), (
            f"Numbers differ at {path}: "
            f"actual={actual!r}, expected={expected!r}, "
            f"rel_tol={rel_tol}, abs_tol={abs_tol}"
        )
        return

    if isinstance(actual, dict) or isinstance(expected, dict):
        assert isinstance(actual, dict) and isinstance(expected, dict), (
            f"Type mismatch at {path}: "
            f"actual={type(actual).__name__}, "
            f"expected={type(expected).__name__}"
        )

        assert actual.keys() == expected.keys(), (
            f"Dictionary keys differ at {path}: "
            f"actual-only={actual.keys() - expected.keys()}, "
            f"expected-only={expected.keys() - actual.keys()}"
        )

        for key in actual:
            assert_json_close(
                actual[key],
                expected[key],
                rel_tol=rel_tol,
                abs_tol=abs_tol,
                path=f"{path}.{key}",
            )
        return

    if isinstance(actual, list) or isinstance(expected, list):
        assert isinstance(actual, list) and isinstance(expected, list), (
            f"Type mismatch at {path}: "
            f"actual={type(actual).__name__}, "
            f"expected={type(expected).__name__}"
        )

        assert len(actual) == len(expected), (
            f"List lengths differ at {path}: "
            f"actual={len(actual)}, expected={len(expected)}"
        )

        for index, (actual_item, expected_item) in enumerate(
            zip(actual, expected)
        ):
            assert_json_close(
                actual_item,
                expected_item,
                rel_tol=rel_tol,
                abs_tol=abs_tol,
                path=f"{path}[{index}]",
            )
        return

    assert actual == expected, (
        f"Values differ at {path}: "
        f"actual={actual!r}, expected={expected!r}"
    )



class TestWorkflow:
   
    @pytest.fixture
    def workflow(self) -> Workflow:
        return Workflow(
            name="Workflow Test",
            makespan=100.0,
            executed_at=str(datetime.now().astimezone().isoformat()),
            runtime_system_name="WfCommons",
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
                    "files": [],
                    "metrics": {
                        "numberOfTasks": 0,
                        "numberOfFiles": 0
                    }
                },
                "execution": {
                    "makespanInSeconds": 100.0,
                    "executedAt": workflow.executed_at,
                    "tasks": []
                }
            },
            "runtimeSystem": {
                "name": "WfCommons",
                "version": f"unknown",
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
        # assert(original_json == written_json)
        # assert(original_json == pytest.approx(written_json))
        assert_json_close(original_json,written_json)

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


def _build_workflow(tasks: list[tuple[str, str]], dependencies: list[tuple[str, str]], **kwargs) -> Workflow:
    """Build a workflow from (name, id) pairs and (parent id, child id) pairs."""
    workflow = Workflow(name="wfformat_test", makespan=10.0, **kwargs)
    for name, task_id in tasks:
        workflow.add_task(Task(name=name, task_id=task_id, runtime=1.0,
                               output_files=[File(f"{task_id}_out", 10)]))
    for parent, child in dependencies:
        workflow.add_dependency(parent, child)
    workflow.generate_json()
    return workflow


class TestWorkflowWfFormat:

    @pytest.mark.unit
    def test_metrics_with_shared_task_names(self) -> None:
        # four tasks sharing the same name fan into a fifth one
        workflow = _build_workflow(
            [("split", f"split_ID000000{i}") for i in range(1, 5)] + [("merge", "merge_ID0000005")],
            [(f"split_ID000000{i}", "merge_ID0000005") for i in range(1, 5)])
        metrics = workflow.workflow_json["workflow"]["specification"]["metrics"]
        assert metrics["numberOfLevels"] == 2
        assert metrics["minimumWidth"] == 1
        assert metrics["maximumWidth"] == 4
        SchemaValidator().validate_instance(workflow.workflow_json)

    @pytest.mark.unit
    def test_metrics_of_a_chain_of_same_name_tasks(self) -> None:
        ids = [f"step_ID000000{i}" for i in range(1, 5)]
        workflow = _build_workflow([("step", task_id) for task_id in ids], list(zip(ids, ids[1:])))
        metrics = workflow.workflow_json["workflow"]["specification"]["metrics"]
        assert (metrics["numberOfLevels"], metrics["minimumWidth"], metrics["maximumWidth"]) == (4, 1, 1)
        SchemaValidator().validate_instance(workflow.workflow_json)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "task_ids",
        [
            ["ID0000001", "ID0000002"],  # as in the WfFormat documentation
            ["ingest:raw", "proc/a"],  # ':' and '/' are allowed from WfFormat 1.6
            ["main.cwl#step-1", "main.cwl#step-2"],
            ["stage_in_local_local_0_0", "get_IDs_ID0000003"],
        ]
    )
    def test_task_ids_allowed_by_wfformat(self, task_ids: list[str]) -> None:
        workflow = _build_workflow([("task", task_id) for task_id in task_ids], [(task_ids[0], task_ids[1])])
        SchemaValidator().validate_instance(workflow.workflow_json)

    @pytest.mark.unit
    def test_byte_metrics_are_integers(self) -> None:
        workflow = Workflow(name="bytes_test")
        workflow.add_task(Task(name="a", task_id="a", runtime=0, bytes_read=2, bytes_written=3))
        workflow.add_task(Task(name="b", task_id="b", runtime=0, bytes_read=4, bytes_written=0.5))
        workflow.generate_json()
        metrics = workflow.workflow_json["workflow"]["execution"]["metrics"]
        # totals do not depend on runtimes, and fractional totals are omitted
        assert metrics == {"totalNumBytesRead": 6}
        assert isinstance(metrics["totalNumBytesRead"], int)
        SchemaValidator().validate_instance(workflow.workflow_json)

    @pytest.mark.unit
    def test_execution_fields(self) -> None:
        task = Task(name="a", task_id="a", runtime=1.0, cores=0, args=["-x", "", 3],
                    executed_at="2025-03-01 10:00:01.000000")
        execution = task.execution_as_dict()
        assert "coreCount" not in execution
        assert execution["command"] == {"arguments": ["-x", "3"]}
        assert datetime.fromisoformat(execution["executedAt"]).tzinfo is not None
        assert "command" not in Task(name="b", task_id="b", runtime=1.0).execution_as_dict()

    @pytest.mark.unit
    def test_specification_only_workflow(self) -> None:
        workflow = _build_workflow([("a", "a"), ("b", "b")], [("a", "b")], include_execution=False)
        assert "execution" not in workflow.workflow_json["workflow"]
        SchemaValidator().validate_instance(workflow.workflow_json)

    @pytest.mark.unit
    def test_instance_properties(self) -> None:
        workflow = _build_workflow([("a", "a")], [], author_name="Jane Doe", author_email="jane@example.org",
                                   author_institution="Example University", author_country="US",
                                   created_at="2026-08-04T15:19:28-08:00", executed_at="2026-08-04T15:20:00-08:00")
        assert workflow.workflow_json["author"] == {"name": "Jane Doe", "email": "jane@example.org",
                                                   "institution": "Example University", "country": "US"}
        assert workflow.workflow_json["createdAt"] == "2026-08-04T15:19:28-08:00"
        assert workflow.workflow_json["workflow"]["execution"]["executedAt"] == "2026-08-04T15:20:00-08:00"

    @pytest.mark.unit
    def test_duplicate_file_ids(self) -> None:
        workflow = Workflow(name="files_test")
        workflow.add_task(Task(name="a", task_id="a", runtime=1.0, output_files=[File("f.txt", 10)]))
        workflow.add_task(Task(name="b", task_id="b", runtime=1.0, input_files=[File("f.txt", 20)]))
        workflow.add_dependency("a", "b")
        workflow.generate_json()
        assert workflow.workflow_json["workflow"]["specification"]["files"] == [{"id": "f.txt", "sizeInBytes": 10}]
        SchemaValidator().validate_instance(workflow.workflow_json)
