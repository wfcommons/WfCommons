#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import pathlib
import pytest

from wfcommons.version import __schema_version__
from wfcommons.wfinstances import Instance, SchemaValidator


class TestInstance:

    @pytest.fixture
    def instance(self) -> pathlib.Path:
        instance = {
            "name": "workflow_test",
            "description": "Instance generate for WfCommons Test",
            "createdAt": "2020-12-30T02:19:01.238077",
            "schemaVersion": "1.6",
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
                            "children": [],
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
    def test_instance_creation(self, instance : pathlib.Path) -> None:
        instance = Instance(instance)
        assert(instance.leaves() == instance.roots())


def _wfformat_example() -> dict:
    """The example instance of the WfFormat 1.6 documentation."""
    return {
        "name": "my-workflow-instance",
        "description": "A two-task workflow executed on a single node",
        "createdAt": "2026-08-04T15:19:28-08:00",
        "schemaVersion": "1.6",
        "author": {"name": "Jane Doe", "email": "jane.doe@example.org",
                   "institution": "Example University", "country": "US"},
        "runtimeSystem": {"name": "pegasus", "version": "5.0.6", "url": "https://pegasus.isi.edu"},
        "workflow": {
            "specification": {
                "tasks": [
                    {"name": "split", "id": "ID0000001", "parents": [], "children": ["ID0000002"],
                     "inputFiles": ["input.txt"], "outputFiles": ["chunk.txt"]},
                    {"name": "analyze", "id": "ID0000002", "parents": ["ID0000001"], "children": [],
                     "inputFiles": ["chunk.txt"], "outputFiles": ["result.txt"]}
                ],
                "files": [
                    {"id": "input.txt", "sizeInBytes": 1024}, {"id": "chunk.txt", "sizeInBytes": 512},
                    {"id": "result.txt", "sizeInBytes": 128}
                ],
                "metrics": {"numberOfTasks": 2, "numberOfFiles": 3, "sumOfFileSizesInBytes": 1664,
                            "numberOfLevels": 2, "minimumWidth": 1, "maximumWidth": 1}
            },
            "execution": {
                "makespanInSeconds": 14.2,
                "executedAt": "2026-08-04T15:20:00-08:00",
                "tasks": [
                    {"id": "ID0000001", "runtimeInSeconds": 4.2, "executedAt": "2026-08-04T15:20:00-08:00",
                     "command": {"program": "split", "arguments": ["-n", "2", "input.txt"]},
                     "coreCount": 1, "avgCPU": 93.78, "readBytes": 1024, "writtenBytes": 512,
                     "memoryInBytes": 4194304, "energyInKWh": 0.25, "avgPowerInW": 120.0,
                     "machines": ["node01.example.org"]},
                    {"id": "ID0000002", "runtimeInSeconds": 10.0, "executedAt": "2026-08-04T15:20:04-08:00",
                     "coreCount": 4, "readBytes": 512, "writtenBytes": 128, "machines": ["node01.example.org"]}
                ],
                "machines": [
                    {"system": "linux", "architecture": "x86_64", "nodeName": "node01.example.org",
                     "release": "5.15.0", "memoryInBytes": 34359738368,
                     "cpu": {"coreCount": 8, "speedInMHz": 3600, "vendor": "AMD"}}
                ],
                "metrics": {"sumTaskRuntimesInSeconds": 14.2, "totalNumBytesRead": 1536,
                            "totalNumBytesWritten": 640}
            }
        }
    }


def _load(instance: dict, tmp_path: pathlib.Path) -> Instance:
    instance_path = tmp_path / "instance.json"
    instance_path.write_text(json.dumps(instance))
    return Instance(instance_path)


class TestInstanceWfFormat:

    @pytest.mark.unit
    def test_round_trip(self, tmp_path: pathlib.Path) -> None:
        original = _wfformat_example()
        output_path = tmp_path / "written.json"
        _load(original, tmp_path).workflow.write_json(output_path)
        written = json.loads(output_path.read_text())
        assert written == original
        SchemaValidator().validate_instance(written)

    @pytest.mark.unit
    def test_specification_only(self, tmp_path: pathlib.Path) -> None:
        original = _wfformat_example()
        del original["workflow"]["execution"]
        instance = _load(original, tmp_path)
        assert len(instance.workflow.tasks) == 2
        instance.workflow.generate_json()
        assert "execution" not in instance.workflow.workflow_json["workflow"]
        SchemaValidator().validate_instance(instance.workflow.workflow_json)

    @pytest.mark.unit
    def test_minimal_instance(self, tmp_path: pathlib.Path) -> None:
        # only name, schemaVersion, and workflow are required, and task IDs may hold ':' and '/'
        minimal = {
            "name": "minimal",
            "schemaVersion": "1.6",
            "workflow": {
                "specification": {"tasks": [
                    {"name": "t", "id": "t:1", "parents": [], "children": ["t/2"]},
                    {"name": "t", "id": "t/2", "parents": ["t:1"], "children": []}
                ]},
                "execution": {"makespanInSeconds": 1.0, "executedAt": "2026-08-04T15:20:00-08:00",
                              "tasks": [{"id": "t:1", "runtimeInSeconds": 0.5}]}
            }
        }
        instance = _load(minimal, tmp_path)
        assert instance.roots() == ["t:1"] and instance.leaves() == ["t/2"]
        instance.workflow.generate_json()
        written = instance.workflow.workflow_json
        SchemaValidator().validate_instance(written)
        # optional properties that the instance does not declare are not made up
        assert not {"description", "createdAt", "author", "runtimeSystem"} & written.keys()
        assert "executedAt" not in written["workflow"]["execution"]["tasks"][1]

    @pytest.mark.unit
    def test_unused_machine(self, tmp_path: pathlib.Path) -> None:
        original = _wfformat_example()
        original["workflow"]["execution"]["machines"].append({"nodeName": "node02.example.org"})
        instance = _load(original, tmp_path)
        instance.workflow.generate_json()
        machines = instance.workflow.workflow_json["workflow"]["execution"]["machines"]
        assert machines == original["workflow"]["execution"]["machines"]

    @pytest.mark.unit
    def test_machine_without_cpu(self, tmp_path: pathlib.Path) -> None:
        original = _wfformat_example()
        del original["workflow"]["execution"]["machines"][0]["cpu"]
        instance = _load(original, tmp_path)
        assert instance.machines["node01.example.org"].cpu_cores is None

    @pytest.mark.unit
    def test_integral_float_file_size(self, tmp_path: pathlib.Path) -> None:
        # JSON Schema draft 2020-12 considers 1024.0 to be an integer
        original = _wfformat_example()
        original["workflow"]["specification"]["files"][0]["sizeInBytes"] = 1024.0
        _load(original, tmp_path)

    @pytest.mark.unit
    def test_invalid_metrics(self, tmp_path: pathlib.Path) -> None:
        original = _wfformat_example()
        original["workflow"]["specification"]["metrics"]["maximumWidth"] = 2
        with pytest.raises(RuntimeError, match="semantic"):
            _load(original, tmp_path)

    @pytest.mark.unit
    def test_previous_schema_version(self, tmp_path: pathlib.Path) -> None:
        original = _wfformat_example()
        original["schemaVersion"] = "1.5"
        with pytest.raises(RuntimeError, match="syntax"):
            _load(original, tmp_path)

    @pytest.mark.unit
    def test_bundled_schema(self, tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
        # a schema file in the working directory (e.g., left by an older WfCommons) is ignored
        (tmp_path / "wfcommons-schema.json").write_text(json.dumps({"type": "object", "required": ["bogus"]}))
        monkeypatch.chdir(tmp_path)
        validator = SchemaValidator()
        assert validator.schema["properties"]["schemaVersion"]["enum"] == [__schema_version__]
        validator.validate_instance(_wfformat_example())
