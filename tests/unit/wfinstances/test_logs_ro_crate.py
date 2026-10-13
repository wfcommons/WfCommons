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

import networkx as nx
import pytest

from wfcommons.utils import WFFORMAT_ID_PATTERN
from wfcommons.wfinstances.logs import ROCrateLogsParser, StreamflowLogsParser


def _write_crate(crate_dir: pathlib.Path, instrument: str = "shell.cwl") -> None:
    """Write a minimal StreamFlow-like RO-Crate: split -> (blast, blast) -> cat."""
    def step(k, name, inputs, outputs, start, end):
        return {"@id": f"#step-{k}", "@type": "CreateAction", "name": f"Run of workflow/main.cwl#{name}",
                "instrument": {"@id": instrument}, "object": inputs, "result": outputs,
                "startTime": f"2025-03-01T10:{start}+00:00", "endTime": f"2025-03-01T10:{end}+00:00"}

    graph = [
        {"@id": "./", "@type": "Dataset", "mainEntity": {"@id": "main.cwl"}},
        {"@id": "main.cwl", "@type": ["File", "ComputationalWorkflow"], "name": "main.cwl"},
        {"@id": "#main", "@type": "CreateAction", "name": "Run of workflow/main.cwl",
         "startTime": "2025-03-01T10:00:00+00:00", "endTime": "2025-03-01T10:05:00+00:00",
         "object": [], "result": []},
        # a single (non-list) input object
        step(0, "split (chr1, sample+1)", {"@id": "#f-a"}, [{"@id": "#f-b"}, {"@id": "#f-c"}], "00:01", "01:01"),
        step(1, "blast", [{"@id": "#f-b"}], [{"@id": "#f-d"}], "01:02", "02:02"),
        step(2, "blast", [{"@id": "#f-c"}], [{"@id": "#f-e"}], "01:02", "02:10"),
        step(3, "cat", [{"@id": "#f-d"}, {"@id": "#f-e"}], [], "02:11", "03:20"),
        {"@id": "#f-a", "@type": "File", "contentSize": "1000", "alternateName": "my input.txt"},
        {"@id": "#f-b", "@type": "File", "contentSize": "2000", "alternateName": "b.txt"},
        {"@id": "#f-c", "@type": "File", "contentSize": "3000", "alternateName": "c.txt"},
        # two distinct files with the same name
        {"@id": "#f-d", "@type": "File", "contentSize": "4000", "alternateName": "out.txt"},
        {"@id": "#f-e", "@type": "File", "contentSize": "5000", "alternateName": "out.txt"},
        # a connection between the parameters of the same instrument
        {"@id": "#conn", "@type": "ParameterConnection", "sourceParameter": {"@id": f"{instrument}#output"},
         "targetParameter": {"@id": f"{instrument}#input"}},
    ]
    (crate_dir / "ro-crate-metadata.json").write_text(json.dumps({"@graph": graph}))


class TestROCrateLogsParser:

    @pytest.mark.unit
    def test_build_workflow(self, tmp_path: pathlib.Path) -> None:
        _write_crate(tmp_path)
        workflow = StreamflowLogsParser(tmp_path, streamflow_version="0.2.0").build_workflow("crate-test")
        workflow.generate_json()
        wf_json = workflow.workflow_json
        spec = {t["id"]: t for t in wf_json["workflow"]["specification"]["tasks"]}

        # the version is passed through
        assert wf_json["runtimeSystem"]["version"] == "0.2.0"

        # task IDs are sanitized, and numbered with a counter incremented by one
        assert list(spec) == ["main.cwl#split__chr1__sample_1__ID0000001", "main.cwl#blast_ID0000002",
                              "main.cwl#blast_ID0000003", "main.cwl#cat_ID0000004"]

        # distinct files with the same name have distinct (sanitized) IDs
        files = {f["id"]: f["sizeInBytes"] for f in wf_json["workflow"]["specification"]["files"]}
        assert files == {"my_input.txt": 1000, "b.txt": 2000, "c.txt": 3000, "out.txt": 4000, "out.txt_2": 5000}
        assert spec["main.cwl#cat_ID0000004"]["inputFiles"] == ["out.txt", "out.txt_2"]
        for identifier in list(spec) + list(files):
            assert WFFORMAT_ID_PATTERN.match(identifier)

        # a single input object is read as a list
        assert spec["main.cwl#split__chr1__sample_1__ID0000001"]["inputFiles"] == ["my_input.txt"]

        # the connection of an instrument to itself creates neither self-loops nor cycles
        graph = nx.DiGraph()
        for task_id, task in spec.items():
            assert task_id not in task["parents"] and task_id not in task["children"]
            graph.add_edges_from((parent, task_id) for parent in task["parents"])
        assert nx.is_directed_acyclic_graph(graph)
        assert sorted(spec["main.cwl#cat_ID0000004"]["parents"]) == ["main.cwl#blast_ID0000002",
                                                                     "main.cwl#blast_ID0000003"]
        assert spec["main.cwl#blast_ID0000002"]["parents"] == ["main.cwl#split__chr1__sample_1__ID0000001"]

        # task start times come from the crate
        tasks = {t["id"]: t for t in wf_json["workflow"]["execution"]["tasks"]}
        assert tasks["main.cwl#blast_ID0000003"]["executedAt"] == "2025-03-01T10:01:02+00:00"
        assert tasks["main.cwl#blast_ID0000003"]["runtimeInSeconds"] == 68.0

    @pytest.mark.unit
    def test_instruments_to_ignore(self, tmp_path: pathlib.Path) -> None:
        # instruments to ignore are matched by their exact name or by their base name
        _write_crate(tmp_path, instrument="file:///wf/bin/shell.cwl")
        for instruments_to_ignore in [["shell.cwl"], ["file:///wf/bin/shell.cwl"]]:
            parser = ROCrateLogsParser(tmp_path, wms_name="Streamflow", instruments_to_ignore=instruments_to_ignore)
            parser.build_workflow("crate-test")
            assert parser._ignore_instrument("file:///wf/bin/shell.cwl")
            assert not parser._self_connected_instruments
        parser = ROCrateLogsParser(tmp_path, wms_name="Streamflow", instruments_to_ignore=["other.cwl"])
        parser.build_workflow("crate-test")
        assert not parser._ignore_instrument("file:///wf/bin/shell.cwl")
        assert parser._self_connected_instruments == {"file:///wf/bin/shell.cwl"}

    @pytest.mark.unit
    def test_missing_timestamps(self, tmp_path: pathlib.Path) -> None:
        _write_crate(tmp_path)
        metadata = tmp_path / "ro-crate-metadata.json"
        crate = json.loads(metadata.read_text())
        for item in crate["@graph"]:
            item.pop("startTime", None)
            item.pop("endTime", None)
        metadata.write_text(json.dumps(crate))

        workflow = ROCrateLogsParser(tmp_path, wms_name="Streamflow").build_workflow("crate-test")
        workflow.generate_json()
        assert workflow.makespan == 0.0
        assert all(t["runtimeInSeconds"] == 0.0 for t in workflow.workflow_json["workflow"]["execution"]["tasks"])
