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

from wfcommons.wfinstances.logs import NextflowLogsParser


def _write_nextflow_run(run_dir: pathlib.Path) -> None:
    """Write a minimal Nextflow execution directory (nf-prov RO-Crate and execution trace)."""
    graph = [
        {"@id": "./", "@type": "Dataset", "mainEntity": {"@id": "main.nf"}},
        {"@id": "main.nf", "@type": ["File", "ComputationalWorkflow"], "name": "main.nf"},
        {"@id": "#run", "@type": "CreateAction", "name": "Nextflow workflow run 1234",
         "startTime": "2025-03-01T10:00:00.000+01:00", "endTime": "2025-03-01T10:10:00.000+01:00",
         "object": [], "result": []},
        {"@id": "data/in.fa", "@type": "File", "contentSize": "10"},
        {"@id": "#task/aa/out.txt", "@type": "File", "contentSize": "20"},
        {"@id": "#aa", "@type": "CreateAction", "name": "NFCORE:MAP_IDS (1)",
         "instrument": {"@id": "#process-NFCORE:MAP_IDS"},
         "object": [{"@id": "data/in.fa"}], "result": [{"@id": "#task/aa/out.txt"}]},
        {"@id": "#bb", "@type": "CreateAction", "name": "NFCORE:CAT (2)",
         "instrument": {"@id": "#process-NFCORE:CAT"},
         "object": {"@id": "#task/aa/out.txt"}, "result": []},
    ]
    (run_dir / "ro-crate-metadata.json").write_text(json.dumps({"@graph": graph}))
    (run_dir / "execution_trace.txt").write_text(
        "task_id\tname\tsubmit\tduration\n"
        "1\tNFCORE:MAP_IDS (1)\t2025-03-01 10:00:01.100\t519ms\n"
        "2\tNFCORE:CAT (2)\t2025-03-01 10:00:02.000\t1d 2h 3m 4s\n")


class TestNextflowLogsParser:

    @pytest.mark.unit
    @pytest.mark.parametrize("duration, seconds", [
        ("519ms", 0.519),
        ("2m 40s", 160.0),
        ("1h 3m 2s", 3782.0),
        ("1d 2h 3m 4s", 93784.0),
        ("4.5s", 4.5),
        ("1m 1ms", 60.001),
        ("-", 0.0),
    ])
    def test_parse_duration(self, duration: str, seconds: float) -> None:
        assert NextflowLogsParser._parse_duration(duration) == pytest.approx(seconds)

    @pytest.mark.unit
    def test_build_workflow(self, tmp_path: pathlib.Path) -> None:
        _write_nextflow_run(tmp_path)
        workflow = NextflowLogsParser(tmp_path, nextflow_version="24.10.4").build_workflow("nextflow-test")
        workflow.generate_json()
        wf_json = workflow.workflow_json

        assert wf_json["runtimeSystem"]["version"] == "24.10.4"
        tasks = {t["id"]: t for t in wf_json["workflow"]["execution"]["tasks"]}
        # a process name containing "_ID" is fine
        assert set(tasks) == {"NFCORE:MAP_IDS__1__ID0000001", "NFCORE:CAT__2__ID0000002"}
        assert tasks["NFCORE:MAP_IDS__1__ID0000001"]["runtimeInSeconds"] == pytest.approx(0.519)
        assert tasks["NFCORE:CAT__2__ID0000002"]["runtimeInSeconds"] == pytest.approx(93784.0)
        # the trace timestamps are in the timezone of the workflow run
        assert tasks["NFCORE:MAP_IDS__1__ID0000001"]["executedAt"] == "2025-03-01T10:00:01.100000+01:00"
        assert wf_json["workflow"]["execution"]["executedAt"] == "2025-03-01T10:00:00+01:00"
        spec = {t["id"]: t for t in wf_json["workflow"]["specification"]["tasks"]}
        assert spec["NFCORE:CAT__2__ID0000002"]["parents"] == ["NFCORE:MAP_IDS__1__ID0000001"]
