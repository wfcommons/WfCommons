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

from wfcommons.wfinstances.logs import MakeflowLogsParser


def _write_makeflow_run(run_dir: pathlib.Path) -> None:
    """Write a minimal Makeflow execution directory with resource monitor summaries."""
    (run_dir / "data").mkdir()
    (run_dir / "monitor").mkdir()
    (run_dir / "data/a.txt").write_bytes(b"x" * 100)
    (run_dir / "workflow.makeflow").write_text(
        "# Makeflow workflow specification\n"
        "URL=http://example.org/data\n"
        "OUT=data/b.txt\n"
        "\n"
        "$(OUT) : data/a.txt\n"
        "\tg++ -o $(OUT) data/a.txt\n"
        "\n"
        "data/c.txt : ${OUT}->b_remote.txt\n"
        "\twfbench --out data/c.txt\n")
    commands = ["g++ -o data/b.txt data/a.txt", "wfbench --out data/c.txt"]
    (run_dir / "workflow.makeflow.makeflowlog").write_text(
        "# STARTED 1700000000000000\n" +
        "".join(f"# NODE\t{k}\t{c}\n" for k, c in enumerate(commands)) +
        "# FILE 1700000001000500 data/b.txt 1 2048\n"
        "# FILE 1700000002000500 data/c.txt 1 4096\n"
        "# COMPLETED 1700000007123456\n")
    for k, command in enumerate(commands):
        summary = {"command": command, "start": [1700000001000000 + k * 1000000, "us"],
                   "wall_time": [2.0, "s"] if k == 0 else [0, "s"], "cpu_time": [1.0, "s"],
                   "cores": [0, "cores"] if k == 1 else [2, "cores"], "memory": [12, "MB"],
                   "bytes_read": [0.75, "MB"], "bytes_written": [2.5, "MB"], "machine_cpus": [8, "cores"],
                   "host": "worker0", "monitor_version": "7.8.2 FINAL"}
        (run_dir / f"monitor/resource-rule-{k}.summary").write_text(json.dumps(summary))


class TestMakeflowLogsParser:

    @pytest.mark.unit
    def test_build_workflow(self, tmp_path: pathlib.Path) -> None:
        _write_makeflow_run(tmp_path)
        workflow = MakeflowLogsParser(tmp_path, tmp_path / "monitor").build_workflow("makeflow-test")
        workflow.generate_json()
        wf_json = workflow.workflow_json

        # the version comes from the resource monitor (makeflow and resource monitor are part of cctools)
        assert wf_json["runtimeSystem"]["version"] == "7.8.2"

        # task IDs are sanitized, variables are expanded, and renamed files use their workflow-level name
        spec = {t["id"]: t for t in wf_json["workflow"]["specification"]["tasks"]}
        assert list(spec) == ["g___ID0000001", "wfbench_ID0000002"]
        assert spec["g___ID0000001"]["outputFiles"] == ["data/b.txt"]
        assert spec["wfbench_ID0000002"]["inputFiles"] == ["data/b.txt"]
        assert spec["wfbench_ID0000002"]["parents"] == ["g___ID0000001"]
        files = {f["id"]: f["sizeInBytes"] for f in wf_json["workflow"]["specification"]["files"]}
        assert files == {"data/a.txt": 100, "data/b.txt": 2048, "data/c.txt": 4096}

        # resource monitor values are converted from MB into bytes
        tasks = {t["id"]: t for t in wf_json["workflow"]["execution"]["tasks"]}
        assert tasks["g___ID0000001"]["memoryInBytes"] == 12 * 1024 ** 2
        assert tasks["g___ID0000001"]["readBytes"] == 786432
        assert tasks["g___ID0000001"]["writtenBytes"] == int(2.5 * 1024 ** 2)
        assert tasks["g___ID0000001"]["coreCount"] == 2
        assert tasks["g___ID0000001"]["avgCPU"] == 50.0
        assert tasks["g___ID0000001"]["executedAt"] == "2023-11-14T22:13:21+00:00"
        # a task runs on at least one core, and a zero runtime does not break the CPU utilization
        assert tasks["wfbench_ID0000002"]["coreCount"] == 1
        assert "avgCPU" not in tasks["wfbench_ID0000002"]

        # machines are emitted (once)
        assert tasks["g___ID0000001"]["machines"] == ["worker0"]
        assert wf_json["workflow"]["execution"]["machines"] == [{"nodeName": "worker0", "cpu": {"coreCount": 8}}]
