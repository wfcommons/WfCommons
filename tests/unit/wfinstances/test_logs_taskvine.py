#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib

import pytest

from wfcommons.wfinstances.logs import TaskVineLogsParser


def _write_vine_logs(logs_dir: pathlib.Path) -> None:
    """Write minimal TaskVine logs: a 'cat' task producing a file read by a 'wc' task."""
    (logs_dir / "debug").write_text(
        "2026/01/31 12:00:01.00 vine_manager[245] vine: worker-1 (10.0.0.1:9123) busy on './cat a.txt b.txt'\n"
        "2026/01/31 12:00:01.00 vine_manager[245] vine: Task 1 state change: READY (1) to RUNNING (2)\n"
        "2026/01/31 12:00:02.00 vine_manager[245] vine: worker-1 (10.0.0.1:9123) busy on 'wc -l \"my file.txt\"'\n"
        "2026/01/31 12:00:02.00 vine_manager[245] vine: Task 2 state change: READY (1) to RUNNING (2)\n")
    (logs_dir / "taskgraph").write_text(
        "# taskgraph\n"
        'FILE file-meta-x "bin/cat" 0\n'
        'FILE file-meta-w "wfbench" 0\n'
        'FILE file-meta-a "data/my input.txt" 0\n'
        'FILE file-rnd-c "concatenated.txt" 0\n'
        'FILE file-rnd-m "my file.txt" 0\n'
        'TASK T1 "./cat a.txt b.txt" INPUTS file-meta-x file-meta-w file-meta-a OUTPUTS file-rnd-c\n'
        'TASK T2 "wc -l "my file.txt" INPUTS" INPUTS file-rnd-c OUTPUTS file-rnd-m\n')
    (logs_dir / "transactions").write_text(
        "# time manager_pid MANAGER manager_pid START|END time_from_origin\n"
        "1769907500000000 245 MANAGER 245 START 0\n"
        "1769907502000000 245 TASK 1 RUNNING worker-1 FIRST_RESOURCES {}\n"
        "1769907502000001 245 WORKER worker-1 TRANSFER INPUT file-meta-a 1007 72409 1769907492000000\n"
        "1769907502000002 245 WORKER worker-1 TRANSFER OUTPUT file-rnd-c 2003 72409 1769907492000000\n"
        "1769907506500000 245 TASK 1 DONE SUCCESS 0 {} {}\n"
        "1769907512000000 245 TASK 2 RUNNING worker-1 FIRST_RESOURCES {}\n"
        "1769907512000002 245 WORKER worker-1 TRANSFER OUTPUT file-rnd-m 11 72409 1769907492000000\n"
        "1769907516500000 245 TASK 2 DONE SUCCESS 0 {} {}\n"
        "1769907532000000 245 MANAGER 245 END 0\n")


class TestTaskVineLogsParser:

    @pytest.mark.unit
    def test_default_filenames_to_ignore(self, tmp_path: pathlib.Path) -> None:
        _write_vine_logs(tmp_path)
        parser = TaskVineLogsParser(tmp_path)
        assert parser.filenames_to_ignore == set()

    @pytest.mark.unit
    def test_build_workflow(self, tmp_path: pathlib.Path, capsys) -> None:
        _write_vine_logs(tmp_path)
        parser = TaskVineLogsParser(tmp_path, filenames_to_ignore=["wfbench"], taskvine_version="7.15.0")
        workflow = parser.build_workflow("taskvine-test")
        workflow.generate_json()
        wf_json = workflow.workflow_json

        # nothing is printed to stdout
        assert capsys.readouterr().out == ""

        assert wf_json["runtimeSystem"]["version"] == "7.15.0"
        files = {f["id"]: f["sizeInBytes"] for f in wf_json["workflow"]["specification"]["files"]}
        # names are matched exactly (the 'cat' executable does not hide 'concatenated.txt'), quoted
        # names may contain spaces, and file IDs are sanitized
        assert files == {"data/my_input.txt": 1007, "concatenated.txt": 2003, "my_file.txt": 11}

        spec = {t["id"]: t for t in wf_json["workflow"]["specification"]["tasks"]}
        assert spec["Task_ID0000001"]["inputFiles"] == ["data/my_input.txt"]
        assert spec["Task_ID0000002"]["inputFiles"] == ["concatenated.txt"]
        assert spec["Task_ID0000002"]["parents"] == ["Task_ID0000001"]

        # timestamps (in microseconds) come from the transactions log
        execution = wf_json["workflow"]["execution"]
        assert execution["executedAt"] == "2026-02-01T00:58:20+00:00"
        assert execution["makespanInSeconds"] == 32.0
        tasks = {t["id"]: t for t in execution["tasks"]}
        assert tasks["Task_ID0000001"]["executedAt"] == "2026-02-01T00:58:22+00:00"
        assert tasks["Task_ID0000001"]["runtimeInSeconds"] == 4.5
