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
import sqlite3

import pytest

from wfcommons.utils import sanitize_wfformat_id, to_rfc3339
from wfcommons.wfinstances.logs import SnakemakeLogsParser


def _write_snkmt_db(run_dir: pathlib.Path) -> pathlib.Path:
    """Write a minimal snkmt database (with the column positions read by the parser) and data files."""
    data = run_dir / "data"
    (data / "sample~A,rep=1").mkdir(parents=True)
    for name, size in [("a.txt", 100), ("sample~A,rep=1/b.txt", 200)]:
        (data / name).write_bytes(b"x" * size)
    db = run_dir / "snkmt.sqlite"
    conn = sqlite3.connect(db)
    conn.execute("CREATE TABLE workflows (id TEXT, snakefile TEXT, started_at TEXT, status TEXT, end_time TEXT)")
    conn.execute("CREATE TABLE rules (id INTEGER, name TEXT, workflow_id TEXT)")
    conn.execute("CREATE TABLE jobs (id INTEGER, workflow_id TEXT, snakemake_id INTEGER, rule_id INTEGER, "
                 "message TEXT, wildcards TEXT, reason TEXT, resources TEXT, shellcmd TEXT, threads INTEGER, "
                 "priority INTEGER, status TEXT, started_at TEXT, end_time TEXT)")
    conn.execute("CREATE TABLE files (id INTEGER, path TEXT, file_type TEXT, job_id INTEGER)")
    conn.execute("INSERT INTO workflows VALUES ('u1', 'Snakefile', '2025-03-01 10:00:00', 'RUNNING', NULL)")
    conn.executemany("INSERT INTO rules VALUES (?, ?, 'u1')", [(1, "split"), (2, "get_IDs")])
    conn.executemany("INSERT INTO jobs VALUES (?, 'u1', ?, ?, NULL, NULL, NULL, NULL, ?, ?, 0, ?, ?, ?)", [
        (1, 1, 1, "wfbench --name split", 2, "SUCCESS", "2025-03-01 10:00:01.000000", "2025-03-01 10:01:01.500000"),
        (2, 2, 2, None, 1, "FAILED", "2025-03-01 10:01:02.000000", None),
    ])
    conn.executemany("INSERT INTO files VALUES (?, ?, ?, ?)", [
        (1, f"{data}/a.txt (access: read-only)", "INPUT", 1),
        (2, f"{data}/sample~A,rep=1/b.txt", "OUTPUT", 1),
        (3, f"{data}/sample~A,rep=1/b.txt", "INPUT", 2),
        (4, f"{run_dir}/logs/1.log", "LOG", 1),
    ])
    conn.commit()
    conn.close()
    return db


class TestSnakemakeLogsParser:

    @pytest.mark.unit
    def test_build_workflow(self, tmp_path: pathlib.Path) -> None:
        db = _write_snkmt_db(tmp_path)
        workflow = SnakemakeLogsParser(tmp_path, snkmt_db=db, snakemake_version="9.1.0").build_workflow("smk-test")
        workflow.generate_json()
        wf_json = workflow.workflow_json

        assert wf_json["runtimeSystem"]["version"] == "9.1.0"

        # file IDs are sanitized
        files = {f["id"]: f["sizeInBytes"] for f in wf_json["workflow"]["specification"]["files"]}
        assert files == {sanitize_wfformat_id(f"{tmp_path}/data/a.txt"): 100,
                         sanitize_wfformat_id(f"{tmp_path}/data/sample~A,rep=1/b.txt"): 200}
        assert any(file_id.endswith("/data/sample_A_rep_1/b.txt") for file_id in files)

        # a rule name containing "_ID" is fine
        spec = {t["id"]: t for t in wf_json["workflow"]["specification"]["tasks"]}
        assert spec["get_IDs_ID0000002"]["parents"] == ["split_ID0000001"]

        # the workflow start time comes from the database, and the makespan from the last job end time
        # (the workflow has no end time)
        execution = wf_json["workflow"]["execution"]
        assert execution["executedAt"] == to_rfc3339("2025-03-01 10:00:00")
        assert execution["makespanInSeconds"] == 61.5

        # a job with no end time does not break the parser
        tasks = {t["id"]: t for t in execution["tasks"]}
        assert tasks["split_ID0000001"]["runtimeInSeconds"] == 60.5
        assert tasks["split_ID0000001"]["executedAt"] == to_rfc3339("2025-03-01 10:00:01")
        assert tasks["split_ID0000001"]["coreCount"] == 2
        assert tasks["get_IDs_ID0000002"]["runtimeInSeconds"] == 0.0
        assert tasks["get_IDs_ID0000002"]["executedAt"] == to_rfc3339("2025-03-01 10:01:02")
