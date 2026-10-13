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

from wfcommons.wfinstances.logs import PegasusLogsParser

_KICKSTART = """\
- invocation: True
  version: 3.0
  start: 2023-11-05T14:4{k}:30.123-08:00
  duration: 60.5
  transformation: "diamond::{name}:4.0"
  mainjob:
    start: 2023-11-05T14:4{k}:30.200-08:00
    duration: {duration}
    usage:
      utime: 59.5
      stime: 0.5
      maxrss: 2048
    procs:
      - rchar: 1024
        wchar: 2048
        rbytes: 0
        wbytes: 4096
  machine:
    uname_system: {system}
    uname_nodename: {node}.example.org
    uname_release: 5.15.0
    uname_machine: x86_64
    ram_total: 16384000
    cpu_count: 8
    cpu_speed: 2400
    cpu_vendor: GenuineIntel
"""


def _write_submit_dir(submit_dir: pathlib.Path) -> None:
    """Write a minimal Pegasus 5 submit directory with two compute jobs and a stage-in job."""
    (submit_dir / "00/00").mkdir(parents=True)
    (submit_dir / "braindump.yml").write_text(
        "pegasus_wf_name: diamond-0\ntimestamp: 20231105T144500-0800\nplanner_version: 5.0.4\n")
    (submit_dir / "workflow.yml").write_text("""\
pegasus: '5.0'
name: diamond
x-pegasus: {apiLang: python, createdBy: user, createdOn: 11-05-23T14:44:58Z}
jobs:
- type: job
  name: preprocess
  id: ID0000001
  arguments: [-a, preprocess, -T, 60, '', -i, sample A=1+x.fa, -o, f.b]
  uses:
  - {lfn: sample A=1+x.fa, type: input}
  - {lfn: f.b, type: output, stageOut: true}
- type: job
  name: analyze
  id: ID0000002
  profiles:
    pegasus: {cores: '4'}
  arguments: [-a, analyze, -i, f.b]
  uses:
  - {lfn: f.b, type: input}
""")
    jobs = ["stage_in_local_local_0_0", "preprocess_ID0000001", "analyze_ID0000002"]
    (submit_dir / "diamond-0.dag").write_text(
        "".join(f"JOB {j} 00/00/{j}.sub\n" for j in jobs) +
        "PARENT stage_in_local_local_0_0 CHILD preprocess_ID0000001\n"
        "PARENT preprocess_ID0000001 CHILD analyze_ID0000002\n")
    (submit_dir / "diamond-0.dag.dagman.out").write_text(
        "11/05/23 14:45:01 ** condor_scheduniv_exec.1.0 (CONDOR_DAGMAN) STARTING UP\n"
        "11/05/23 14:52:33 **** condor_scheduniv_exec.1.0 (condor_DAGMAN) EXITING WITH STATUS 0\n")
    (submit_dir / "00/00/preprocess_ID0000001.meta").write_text(json.dumps(
        [{"_id": "sample A=1+x.fa", "_attributes": {"size": "1024"}}, {"_id": "f.b", "_attributes": {"size": "2048"}}]))
    (submit_dir / "00/00/preprocess_ID0000001.out.000").write_text(
        _KICKSTART.format(k=6, name="preprocess", duration=60.0, system="linux", node="worker0"))
    (submit_dir / "00/00/analyze_ID0000002.out.000").write_text(
        _KICKSTART.format(k=7, name="analyze", duration=0.0, system="darwin", node="worker1"))


class TestPegasusLogsParser:

    @pytest.mark.unit
    def test_build_workflow(self, tmp_path: pathlib.Path, monkeypatch) -> None:
        _write_submit_dir(tmp_path / "submit")
        monkeypatch.chdir(tmp_path)
        workflow = PegasusLogsParser(tmp_path / "submit").build_workflow()
        workflow.generate_json()
        wf_json = workflow.workflow_json

        # workflow name, planner version, and DAGMan start time
        assert wf_json["name"] == "diamond-0"
        assert wf_json["runtimeSystem"]["version"] == "5.0.4"
        assert wf_json["runtimeSystem"]["url"] == "https://pegasus.isi.edu"
        execution = wf_json["workflow"]["execution"]
        assert execution["executedAt"] == "2023-11-05T14:45:01-08:00"
        assert execution["makespanInSeconds"] == 452.0

        # file IDs are sanitized, and their sizes come from the meta files
        files = {f["id"]: f["sizeInBytes"] for f in wf_json["workflow"]["specification"]["files"]}
        assert files == {"sample_A_1_x.fa": 1024, "f.b": 2048}

        tasks = {t["id"]: t for t in execution["tasks"]}
        assert set(tasks) == {"preprocess_ID0000001", "analyze_ID0000002"}
        preprocess = tasks["preprocess_ID0000001"]
        # at least one core, or the 'cores' profile
        assert preprocess["coreCount"] == 1
        assert tasks["analyze_ID0000002"]["coreCount"] == 4
        # arguments are non-empty strings
        assert preprocess["command"]["arguments"] == ["-a", "preprocess", "-T", "60", "-i", "sample A=1+x.fa",
                                                      "-o", "f.b"]
        # kickstart start time and maxrss (in KB)
        assert preprocess["executedAt"] == "2023-11-05T14:46:30.200000-08:00"
        assert preprocess["memoryInBytes"] == 2048 * 1024
        assert preprocess["avgCPU"] == 100.0
        # a zero duration does not break the CPU utilization
        assert tasks["analyze_ID0000002"]["runtimeInSeconds"] == 0.0
        assert "avgCPU" not in tasks["analyze_ID0000002"]

        # machines are emitted, and the total RAM (in KB) is converted into bytes
        assert preprocess["machines"] == ["worker0.example.org"]
        machines = {m["nodeName"]: m for m in execution["machines"]}
        assert set(machines) == {"worker0.example.org", "worker1.example.org"}
        assert machines["worker0.example.org"]["memoryInBytes"] == 16384000 * 1024
        assert machines["worker0.example.org"]["system"] == "linux"
        assert machines["worker1.example.org"]["system"] == "macos"

    @pytest.mark.unit
    def test_auxiliary_jobs(self, tmp_path: pathlib.Path, monkeypatch) -> None:
        _write_submit_dir(tmp_path / "submit")
        monkeypatch.chdir(tmp_path)
        workflow = PegasusLogsParser(tmp_path / "submit", ignore_auxiliary=False).build_workflow("pegasus-test")
        workflow.generate_json()
        spec = {t["id"]: t for t in workflow.workflow_json["workflow"]["specification"]["tasks"]}
        assert spec["preprocess_ID0000001"]["parents"] == ["stage_in_local_local_0_0"]
        tasks = {t["id"]: t for t in workflow.workflow_json["workflow"]["execution"]["tasks"]}
        assert tasks["stage_in_local_local_0_0"]["coreCount"] == 1

    @pytest.mark.unit
    def test_invalid_kickstart_record(self, tmp_path: pathlib.Path, monkeypatch) -> None:
        _write_submit_dir(tmp_path / "submit")
        (tmp_path / "submit/00/00/preprocess_ID0000001.out.000").write_text("- invocation: True\n  mainjob: [x\n")
        monkeypatch.chdir(tmp_path)
        with pytest.raises(ValueError):
            PegasusLogsParser(tmp_path / "submit").build_workflow("pegasus-test")
        # the temporary file is removed
        assert not list(tmp_path.glob(".pegasus-parser-tmp-*"))
