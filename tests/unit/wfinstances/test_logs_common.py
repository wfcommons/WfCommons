#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging
import pathlib
import sqlite3

import pytest

from wfcommons.utils import WFFORMAT_ID_PATTERN
from wfcommons.wfinstances.logs import MakeflowLogsParser, NextflowLogsParser, PegasusLogsParser, \
    ROCrateLogsParser, SnakemakeLogsParser, StreamflowLogsParser, TaskVineLogsParser
from wfcommons.wfinstances.logs.abstract_logs_parser import LogsParser


class _DummyLogsParser(LogsParser):
    def build_workflow(self, workflow_name=None):
        return None


class TestLogsParsersCommon:

    @pytest.mark.unit
    def test_exports(self) -> None:
        import wfcommons.wfinstances as wfinstances
        for name in ["MakeflowLogsParser", "NextflowLogsParser", "PegasusLogsParser", "ROCrateLogsParser",
                     "SnakemakeLogsParser", "StreamflowLogsParser", "TaskVineLogsParser"]:
            assert hasattr(wfinstances, name)

    @pytest.mark.unit
    def test_make_file_id(self) -> None:
        parser = _DummyLogsParser("dummy")

        # invalid characters are replaced, valid IDs are kept
        assert parser._make_file_id("data/f.txt") == "data/f.txt"
        assert parser._make_file_id("sample A=1+x.fa") == "sample_A_1_x.fa"
        assert parser._make_file_id("https://example.org/get?id=1&x=%20") == "https://example.org/get_id_1_x__20"
        assert parser._make_file_id("data/d.txt->d_remote.txt") == "data/d.txt-_d_remote.txt"
        for file_id in parser._file_keys:
            assert WFFORMAT_ID_PATTERN.match(file_id)

        # the same file always gets the same ID
        assert parser._make_file_id("sample A=1+x.fa") == "sample_A_1_x.fa"

        # distinct files whose names collide get distinct IDs
        assert parser._make_file_id("out.txt", key="#f-d") == "out.txt"
        assert parser._make_file_id("out.txt", key="#f-e") == "out.txt_2"
        assert parser._make_file_id("out.txt", key="#f-d") == "out.txt"
        assert parser._make_file_id("sample A+1+x.fa") == "sample_A_1_x.fa_2"

    @pytest.mark.unit
    def test_constructors_keyword_wiring(self, tmp_path: pathlib.Path) -> None:
        logger = logging.getLogger("test_logs_common")

        # Pegasus
        (tmp_path / "pegasus").mkdir()
        parser = PegasusLogsParser(tmp_path / "pegasus", description="pegasus run", logger=logger)
        assert (parser.wms_name, parser.wms_version, parser.wms_url) == ("Pegasus", None, "https://pegasus.isi.edu")
        assert parser.description == "pegasus run" and parser.logger is logger

        # Makeflow
        mf = tmp_path / "makeflow"
        (mf / "monitor").mkdir(parents=True)
        (mf / "wf.makeflow").write_text("b.txt : a.txt\n\tcat a.txt > b.txt\n")
        (mf / "wf.makeflow.makeflowlog").write_text("# NODE\t0\tcat a.txt > b.txt\n")
        parser = MakeflowLogsParser(mf, mf / "monitor", description="makeflow run", logger=logger)
        assert (parser.wms_name, parser.wms_version) == ("Makeflow", None)
        assert parser.wms_url.startswith("https://ccl.cse.nd.edu/")
        assert parser.description == "makeflow run" and parser.logger is logger

        # Snakemake
        (tmp_path / "snakemake").mkdir()
        sqlite3.connect(tmp_path / "snakemake" / "snkmt.sqlite").close()
        parser = SnakemakeLogsParser(tmp_path / "snakemake", tmp_path / "snakemake" / "snkmt.sqlite",
                                     description="snakemake run", logger=logger, snakemake_version="9.1.0")
        assert (parser.wms_name, parser.wms_version) == ("Snakemake", "9.1.0")
        assert parser.wms_url == "https://snakemake.readthedocs.io/"
        assert parser.description == "snakemake run" and parser.logger is logger

        # TaskVine
        tv = tmp_path / "vine-logs"
        tv.mkdir()
        for name in ["debug", "taskgraph", "transactions"]:
            (tv / name).write_text("")
        parser = TaskVineLogsParser(tv, ["wfbench"], "taskvine run", logger, taskvine_version="7.15.0")
        assert (parser.wms_name, parser.wms_version) == ("TaskVine", "7.15.0")
        assert parser.wms_url == "https://ccl.cse.nd.edu/software/taskvine/"
        assert parser.description == "taskvine run" and parser.logger is logger

        # Nextflow, Streamflow, and RO-Crate
        nf = tmp_path / "nextflow"
        nf.mkdir()
        (nf / "execution_trace.txt").write_text("")
        parser = NextflowLogsParser(nf, "24.10.4", description="nextflow run", logger=logger)
        assert (parser.wms_name, parser.wms_version, parser.wms_url) == ("Nextflow", "24.10.4", "https://www.nextflow.io")
        assert parser.description == "nextflow run" and parser.logger is logger
        parser = StreamflowLogsParser(nf, "0.2.0", description="streamflow run", logger=logger)
        assert (parser.wms_name, parser.wms_version) == ("Streamflow", "0.2.0")
        assert parser.description == "streamflow run" and parser.logger is logger
        (nf / "ro-crate-metadata.json").write_text("{}")
        parser = ROCrateLogsParser(nf, "Other", "1.0", "https://example.org", "crate run", logger)
        assert (parser.wms_name, parser.wms_version, parser.wms_url) == ("Other", "1.0", "https://example.org")
        assert parser.description == "crate run" and parser.logger is logger
