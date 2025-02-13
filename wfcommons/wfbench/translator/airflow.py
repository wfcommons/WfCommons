#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import re

from logging import Logger
from typing import Optional, Union

from .abstract_translator import Translator
from ...common import Workflow

class AirflowTranslator(Translator):
    """
    A WfFormat parser for creating Nextflow workflow applications.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path],
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 logger: Optional[Logger] = None,
                 input_file_directory: pathlib.Path = pathlib.Path("/")) -> None:
        """Create an object of the translator."""
        super().__init__(workflow, logger)

        self.input_file_directory = input_file_directory
        self._prep_commands()

        self.script = f"""
from __future__ import annotations

from datetime import datetime
from airflow.models.dag import DAG
from airflow.operators.bash import BashOperator

with DAG(
    "{self.workflow.name}_wfcommons",
    description="airflow translation of a wfcommons instance",
    schedule="0 0 * * *",
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=["wfcommons"],
) as dag:
"""

    def translate(self, output_file_path: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description(WfFormat) into a Airflow workflow application.

        : param output_file_path: The name of the output file(e.g., workflow.py).
        : type output_file_path: pathlib.Path
        """

        for task in self.tasks:
            self.script += f"""
    {task} = BashOperator(
        task_id="{task}",
        depends_on_past=False,
        bash_command="{self.task_commands[task]}",
        retries=3,
        )
"""
        for task in self.tasks:
            parents = ", ".join(self.task_parents[task])
            if parents:
                self.script += f"""
    [{parents}] >> {task}
"""

        self._write_output_file(self.script, output_file_path)

    def _prep_commands(self):
        self.task_commands = {}
        for task in self.workflow.workflow_json["workflow"]["execution"]["tasks"]:
            command_str = " ".join([task["command"]["program"]] + task["command"]["arguments"])
            # Prepends { and } with \" (i.e. {hi} -> \"{hi\"}
            command_str = re.sub(r"(\{|\})", r"\"\1", command_str)
            # Prepends txt filenames with absolute path
            command_str = re.sub(r"([\w\-]+\.txt)",
                                 lambda m: f"{self.input_file_directory.absolute().as_posix()}/{m.group(1)}",
                                 command_str)
            self.task_commands[task["id"]] = command_str