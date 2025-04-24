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
import ast
import json

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
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow, logger)

        self.script = f"""
from __future__ import annotations

import os
from datetime import datetime
from airflow.models.dag import DAG
from airflow.operators.bash import BashOperator

with DAG(
    "{self.workflow.name}",
    description="airflow translation of a wfcommons instance",
    schedule="0 0 * * *",
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=["wfcommons"],
) as dag:
"""

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description(WfFormat) into an Airflow workflow application.

        :param output_folder: The name of the output folder.
        :type output_folder: pathlib.Path
        """

        self._prep_commands(output_folder)

        for task in self.tasks.values():
            self.script += f"""
    {task.task_id} = BashOperator(
        task_id="{task.task_id}",
        depends_on_past=False,
        bash_command='{self.task_commands[task.task_id]}',
        env={{"AIRFLOW_HOME": os.environ["AIRFLOW_HOME"]}},
        retries=3,
        )
"""
        for task in self.tasks.values():
            parents = ", ".join(self.task_parents[task.task_id])
            if parents:
                self.script += f"""
    [{parents}] >> {task.task_id}
"""
        # write benchmark files
        output_folder.mkdir(parents=True)
        with open(output_folder.joinpath("workflow.py"), "w") as fp:
            fp.write(self.script)

        # additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)

    def _prep_commands(self, output_folder: pathlib.Path) -> None:
        """
        Prepares the bash_command strings for the BashOperators.

        :param output_folder: The name of the output folder.
        :type output_folder: pathlib.Path
        """
        self.task_commands = {}

        for task in self.tasks.values():
            program = task.program
            args = []
            for a in task.args:
                if "--output-files" in a:
                    flag, output_files_dict = a.split(" ", 1)
                    output_files_dict = {str(f"${{AIRFLOW_HOME}}/dags/{output_folder.name}/data/{key}"): value for key, value in ast.literal_eval(output_files_dict).items()}
                    a = f"{flag} {json.dumps(output_files_dict)}"
                elif "--input-files" in a:
                    flag, input_files_arr = a.split(" ", 1)
                    input_files_arr = [str(f"${{AIRFLOW_HOME}}/dags/{output_folder.name}/data/{file}") for file in ast.literal_eval(input_files_arr)]
                    a = f"{flag} {json.dumps(input_files_arr)}"
                else:
                    a = a.replace("'", "\"")
                args.append(a)

            command_str = " ".join([str(program)] + args)

            # Escapes all double quotes
            command_str = command_str.replace('"', '\\\\"')

            # Wraps --output-files and --input-files arguments in double quotes
            command_str = re.sub(
                r'(--output-files) (\{.*\}) (--input-files) (\[.*?\])',
                lambda m: f'{m.group(1)} "{m.group(2)}" {m.group(3)} "{m.group(4)}"',
                command_str
            )

            self.task_commands[task.task_id] = command_str

