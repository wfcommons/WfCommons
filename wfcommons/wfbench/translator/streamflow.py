#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import shutil
import logging
import pathlib
import ast
import json
from typing import Union, Optional
from collections import defaultdict, deque


from .abstract_translator import Translator
from ...common import Workflow, Task, File

this_dir = pathlib.Path(__file__).resolve().parent


class StreamflowTranslator(Translator):
    """
    A WfFormat parser for creating Streamflow workflow benchmarks.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path],
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """
    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 logger: Optional[logging.Logger] = None) -> None:
        super().__init__(workflow, logger)

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        # Perform the CWL translation (which will create the output folder)
        from wfcommons.wfbench import CWLTranslator
        cwl_translator = CWLTranslator(workflow=self.workflow, logger=self.logger)
        cwl_translator.translate(output_folder)

        # Generate the streamflow.yml file
        self._generate_streamflow_file(output_folder)

        # Generate the inputs.yml file
        self._generate_inputs_file(output_folder)

        # Generate the README file
        self._write_readme_file(output_folder)

    def _generate_streamflow_file(self, output_folder: pathlib.Path) -> None:
        shutil.copy(this_dir.joinpath("templates/streamflow/streamflow.yml"), output_folder)

    def _generate_inputs_file(self, output_folder: pathlib.Path) -> None:

        file_path = output_folder.joinpath("inputs.yml")
        with open(file_path, "w") as out:
            for task_id in self.workflow.roots():
                task = self.workflow.tasks[task_id]
                out.write(f"{task.task_id}_input:\n")
                for f in task.input_files:
                    out.write(f"  - class: File\n")
                    out.write(f"    path: ./data/{f}")




    def _write_readme_file(self, output_folder: pathlib.Path) -> None:
        readme_file_path = output_folder.joinpath("README")
        with open(readme_file_path, "w") as out:
            out.write(f"In directory {str(output_folder)}: streamflow run ./streamflow.yml\n")