#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import shutil

from logging import Logger
from typing import Optional, Union

from .abstract_translator import Translator
from ...common import Workflow

this_dir = pathlib.Path(__file__).resolve().parent


class MakeflowTranslator(Translator):
    """
    A WfFormat parser for creating Makeflow workflow applications.

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
        self._script = ""

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """

        # Generate code
        self._generate_code()

        # write benchmark files
        output_folder.mkdir(parents=True)
        with open(output_folder.joinpath("workflow.makeflow"), "w") as fp:
            fp.write(self._script)

        # additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)

        # README file
        self._write_readme_file(output_folder)

    def _generate_code(self):
        """
        Generate the Makeflow code

        :return: the code
        :rtype: str
        """
        self._script = "# Makeflow workflow specification\n\n"
        for task_name, task in self.workflow.tasks.items():
            make_clause = ""
            # output files
            for output_file in task.output_files:
                make_clause += f"{output_file.file_id} "
            make_clause += ": "
            # input files
            for input_file in task.input_files:
                make_clause += f"{input_file.file_id} "
            make_clause += "\n"
            make_clause += "\t"
            make_clause += task.program + " " + " ".join(task.args)
            make_clause += "\n"
            self._script += make_clause + "\n\n"
        return

        # OLD CODE FOR TASK VINE...
        # args = []
        # for a in task.args:
        #     if "--output-files" in a:
        #         args.append(f"--output-files {output_spec}")
        #     elif "--input-files" in a:
        #         args.append(f"--input-files {input_spec}")
        #     else:
        #         args.append(a)
        # args = " ".join(f"{a}" for a in args)


    def _write_readme_file(self, output_folder: pathlib.Path) -> None:
        """
        Write the README  file.

        :param output_folder: The path of the output folder.
        :type output_folder: pathlib.Path
        """
        readme_file_path = output_folder.joinpath("README")
        with open(readme_file_path, "w") as out:
            out.write(f"In directory {str(output_folder)}:\n")
            out.write(f"  - The Makeflow input file: "
                      f"        workflow.makeflow\n")
            out.write(f"  - Run the workflow: "
                      f"        makeflow workflow.makeflow\n")
