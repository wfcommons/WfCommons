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
from ...common import File, Workflow

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
            # Makeflow infers task dependencies from files: the input/output files include
            # dependency marker files when some dependencies are not expressed by files
            input_files = self._task_input_files(task)
            output_files = self._task_output_files(task)

            make_clause = ""
            # output files
            for output_file in output_files:
                make_clause += f"{self._makeflow_file_name(output_file)} "
            make_clause += ": "
            # input files
            for input_file in input_files:
                make_clause += f"{self._makeflow_file_name(input_file)} "
            make_clause += "\n"
            # Command
            make_clause += "\t"
            make_clause += task.program + " "

            input_spec = "\"["
            for file in input_files:
                input_spec += f"\\\\\"data/{file.file_id}\\\\\","
            input_spec = input_spec.rstrip(",") + "]\""

            output_spec = "\"{"
            for file in output_files:
                output_spec += f"\\\\\"data/{file.file_id}\\\\\":{str(file.size)},"
            output_spec = output_spec.rstrip(",") + "}\""

            args = []
            for a in task.args:
                if "--output-files" in a:
                    args.append(f"--output-files {output_spec}")
                elif "--input-files" in a:
                    args.append(f"--input-files {input_spec}")
                else:
                    args.append(self._makeflow_argument(a))
            # output files (e.g., dependency marker files) are created by wfbench
            if output_files and not any("--output-files" in a for a in task.args):
                args.append(f"--output-files {output_spec}")

            args = " ".join(f"{a}" for a in args)
            make_clause += args + "\n"
            self._script += make_clause + "\n\n"
        return

    @staticmethod
    def _makeflow_file_name(file: File) -> str:
        """
        Get the name of a file in a Makeflow rule. An unquoted '#' starts a comment anywhere in a
        Makeflow file, and ':' separates the targets of a rule from its sources.

        :param file: The file.
        :type file: File

        :return: The (quoted, if needed) file name.
        :rtype: str
        """
        file_name = f"data/{file.file_id}"
        return f"\"{file_name}\"" if "#" in file_name or ":" in file_name else file_name

    def _makeflow_argument(self, argument: str) -> str:
        """
        Get a task argument for a Makeflow command. An unquoted '#' starts a comment anywhere in a
        Makeflow command (e.g., in a '--name' value), and quotes are passed to the shell.

        :param argument: A task argument.
        :type argument: str

        :return: The (quoted, if needed) argument.
        :rtype: str
        """
        argument = self._shell_safe_argument(argument)
        if "#" in argument:
            flag, separator, value = argument.partition(" ")
            argument = f"{flag} '{value}'" if separator else f"'{argument}'"
        return argument

    def _write_readme_file(self, output_folder: pathlib.Path) -> None:
        """
        Write the README  file.

        :param output_folder: The path of the output folder.
        :type output_folder: pathlib.Path
        """
        readme_file_path = output_folder.joinpath("README")
        with open(readme_file_path, "w") as out:
            out.write(f"In directory {str(output_folder)}:\n")
            out.write(f"  - The Makeflow input file:  workflow.makeflow\n")
            out.write(f"  - Run the workflow:         makeflow workflow.makeflow\n")
