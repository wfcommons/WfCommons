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
import json
import ast

from collections import defaultdict
from math import ceil
from logging import Logger
from typing import Dict, List, Optional, Union, MutableSet

from pyparsing import empty

from .abstract_translator import Translator
from ...common import File, FileLink, Workflow
from ...common.task import Task


class NextflowTranslator(Translator):
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

        self.script = ""

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description(WfFormat) into a Nextflow workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """

        # Create the output folder
        output_folder.mkdir(parents=True)

        # Create a topological order of the tasks
        sorted_tasks = self._get_tasks_in_topological_order()
        # print([t.task_id for t in sorted_tasks])

        # Output the code for each task
        for task in sorted_tasks:
            self.script += self._generate_task_code(task)

        # Output the code for the workflow
        self.script += self._generate_workflow_code(sorted_tasks)


        # Output the code to the workflow file
        self._write_output_file(self.script, output_folder.joinpath("workflow.nf"))

        # Create the README file
        self._write_readme_file(output_folder)

        # Create additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)

        return


    def _get_tasks_in_topological_order(self) -> List[Task]:
        levels = {0: self._find_root_tasks()}
        sorted_tasks: List[Task] = levels[0]
        current_level = 1
        while (True):
            tasks_in_current_level = []
            all_children = [self._find_children(p.task_id) for p in levels[current_level-1]]
            all_children = [item for sublist in all_children for item in sublist]
            all_children = list(set(all_children))
            if not all_children:
                break
            for potential_task in all_children:
                if all(parent in sorted_tasks for parent in self._find_parents(potential_task.task_id)):
                    tasks_in_current_level.append(potential_task)
            levels[current_level] = tasks_in_current_level
            sorted_tasks += tasks_in_current_level
            current_level += 1
        return sorted_tasks


    def _generate_task_code(self, task: Task) -> str:
        code = f"process {task.task_id}()" + "{\n"
        code += f"\tinput:\n"
        if self._find_parents(task.task_id):
            for f in task.input_files:
                code += f"\t\tval {f.file_id}\n"
        code += "\n"

        code += f"\toutput:\n"
        if self._find_children(task.task_id):
            for f in task.output_files:
                code += f"\t\tval {f.file_id}\n"
            code += "\n"

        code += "\tscript:\n"

        # Generate output variables
        if self._find_children(task.task_id):
            for f in task.output_files:
                code += "\t\t" + f.file_id + " = \"${params.pwd}/data/" + f.file_id + "\"\n"

        # Generate command
        code += "\t\t\"\"\"\n"
        code += "\t\t\"${params.pwd}/bin/" + task.program + "\""
        for a in task.args:
            if "--output-files" in a:
                flag, output_files_dict = a.split(" ", 1)
                output_files_dict = {str("${params.pwd}/data/" + key): value for key, value in
                                     ast.literal_eval(output_files_dict).items()}
                a = f"{flag} '{json.dumps(output_files_dict)}'"
            elif "--input-files" in a:
                flag, input_files_arr = a.split(" ", 1)
                input_files_arr = [str("${params.pwd}/data/" + file) for file in
                                   ast.literal_eval(input_files_arr)]
                a = f"{flag} '{json.dumps(input_files_arr)}'"

            code += " " + a
        code += "\n"
        code += "\t\t\"\"\"\n"

        code += "}\n\n"
        return code

    def _generate_workflow_code(self, sorted_tasks: List[Task]) -> str:
        code = "workflow {\n"
        for task in sorted_tasks:
            code += self._generate_task_invocation_code(task)
        code += "}\n"
        return code

    def _generate_task_invocation_code(self, task: Task) -> str:

        # Figure out task output values
        if task.output_files and self._find_children(task.task_id):
            output_values = "(" + "\t,\n".join([f.file_id for f in task.output_files]) + ")"
        else:
            output_values = "_"

        # Figure out task input values
        if task.input_files and self._find_parents(task.task_id):
            input_values = "\t,\n".join([f.file_id for f in task.input_files])
        else:
            input_values = ""

        code = output_values + " = " + task.task_id + "(" + input_values + ")\n\n"

        return code


    def _write_readme_file(self, output_folder: pathlib.Path) -> None:
        """
        Write the README  file.

        :param output_folder: The path of the output folder.
        :type output_folder: pathlib.Path
        """
        readme_file_path = output_folder.joinpath("README")
        with open(readme_file_path, "w") as out:
            out.write(f"Run the workflow in directory {str(output_folder)} using the following command:\n")
            out.write(f"nextflow run ./workflow.nf --pwd `pwd`\n")

