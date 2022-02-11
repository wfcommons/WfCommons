#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import pathlib

from logging import Logger
from typing import Optional
from unicodedata import category

from .abstract_translator import Translator
from ...common.file import FileLink


class SwiftTTranslator(Translator):
    """
    A WfFormat parser for creating Swift/T workflow applications.

    :param workflow_json_file_path: Path to the workflow benchmark JSON instance.
    :type workflow_json_file_path: pathlib.Path
    :param work_dir: Path to the workflow working directory.
    :type work_dir: pathlib.Path
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 workflow_json_file_path: pathlib.Path,
                 work_dir: pathlib.Path,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow_json_file_path, logger)

        self.work_dir = work_dir
        self.parsed_tasks = []
        self.out_counter = 1
        self.files_map = {}
        self.tasks_map = {}
        self.script = "import files;\nimport io;\nimport json;\nimport unix;\n\n"

        # find applications
        self.apps = []
        for task in self.tasks.values():
            self.tasks_map[task.name] = task.category

            if task.category not in self.apps:
                self.apps.append(task.category)

    def translate(self, output_file_name: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into a Swift/T workflow application.

        :param output_file_name: The name of the output file (e.g., workflow.swift).
        :type output_file_name: pathlib.Path
        """
        # creating apps
        for app in self.apps:
            self.script += "@suppress=unused_output\n"
            self.script += f"app (file output) {app} (float percent_cpu, int cpu_work, string out_name, file inputs[]) "
            self.script += "{\n" \
                "  \"/sw/summit/python/3.8/anaconda3/2020.07-rhel8/bin/python3\" \\\n" \
                f"  \"{self.work_dir}/wfbench.py\" \\\n" \
                f"  \"{app}\" \\\n" \
                "  \"--percent-cpu\" percent_cpu \\\n" \
                "  \"--cpu-work\" cpu_work \\\n" \
                f"  \"--out\" out_name \\\n" \
                f"  inputs \n" \
                "}\n\n"

        # defining input files
        in_count = 0
        for task_name in self.parent_task_names:
            task = self.tasks[task_name]
            for file in task.files:
                if file.link == FileLink.INPUT:
                    self.files_map[file.name] = f"ins[{in_count}]"
                    in_count += 1
        self.script += f"file ins[] = glob(\"{self.work_dir}/*_input.txt\");\n"
        self.script += "\n"

        # adding tasks
        for task_name in self.parent_task_names:
            self._add_task(task_name)

        # write script to file
        self._write_output_file(self.script, output_file_name)

    def _add_task(self, task_name: str, parent_task: Optional[str] = None) -> None:
        """
        Add a task and its dependencies to the workflow.

        :param task_name: name of the task
        :type task_name: str
        :param parent_task: name of the parent task
        :type parent_task: Optional[str]
        """
        # check dependencies
        for parent in self._find_parents(task_name):
            if parent not in self.parsed_tasks:
                return

        if task_name not in self.parsed_tasks:
            task = self.tasks[task_name]

            # find children
            children = self._find_children(task_name)

            # in/output files
            input_files = []
            for file in task.files:
                if file.link == FileLink.OUTPUT:
                    out_file = file.name
                elif file.link == FileLink.INPUT:
                    input_files.append(self.files_map[file.name])

            # arguments
            args = ", ".join([a.split()[1] for a in task.args[1:3]])
            output_name = task.args[3].replace(
                "--out ", "").replace("{", "json_objectify(\"").replace("}", "\")")
            args += f", {output_name}"
            if len(input_files) > 0:
                self.script += f"file in_{self.out_counter}[];"
                f_i = 0
                for f in input_files:
                    self.script += f" in_{self.out_counter}[{f_i}] = {f};"
                    f_i += 1
                self.script += "\n"
                args += f", in_{self.out_counter}"

            self.script += f"file out_{self.out_counter} <\"{self.work_dir}/{out_file}\"> = {self.tasks_map[task_name]}({args});\n"
            self.files_map[out_file] = f"out_{self.out_counter}"
            self.out_counter += 1

            self.parsed_tasks.append(task_name)

            for child_task_name in children:
                self._add_task(child_task_name)
