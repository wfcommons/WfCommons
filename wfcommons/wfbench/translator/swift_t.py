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
from nis import cat
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
        self.categories_list = []
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
            input = "inputs[]"
            for root_task in self.root_task_names:
                if root_task.startswith(app):
                    input = "inputs"
                    break
            self.script += "@suppress=unused_output\n" \
                f"app (file output) {app} (float percent_cpu, int cpu_work, string out_name, file {input}) " \
                "{\n" \
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
        for task_name in self.root_task_names:
            task = self.tasks[task_name]
            out_count = 0
            for file in task.files:
                if file.link == FileLink.INPUT:
                    self.files_map[file.name] = f"ins[{in_count}]"
                    in_count += 1

                elif file.link == FileLink.OUTPUT:
                    out_count += 1

                if out_count > 1:
                    self.logger.error(
                        "Swift/T does not allow an application to have multiple outputs.")
                    exit(1)

        self.script += f"file ins[] = glob(\"{self.work_dir}/*_input.txt\");\n\n"

        # adding tasks
        for task_name in self.root_task_names:
            self._find_categories_list(task_name)
            # self._add_task(task_name)

        for category in self.categories_list:
            self._add_tasks(category)

        # write script to file
        self._write_output_file(self.script, output_file_name)

    def _find_categories_list(self, task_name: str, parent_task: Optional[str] = None) -> None:
        """"
        Find list of task categories ordered by task dependencies.

        :param task_name: name of the task
        :type task_name: str
        :param parent_task: name of the parent task
        :type parent_task: Optional[str]
        """
        # check dependencies
        for parent in self._find_parents(task_name):
            if parent not in self.parsed_tasks:
                return

        self.parsed_tasks.append(task_name)
        category = self.tasks_map[task_name]
        if category not in self.categories_list:
            self.categories_list.append(category)

        # find children
        children = self._find_children(task_name)

        for child_task_name in children:
            self._find_categories_list(child_task_name)

    def _add_tasks(self, category: str) -> None:
        """
        Add all tasks for a specific category.

        :param category: category name
        :type category: str
        """
        num_tasks = 0
        defined = False
        self.script += f"file {category}_out[];\n"

        for task_name in self.tasks:
            task = self.tasks[task_name]

            if task.category == category:

                # in/output files
                input_files = []
                prefix = ""

                for file in task.files:
                    if file.link == FileLink.OUTPUT:
                        out_file = file.name
                        file_size = file.size
                    elif file.link == FileLink.INPUT:
                        input_files.append(self.files_map[file.name])
                        if not prefix:
                           prefix = self.files_map[file.name].split("_out")[0]
                        elif self.files_map[file.name].split("_out")[0] != prefix:
                           prefix = "Diff"

                # arguments
                args = ", ".join([a.split()[1] for a in task.args[1:3]])
                # args += f", json_objectify(printf(\"'{category}_%i_output.txt': {file_size}\", i))"
                args += f", printf(\"{{'{category}_%i_output.txt': {file_size}}}\", i)"
                if len(input_files) > 0:
                    if prefix.startswith("ins["):
                        args += ", ins[i]"
                    elif len(self._find_parents(task.name)) == 1:
                        args += f", {prefix}_out"
                    else:
                        if not defined:
                            self.script += f"file {category}_in[][];\n"
                            defined = True

                        f_i = 0
                        for f in input_files:
                            self.script += f"{category}_in[{num_tasks}][{f_i}] = {f}; "
                            f_i += 1
                        self.script += "\n"
                        args += f", {category}_in[i]"

                self.files_map[out_file] = f"{category}_out[{num_tasks}]"

                num_tasks += 1

        self.script += f"foreach i in [0:{num_tasks}] {{\n" \
            f"  {category}_out[i] = {category}({args});\n" \
            "}\n\n"

    # def _add_task(self, task_name: str, parent_task: Optional[str] = None) -> None:
    #     """
    #     Add a task and its dependencies to the workflow.

    #     :param task_name: name of the task
    #     :type task_name: str
    #     :param parent_task: name of the parent task
    #     :type parent_task: Optional[str]
    #     """
    #     # check dependencies
    #     for parent in self._find_parents(task_name):
    #         if parent not in self.parsed_tasks:
    #             return

    #     if task_name not in self.parsed_tasks:
    #         task = self.tasks[task_name]

    #         # in/output files
    #         input_files = []
    #         for file in task.files:
    #             if file.link == FileLink.OUTPUT:
    #                 out_file = file.name
    #             elif file.link == FileLink.INPUT:
    #                 input_files.append(self.files_map[file.name])

    #         # arguments
    #         args = ", ".join([a.split()[1] for a in task.args[1:3]])
    #         output_name = task.args[3].replace(
    #             "--out ", "").replace("{", "json_objectify(\"").replace("}", "\")")
    #         args += f", {output_name}"
    #         if len(input_files) > 0:
    #             self.script += f"file in_{self.out_counter}[];"
    #             f_i = 0
    #             for f in input_files:
    #                 self.script += f" in_{self.out_counter}[{f_i}] = {f};"
    #                 f_i += 1
    #             self.script += "\n"
    #             args += f", in_{self.out_counter}"

    #         self.script += f"file out_{self.out_counter} <\"{self.work_dir}/{out_file}\"> = {self.tasks_map[task_name]}({args});\n"
    #         self.files_map[out_file] = f"out_{self.out_counter}"
    #         self.out_counter += 1

    #         self.parsed_tasks.append(task_name)

    #         # find children
    #         children = self._find_children(task_name)

    #         for child_task_name in children:
    #             self._add_task(child_task_name)
