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
        self.script = "import files;\nimport io;\nimport unix;\n\n"

        # find applications
        self.apps = {}
        for task in self.tasks.values():
            inputs_count = 0
            outputs_count = 0
            for file in task.files:
                if file.link == FileLink.INPUT:
                    inputs_count += 1
                else:
                    outputs_count += 1

            self.tasks_map[task.name] = f"{task.category}_{inputs_count}_{outputs_count}"

            if task.category not in self.apps:
                self.apps[task.category] = {}

            if f"{inputs_count}_{outputs_count}" not in self.apps[task.category]:
                self.apps[task.category][f"{inputs_count}_{outputs_count}"] = {
                    "name": task.category,
                    "inputs": inputs_count,
                    "outputs": outputs_count
                }

    def translate(self, output_file_name: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into a Swift/T workflow application.

        :param output_file_name: The name of the output file (e.g., workflow.swift).
        :type output_file_name: pathlib.Path
        """
        # creating apps
        for a in self.apps:
            for io in self.apps[a]:
                app = self.apps[a][io]
                outputs = ", ".join(
                    [f"file out_{i + 1}" for i in range(0, app["outputs"])])
                outputs_o = outputs.replace("file ", "").replace(",", "")
                inputs = ", file inputs[]" if app["inputs"] > 0 else ""
                inputs_o = "\\\n  inputs\n" if app["inputs"] > 0 else "\n"

                self.script += f"app ({outputs}) {app['name']}_{io} (string path_lock, string path_cores, float percent_cpu, int cpu_work, int file_size{inputs}) "
                self.script += "{\n" \
                    "  \"/sw/summit/python/3.8/anaconda3/2020.07-rhel8/bin/python3\" \\\n" \
                    f"  \"{self.work_dir}/wfperf_benchmark.py\" \\\n" \
                    f"  \"{app['name']}_{io}\" \\\n" \
                    "  \"--path-lock\" path_lock \\\n" \
                    "  \"--path-cores\" path_cores \\\n" \
                    "  \"--percent-cpu\" percent_cpu \\\n" \
                    "  \"--cpu-work\" cpu_work \\\n" \
                    "  \"--data\" \\\n" \
                    "  \"--outputs-file-size\" file_size \\\n" \
                    f"  \"--out\" {outputs_o} {inputs_o}" \
                    "}\n\n"

        # defining input files
        in_count = 0
        for task_name in self.parent_task_names:
            task = self.tasks[task_name]
            for file in task.files:
                if file.link == FileLink.INPUT:
                    self.files_map[file.name] = f"ins[{in_count}]"
                    in_count += 1
        self.script += f"file ins[] = glob(\"{self.work_dir}/sys_input_*.txt\");\n"
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
            args = ", ".join([f"\"{a.split('=')[1]}\"" for a in task.args[1:3]])
            args += f', {", ".join([a.split("=")[1] for a in task.args[3:5]])}'
            args += f", {task.args[6].split('=')[1]}"
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
