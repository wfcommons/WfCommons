#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib

from logging import Logger
from typing import Optional, Union

from .abstract_translator import Translator
from ...common import Workflow


class SwiftTTranslator(Translator):
    """
    A WfFormat parser for creating Swift/T workflow applications.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path]
    :param stress_path: Path to the stress-ng command.
    :type stress_path: pathlib.Path
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 stress_path: pathlib.Path = pathlib.Path("stress-ng"),
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow, logger)

        self.stress_path = stress_path
        self.categories_list = []
        self.categories_input = {}
        self.parsed_tasks = set()
        self.files_map = {}
        self.tasks_map = {}
        self.cmd_counter = 1
        self.out_files = set()

        # find applications
        self.apps = set()
        for task in self.tasks.values():
            self.tasks_map[task.task_id] = task.name
            self.apps.add(task.name)
            
            out_count = 0
            for file in task.output_files:
                self.files_map[file.file_id] = f"{task.name}__out"
                out_count += 1

            if out_count > 1:
                self.logger.error(
                    "Swift/T does not allow an application to have multiple outputs.")
                exit(1)

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into a Swift/T workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        self.logger.info("Translating workflow into Swift/T")

        # defining input files
        self.logger.debug("Defining input files")
        in_count = 0
        self.output_folder = output_folder
        self.cpu_benchmark = output_folder.joinpath("./bin/cpu-benchmark").absolute()
        self.script = f"string fs = sprintf(flowcept_start, \"{self.workflow.workflow_id}\");\nstring fss = python_persist(fs);\n\n" if self.workflow.workflow_id else ""
        self.script += "string root_in_files[];\n"

        for task_name in self.root_task_names:
            task = self.tasks[task_name]
            for file in task.input_files:
                if task.name not in self.categories_input.keys():
                    self.categories_input[task.name] = in_count
                    in_file = output_folder.joinpath(f"./data/{file.file_id}").absolute()
                    self.script += f"root_in_files[{in_count}] = \"{in_file}\";\n"
                    in_count += 1
                self.files_map[file.file_id] = f"ins[{in_count}]"
        
        self.script += "\n"

        # adding tasks
        self.logger.info("Finding categories list")
        for task_name in self.root_task_names:
            self._find_categories_list(task_name)

        self.logger.info("Adding tasks")
        for category in self.categories_list:
            self._add_tasks(category)

        # flowcept end
        if self.workflow.workflow_id:
            out_files = ", ".join(f"'{item}'" for item in self.out_files)
            self.script += f"string fc = sprintf(flowcept, \"{self.workflow.workflow_id}\", \"{self.workflow.name}\", \"{out_files}\");\n" \
                            "python_persist(fc);\n"

        run_workflow_code = self._merge_codelines("templates/swift_t/workflow.swift", self.script)

        # write benchmark files
        output_folder.mkdir(parents=True)
        with open(output_folder.joinpath("workflow.swift"), "w") as fp:
            fp.write(run_workflow_code)

        # additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)

    def _find_categories_list(self, task_name: str, parent_task: Optional[str] = None) -> None:
        """"
        Find list of task categories ordered by task dependencies.

        :param task_name: name of the task
        :type task_name: str
        :param parent_task: name of the parent task
        :type parent_task: Optional[str]
        """
        task = self.tasks[task_name]
        if task.name in self.categories_list:
            return

        # check dependencies
        for parent in self.task_parents[task_name]:
            parent_task = self.tasks[parent]
            if parent_task.name not in self.categories_list:
                return

        self.parsed_tasks.add(task_name)
        category = self.tasks_map[task_name]
        if category not in self.categories_list:
            self.categories_list.append(category)

        # find children
        children = self.task_children[task_name]

        for child_task_name in children:
            self._find_categories_list(child_task_name)

    def _add_tasks(self, category: str) -> None:
        """
        Add all tasks for a specific category.

        :param category: category name
        :type category: str
        """
        num_tasks = 0
        num_children = 0
        input_files_cat = {}
        parsed_input_files = set()
        self.script += f"int {category}__out[];\n"

        for task in self.tasks.values():
            if task.name == category:
                # in/output files
                num_children = len(self.task_children[task.task_id])
                input_files = set()
                prefix = ""

                for file in task.output_files:
                    file_size = file.size
                
                for file in task.input_files:
                    cat_prefix = self.files_map[file.file_id].split("__out")[0]
                    if file.file_id not in parsed_input_files:
                        input_files_cat.setdefault(cat_prefix, 0)
                        input_files_cat[cat_prefix] += 1
                        parsed_input_files.add(file.file_id)
                    input_files.add(self.files_map[file.file_id])
                    if not prefix:
                        prefix = cat_prefix

                # arguments
                if num_tasks == 0:
                    args = ""
                    if len(input_files) > 0:
                        if prefix.startswith("ins["):
                            args += f"root_in_files[{self.categories_input[task.name]}], "
                        else:
                            args += f"{category}_in, "

                    percent_cpu = "0.1"
                    cpu_work = "0"
                    gpu_work = "0"
                    for arg in task.args:
                        if arg.startswith("--percent-cpu"):
                            percent_cpu = arg.split()[1]
                        elif arg.startswith("--cpu-work"):
                            cpu_work = arg.split()[1]
                        elif arg.startswith("--gpu-work"):
                            gpu_work = arg.split()[1]

                    args += ", ".join([gpu_work, cpu_work, percent_cpu])
                    args += f", of, {file_size}"

                num_tasks += 1

        cats = " + ".join(f"{k}__out[{v - 1}]" for k, v in input_files_cat.items())
        in_str = ", ".join(f"{k}_{v - 1}_output.txt" for k, v in input_files_cat.items())
        if "ins[" in cats:
            cats = "0"
            in_str = ""
        self.script += f"int dep_{self.cmd_counter} = {cats};\n"
        args += f", dep_{self.cmd_counter}"
        args += f", \"{self.workflow.workflow_id}\", fss" if self.workflow.workflow_id else ", \"\""
        self.script += f"string {category}_in = \"{self.output_folder.absolute()}/data/{in_str}\";\n"

        if num_tasks > 1:
            self.script += f"foreach i in [0:{num_tasks - 1}] {{\n" \
                f"  string of = sprintf(\"{self.output_folder.absolute()}/data/{category}_%i_output.txt\", i);\n" \
                f"  string task_id = \"{category}_\" + i;\n" \
                f"  string cmd_{self.cmd_counter} = sprintf(command, \"{self.cpu_benchmark}\", task_id, {args});\n" \
                f"  string co_{self.cmd_counter} = python_persist(cmd_{self.cmd_counter});\n" \
                f"  string of_{self.cmd_counter} = sprintf(\"0%s\", co_{self.cmd_counter});\n" \
                f"  {category}__out[i] = string2int(of_{self.cmd_counter});\n" \
                "}\n\n"
            if not num_children:
                for i in range(num_tasks):
                    self.out_files.add(f"{self.output_folder.absolute()}/data/{category}_{i}_output.txt")
        else:
            out_file = f"{self.output_folder.absolute()}/data/{category}_0_output.txt"
            if not num_children:
                self.out_files.add(out_file)
            args = args.replace(
                ", of", f", \"{out_file}\"").replace("[i]", "[0]")
            self.script += f"string cmd_{self.cmd_counter} = sprintf(command, \"{self.cpu_benchmark}\", \"{category}_{self.cmd_counter}\", {args});\n" \
                f"string co_{self.cmd_counter} = python_persist(cmd_{self.cmd_counter});\n" \
                f"string of_{self.cmd_counter} = sprintf(\"0%s\", co_{self.cmd_counter});\n" \
                f"{category}__out[0] = string2int(of_{self.cmd_counter});\n\n"
        
        self.cmd_counter += 1
