#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib

from logging import Logger
from typing import Optional

from .abstract_translator import Translator
from ...common.file import FileLink


class SwiftTTranslator(Translator):
    """
    A WfFormat parser for creating Swift/T workflow applications.

    :param workflow_json_file_path: Path to the workflow benchmark JSON instance.
    :type workflow_json_file_path: pathlib.Path
    :param work_dir: Path to the workflow working directory.
    :type work_dir: pathlib.Path
    :param stress_path: Path to the stress-ng command.
    :type stress_path: pathlib.Path
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 workflow_json_file_path: pathlib.Path,
                 work_dir: pathlib.Path,
                 stress_path: pathlib.Path = pathlib.Path("stress-ng"),
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow_json_file_path, logger)

        self.work_dir = work_dir
        self.stress_path = stress_path
        self.categories_list = []
        self.parsed_tasks = []
        self.files_map = {}
        self.tasks_map = {}
        self.cmd_counter = 1
        self.script = "import files;\nimport io;\nimport python;\nimport string;\nimport unix;\n\n"

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
        self.script += "string command = \n" \
                "\"\"\"\n" \
                "import json\n" \
                "import os\n" \
                "import pathlib\n" \
                "import socket\n" \
                "import subprocess\n" \
                "\n" \
                f"this_dir = pathlib.Path(\"{self.work_dir}\").absolute()\n" \
                "\n" \
                "print(f\"[WfBench] Starting Benchmark on {socket.gethostname()}\")\n" \
                "\n" \
                "print(\"[WfBench] Starting IO Read Benchmark...\")\n" \
                "files_list = \"%s\"\n" \
                "if \"__\" not in files_list:\n" \
                "    with open(this_dir.joinpath(files_list), \"rb\") as fp:\n" \
                "        print(f\"[WfBench]   Reading '{files_list}'\")\n" \
                "        fp.readlines()\n" \
                "else:\n" \
                "    files = files_list.split(\", \")\n" \
                "    for file in files:\n" \
                "        counter = 0\n" \
                "        fd = file.split(\"__\")\n" \
                "        for f in this_dir.glob(f\"{fd[0]}_*_output.txt\"):\n" \
                "            if counter >= int(fd[1]):\n" \
                "                break\n" \
                "            with open(f, \"rb\") as fp:\n" \
                "                print(f\"[WfBench]   Reading '{f}'\")\n" \
                "                fp.readlines()\n" \
                "            counter += 1\n" \
                "print(\"[WfBench] Completed IO Read Benchmark\")\n" \
                "\n" \
                "print(\"[WfBench] Starting CPU and Memory Benchmarks...\")\n" \
                "cpu_threads=int(10 * %f)\n" \
                "mem_threads=10 - cpu_threads\n" \
                "cpu_work=int(%i)\n" \
                "total_mem_bytes = 0.05\n" \
                "cpu_work_per_thread = int(cpu_work / cpu_threads)\n" \
                "\n" \
                "cpu_procs = []\n" \
                "cpu_prog = [\n" \
                "    f\"{this_dir.joinpath('cpu-benchmark')}\", f\"{cpu_work_per_thread}\"]\n" \
                f"mem_prog = [\"{self.stress_path}\", \"--vm\", f\"{{mem_threads}}\",\n" \
                "            \"--vm-bytes\", f\"{total_mem_bytes}%%\", \"--vm-keep\"]\n" \
                "\n" \
                "for i in range(cpu_threads):\n" \
                "    cpu_proc = subprocess.Popen(cpu_prog)\n" \
                "    cpu_procs.append(cpu_proc)\n" \
                "\n" \
                "mem_proc = subprocess.Popen(mem_prog)\n" \
                "\n" \
                "for proc in cpu_procs:\n" \
                "    proc.wait()\n" \
                "mem_kill = subprocess.Popen([\"killall\", \"stress-ng\"])\n" \
                "mem_kill.wait()\n" \
                "print(\"[WfBench] Completed CPU and Memory Benchmarks\")\n" \
                "\n" \
                "print(f\"[WfBench] Writing output file\")\n" \
                "with open(this_dir.joinpath(\"%s\"), \"wb\") as fp:\n" \
                "    fp.write(os.urandom(int(%i)))\n" \
                "\n" \
                "print(\"[WfBench] Benchmark completed!\")\n" \
                "dep = \"%s\"\n" \
                "\"\"\";\n\n"

        # defining input files
        in_count = 0
        self.script += f"string root_in_files[];\n"

        for task_name in self.root_task_names:
            task = self.tasks[task_name]
            out_count = 0
            for file in task.files:
                if file.link == FileLink.INPUT:
                    self.script += f"root_in_files[{in_count}] = \"{file.name}\";\n"
                    self.files_map[file.name] = f"ins[{in_count}]"
                    in_count += 1

                elif file.link == FileLink.OUTPUT:
                    out_count += 1

                if out_count > 1:
                    self.logger.error(
                        "Swift/T does not allow an application to have multiple outputs.")
                    exit(1)
        
        self.script += "\n"

        # adding tasks
        for task_name in self.root_task_names:
            self._find_categories_list(task_name)

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
        input_files_cat = {}
        parsed_input_files = []
        # defined = False
        self.script += f"string {category}_out[];\n"

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
                        cat_prefix = self.files_map[file.name].split("_out")[0]
                        if file.name not in parsed_input_files:
                            input_files_cat.setdefault(cat_prefix, 0)
                            input_files_cat[cat_prefix] += 1
                            parsed_input_files.append(file.name)
                        input_files.append(self.files_map[file.name])
                        if not prefix:
                           prefix = cat_prefix

                # arguments
                args = ""
                if len(input_files) > 0:
                    # print(f"{task.name}: {len(self._find_parents(task.name))}")
                    if prefix.startswith("ins["):
                        args += "root_in_files[i], "
                    # elif len(self._find_parents(task.name)) == 1:
                    #     args += f"{prefix}_out[0], "
                    else:
                        args += f"{category}_in, "
                    #     if not defined:
                    #         self.script += f"string {category}_in;\n"
                    #         defined = True
                        
                    #     # break input files into several strings
                    #     concat_limit = 50
                    #     start = 0
                    #     end = concat_limit if len(input_files) > concat_limit else len(input_files)
                    #     count = 0
                    #     counter = 0
                    #     while count < len(input_files):
                    #         in_format = ", ".join("%s" for f in input_files[start:end])
                    #         in_args = ", ".join(f for f in input_files[start:end])
                    #         self.script += f"string {category}_inf_{counter} = sprintf(\"{in_format}\", {in_args});\n"

                    #         # in_args = " + \", \" + ".join(f for f in input_files[start:end])
                    #         # self.script += f"string {category}_inf_{counter} = {in_args};\n"
                            
                    #         count = end
                    #         end = end + (concat_limit if len(input_files) > end + concat_limit else len(input_files))
                    #         start = count
                    #         counter += 1
                        
                    #     # in_args = " + \", \" + ".join(f for f in input_files)
                    #     # self.script += f"string {category}_inf = {in_args};\n"
                        
                    #     args += f"{category}_in, "
                    #     in_args = " + \", \" + ".join(f"{category}_inf_{c}" for c in range(0, counter))
                    #     self.script += f"{category}_inf = {in_args};\n"
                    #     self.script += f"{category}_in[{num_tasks}] = {category}_inf;\n"

                args += ", ".join([a.split()[1] for a in task.args[1:3]])
                args += f", of, {file_size}"

                self.files_map[out_file] = f"{category}_out[{num_tasks}]"
                num_tasks += 1

        cats = " + \", \" + ".join(f"{cat}_out" for cat in input_files_cat)
        cats = f"repr({cats})"
        in_str = ", ".join(f"{k}__{v}" for k, v in input_files_cat.items())
        if "ins[" in cats:
            cats = "\"\""
            in_str = ""
        self.script += f"string dep_{self.cmd_counter} = {cats};\n"
        args += f", dep_{self.cmd_counter}"
        self.script += f"string {category}_in = \"{in_str}\";\n"

        if num_tasks > 1:
            self.script += f"foreach i in [0:{num_tasks - 1}] {{\n" \
                f"  string of = sprintf(\"{category}_%i_output.txt\", i);\n" \
                f"  string cmd_{self.cmd_counter} = sprintf(command, {args});\n" \
                f"  string co_{self.cmd_counter} = python(cmd_{self.cmd_counter});\n" \
                f"  string of_{self.cmd_counter} = sprintf(\"'{category}_%i_output.txt%s'\", i, co_{self.cmd_counter});\n" \
                f"  {category}_out[i] = of_{self.cmd_counter};\n" \
                "}\n\n"
        else:
            args = args.replace(
                ", of", f", \"{category}_0_output.txt\"").replace("[i]", "[0]")
            self.script += f"string cmd_{self.cmd_counter} = sprintf(command, {args});\n" \
                f"string co_{self.cmd_counter} = python(cmd_{self.cmd_counter});\n" \
                f"string of_{self.cmd_counter} = sprintf(\"'{category}_0_output.txt%s'\", co_{self.cmd_counter});\n" \
                f"{category}_out[0] = of_{self.cmd_counter};\n\n"
        
        self.cmd_counter += 1
