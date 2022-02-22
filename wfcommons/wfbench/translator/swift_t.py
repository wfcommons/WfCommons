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
        self.logger.info("Translating workflow into Swift/T")
        self.script += "string command = \n" \
                "\"\"\"\n" \
                "import json\n" \
                "import os\n" \
                "import pathlib\n" \
                "import socket\n" \
                "import subprocess\n" \
                "import time\n" \
                "\n" \
                f"this_dir = pathlib.Path(\"{self.work_dir}\").absolute()\n" \
                "\n" \
                "task_name = \"%s\"\n" \
                "print(f\"[WfBench] [{task_name}] Starting Benchmark on {socket.gethostname()}\")\n" \
                "\n" \
                "print(f\"[WfBench] [{task_name}] Starting IO Read Benchmark...\")\n" \
                "files_list = \"%s\"\n" \
                "if \"__\" not in files_list:\n" \
                "    with open(this_dir.joinpath(files_list), \"rb\") as fp:\n" \
                "        start = time.perf_counter()\n" \
                "        print(f\"[WfBench]   Reading '{files_list}'\")\n" \
                "        fp.readlines()\n" \
                "        end = time.perf_counter()\n" \
                "        print(f\"[WfBench] [{task_name}] Metrics (read) [time,size]: {end - start},{this_dir.joinpath(files_list).stat().st_size}\")\n" \
                "else:\n" \
                "    files = files_list.split(\", \")\n" \
                "    for file in files:\n" \
                "        counter = 0\n" \
                "        fd = file.split(\"__\")\n" \
                "        start = time.perf_counter()\n" \
                "        file_size = 0\n" \
                "        for f in this_dir.glob(f\"{fd[0]}_*_output.txt\"):\n" \
                "            if counter >= int(fd[1]):\n" \
                "                break\n" \
                "            file_size += os.stat(f).st_size\n" \
                "            with open(f, \"rb\") as fp:\n" \
                "                print(f\"[WfBench]   Reading '{f}'\")\n" \
                "                fp.readlines()\n" \
                "            counter += 1\n" \
                "        end = time.perf_counter()\n" \
                "        print(f\"[WfBench] [{task_name}] Metrics (read) [time,size]: {end - start},{file_size}\")\n" \
                "print(f\"[WfBench] [{task_name}] Completed IO Read Benchmark\")\n" \
                "\n" \
                "print(f\"[WfBench] [{task_name}] Starting CPU and Memory Benchmarks...\")\n" \
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
                "start = time.perf_counter()\n" \
                "for i in range(cpu_threads):\n" \
                "    cpu_proc = subprocess.Popen(cpu_prog)\n" \
                "    cpu_procs.append(cpu_proc)\n" \
                "\n" \
                "mem_proc = subprocess.Popen(mem_prog, stderr=subprocess.DEVNULL)\n" \
                "\n" \
                "for proc in cpu_procs:\n" \
                "    proc.wait()\n" \
                "mem_kill = subprocess.Popen([\"killall\", \"stress-ng\"])\n" \
                "mem_kill.wait()\n" \
                "end = time.perf_counter()\n" \
                "print(f\"[WfBench] [{task_name}] Metrics (compute) [time,work]: {end - start},{cpu_work}\")\n" \
                "print(f\"[WfBench] [{task_name}] Completed CPU and Memory Benchmarks\")\n" \
                "\n" \
                "print(f\"[WfBench] [{task_name}] Writing output file\")\n" \
                "start = time.perf_counter()\n" \
                "with open(this_dir.joinpath(\"%s\"), \"wb\") as fp:\n" \
                "    file_size = int(%i)\n" \
                "    fp.write(os.urandom(file_size))\n" \
                "end = time.perf_counter()\n" \
                "print(f\"[WfBench] [{task_name}] Metrics (write) [time,size]: {end - start},{file_size}\")\n" \
                "\n" \
                "print(f\"[WfBench] [{task_name}] Benchmark completed!\")\n" \
                "dep = \"%s\"\n" \
                "\"\"\";\n\n"

        # defining input files
        self.logger.debug("Defining input files")
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
        self.logger.info("Finding categories list")
        for task_name in self.root_task_names:
            self._find_categories_list(task_name)

        self.logger.info("Adding tasks")
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
        task = self.tasks[task_name]
        if task.category in self.categories_list:
            return

        # check dependencies
        for parent in self.task_parents[task_name]:
            parent_task = self.tasks[parent]
            if parent_task.category not in self.categories_list:
                return

        self.parsed_tasks.append(task_name)
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
        input_files_cat = {}
        parsed_input_files = []
        self.script += f"string {category}__out[];\n"

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
                        cat_prefix = self.files_map[file.name].split("__out")[0]
                        if file.name not in parsed_input_files:
                            input_files_cat.setdefault(cat_prefix, 0)
                            input_files_cat[cat_prefix] += 1
                            parsed_input_files.append(file.name)
                        input_files.append(self.files_map[file.name])
                        if not prefix:
                           prefix = cat_prefix

                # arguments
                if num_tasks == 0:
                    args = ""
                    if len(input_files) > 0:
                        if prefix.startswith("ins["):
                            args += "root_in_files[i], "
                        else:
                            args += f"{category}_in, "

                    args += ", ".join([a.split()[1] for a in task.args[1:3]])
                    args += f", of, {file_size}"

                self.files_map[out_file] = f"{category}__out[{num_tasks}]"
                num_tasks += 1

        cats = " + \", \" + ".join(f"repr({cat}__out)" for cat in input_files_cat)
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
                f"  string cmd_{self.cmd_counter} = sprintf(command, \"{category}\", {args});\n" \
                f"  string co_{self.cmd_counter} = python(cmd_{self.cmd_counter});\n" \
                f"  string of_{self.cmd_counter} = sprintf(\"'{category}_%i_output.txt%s'\", i, co_{self.cmd_counter});\n" \
                f"  {category}__out[i] = of_{self.cmd_counter};\n" \
                "}\n\n"
        else:
            args = args.replace(
                ", of", f", \"{category}_0_output.txt\"").replace("[i]", "[0]")
            self.script += f"string cmd_{self.cmd_counter} = sprintf(command, \"{category}\", {args});\n" \
                f"string co_{self.cmd_counter} = python(cmd_{self.cmd_counter});\n" \
                f"string of_{self.cmd_counter} = sprintf(\"'{category}_0_output.txt%s'\", co_{self.cmd_counter});\n" \
                f"{category}__out[0] = of_{self.cmd_counter};\n\n"
        
        self.cmd_counter += 1
