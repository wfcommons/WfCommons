#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import itertools
import math
import pathlib
import sys

from datetime import datetime, timezone
from logging import Logger
from typing import List, Optional

from .abstract_logs_parser import LogsParser
from ...common.file import File, FileLink
from ...common.machine import Machine
from ...common.task import Task, TaskType
from ...common.workflow import Workflow


class TaskVineLogsParser(LogsParser):
    """
    Parse TaskVine logs to generate workflow instance.

    :param vine_run_info_dir: TaskVine's  vine-run-info directory.
    :type vine_run_info_dir: pathlib.Path
    :param filenames_to_ignore: TaskVine considers that executables and package files (e.g., poncho package.tgz)
                                are input to tasks. This argument is the list of names of files that should be
                                ignored in the reconstructed instances, which typically do not include such
                                files at task input.
    :type filenames_to_ignore: List[str]
    :param description: Workflow instance description.
    :type description: Optional[str]
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Optional[Logger]
    """
    def __init__(self,
                 vine_run_info_dir: pathlib.Path,
                 filenames_to_ignore: Optional[List[str]] = None,
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the makeflow log parser."""
        super().__init__('TaskVine', 'http://https://ccl.cse.nd.edu/software/taskvine/', description, logger)

        # Sanity check
        if not vine_run_info_dir.is_dir():
            raise OSError(f'The provided path does not exist or is not a folder: {vine_run_info_dir}')

        debug_file:  pathlib.Path = vine_run_info_dir / "most-recent/vine-logs/debug"
        if not debug_file.is_file():
            raise OSError(f'Cannot find file: {debug_file}')
        taskgraph_file: pathlib.Path = vine_run_info_dir / "most-recent/vine-logs/taskgraph"
        if not taskgraph_file.is_file():
            raise OSError(f'Cannot find file: {taskgraph_file}')
        transactions_file: pathlib.Path = vine_run_info_dir / "most-recent/vine-logs/transactions"
        if not transactions_file.is_file():
            raise OSError(f'Cannot find file: {transactions_file}')

        self.debug_file: pathlib.Path = debug_file
        self.taskgraph_file: pathlib.Path = taskgraph_file
        self.transactions_file: pathlib.Path = transactions_file

        self.filenames_to_ignore: set[str] = set(filenames_to_ignore) or set({})

        self.files_map = {}
        self.task_command_lines = {}
        self.task_runtimes = {}
        self.task_input_files = {}
        self.task_output_files = {}

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow instance based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: Optional[str]

        :return: A workflow instance object.
        :rtype: Workflow
        """
        self.workflow_name = workflow_name

        # create base workflow instance object
        self.workflow = Workflow(name=self.workflow_name,
                                 description=self.description,
                                 runtime_system_name=self.wms_name,
                                 runtime_system_url=self.wms_url)

        # Construct the task command-line array
        self._construct_task_command_lines()
        # sys.stderr.write(str(self.task_command_lines) + "\n")

        # Construct file map
        self._construct_file_map()
        # sys.stderr.write("FILEMAP: " + str(self.files_map) + "\n")

        # Construct the task runtimes
        self._construct_task_runtimes()
        # sys.stderr.write("TASK RUN TIMES: " + str(self.task_runtimes) + "\n")

        # Construct the input and output file for each task
        self._construct_task_input_output_files()
        # print("TASK INPUT FILES: " + str(self.task_input_files))
        # print("TASK OUTPUT FILES: " + str(self.task_output_files))

        # Construct the workflow
        self._construct_workflow()

        return self.workflow


    def _construct_task_command_lines(self) -> None:
        with open(self.debug_file) as f:
            for line in f:
                if "state change: READY (1) to RUNNING (2)" in line:
                    [task_index] = line[line.find("Task ") + len("Task "):].split()[0:1]
                    command_line = previous_line[previous_line.find("busy on '") + len("busy on '"):-2]
                    self.task_command_lines[int(task_index)] = command_line
                    executable = command_line.split()[0]
                    self.filenames_to_ignore.add(executable)
                previous_line = line


    def _construct_file_map(self) -> None:

        filename_to_key_map = {}
        # One pass through the debug file to create the initial file key -> filename mapping
        with open(self.debug_file) as f:
            for line in f:
                if "__vine_env_task" in line: # Ignore that weird task/file
                    continue
                # 2025/09/09 21:12:48.02 vine_manager[239]vine: tx to dab178765b01 (127.0.0.1:34382): infile file-rnd-fmtpwpiobiumeze blastall_00000016_outfile_0016 0
                if "infile " in line:
                    [file_key, filename] = line[line.find("infile ") + len("infile "):].split()[0:2]
                # 2025/09/09 21:12:42.12 vine_manager[239]vine: tx to dab178765b01 (127.0.0.1:34382): outfile file-rnd-jpnzrjrjnxxqhej blastall_00000017_outfile_0017 0
                elif "outfile " in line:
                    [file_key, filename] = line[line.find("outfile ") + len("outfile "):].split()[0:2]
                else:
                    continue
                if filename in self.filenames_to_ignore:
                    continue
                self.files_map[file_key] = {"filename": filename}
                filename_to_key_map[filename] = file_key

        # Another pass through the debug file to get the actual file paths
        with open(self.debug_file) as f:
            for line in f:
                # 2025/09/09 21:12:48.01 vine_manager[239]vine: dab178765b01 (127.0.0.1:34382) needs file data/blastall_00000003_outfile_0003 as blastall_00000003_outfile_0003
                if "needs file " in line:
                    [full_path, ignore, filename] = line[line.find("needs file ") + len("needs file "):].split()[0:3]
                    file_key = filename_to_key_map.get(filename)
                # 2025/09/09 21:12:47.92 vine_manager[239]vine: dab178765b01 (127.0.0.1:34382) sending back file-rnd-jajwzwsrtyzbkfs to data/blastall_00000020_outfile_0020
                elif "sending back " in line:
                    [file_key, ignore, full_path] = line[line.find("sending back ") + len("sending back "):].split()[0:3]
                    filename = self.files_map[file_key]["filename"]
                else:
                    continue
                if filename in self.filenames_to_ignore:
                    continue
                self.files_map[file_key]["path"] = full_path

        # Pass through the transactions file to get the file sizes
        with open(self.transactions_file) as f:
            for line in f:
                # 1757452362084671 239 WORKER worker-50dc215f08057f4005f3b65dee08592f TRANSFER OUTPUT file-rnd-wzkjcrgiivvzbci 227273 1327 1757452362083301
                # 1757452358704968 239 WORKER worker-50dc215f08057f4005f3b65dee08592f TRANSFER INPUT file-meta-9b84b334875319e856f72be634aae964 17648 1129 1757452358703835
                if "TRANSFER INPUT " in line:
                    [file_key, file_size] = line[line.find("TRANSFER INPUT ") + len("TRANSFER INPUT "):].split()[0:2]
                elif "TRANSFER OUTPUT " in line:
                    [file_key, file_size] = line[line.find("TRANSFER OUTPUT ") + len("TRANSFER OUTPUT "):].split()[0:2]
                elif "CACHE_UPDATE " in line:
                    [file_key, file_size] = line[line.find("CACHE_UPDATE ") + len("CACHE_UPDATE "):].split()[0:2]
                else:
                    continue
                if file_key in self.files_map:
                    self.files_map[file_key]["size"] = int(file_size)

        # print(str(self.files_map))


    def _construct_task_runtimes(self) -> None:
        task_start_times = {}
        task_end_times = {}

        with open(self.transactions_file) as f:
            for line in f:
                if line[0] == "#":
                    continue
                if "RUNNING" in line:
                    [start_date, ignore, ignore, task_index] = line.split()[0:4]
                    task_start_times[int(task_index)] = int(start_date)
                elif "DONE" in line:
                    [end_date, ignore, ignore, task_index] = line.split()[0:4]
                    task_end_times[int(task_index)] = int(end_date)

        for task_index in task_start_times:
            self.task_runtimes[task_index] = (
                    float(task_end_times[task_index] - task_start_times[task_index]) / 1_000_000.0)

    def _construct_task_input_output_files(self) -> None:

        # Initialize all entries
        for task_id in self.task_runtimes.keys():
            self.task_input_files[task_id] = []
            self.task_output_files[task_id] = []

        with open(self.taskgraph_file) as f:
            for line in f:
                if "->" not in line:
                    continue
                if "file-task" in line:  # Ignoring what I think are taskvine internal/specific things
                    continue
                line = line[:-1]
                print(f"LINE: {line}")
                [source, ignore, destination] = line.split()
                # Remove quotes
                source = source [1:-1]
                destination = destination [1:-2]
                # Remove weird file- prefix
                if source.startswith("file-"):
                    source = source[len("file-"):]
                if destination.startswith("file-"):
                    destination = destination[len("file-"):]

                if "task" in source and "file" not in source:
                    task_id = int(source.split("-")[1])
                    if task_id not in self.task_runtimes:
                        continue
                    file_key = destination
                    if file_key not in self.files_map:
                        continue
                    output_file = self.files_map[file_key]["filename"]
                    self.task_output_files[task_id].append(output_file)
                elif "task" in destination and "file" not in destination:
                    task_id = int(destination.split("-")[1])
                    if task_id not in self.task_runtimes:
                        continue
                    file_key = source
                    if file_key not in self.files_map:
                        continue
                    input_file = self.files_map[file_key]["filename"]
                    self.task_input_files[task_id].append(input_file)
                else:
                    raise ValueError("Error in the taskgraph file")


    def _construct_workflow(self) -> None:
        # Create files and put them in a map
        file_object_map = {}
        for file_key in self.files_map:
            filename = self.files_map[file_key]["filename"]
            file_size = self.files_map[file_key]["size"]
            # file_path = self.files_map[file_key]["path"]
            file_object_map[filename] = File(file_id=filename,
                                    size=file_size,
                                    link=FileLink.INPUT, # TODO: REMOVE THAT ONCE NOLINK HAS BEEN MERGED
                                    logger=self.logger)

        # Create all tasks
        task_map = {}
        for task_id in self.task_runtimes:
            task_name = "Task %d" % task_id
            task = Task(name=task_name,
                        task_id=task_name,
                        task_type=TaskType.COMPUTE,
                        runtime=self.task_runtimes[task_id],
                        program=self.task_command_lines[task_id].split()[0],
                        args=self.task_command_lines[task_id].split()[1:],
                        cores=1,
                        input_files=[file_object_map[filename] for filename in self.task_input_files[task_id]],
                        output_files=[file_object_map[filename] for filename in self.task_output_files[task_id]],
                        logger=self.logger)
            task_map[task_id] = task
            self.workflow.add_task(task)
            # sys.stderr.write(f"Added task {task_name}: {len(self.workflow.tasks)}\n")


        # Adding all edges, which is pretty inefficiently done for now, by looking at all pairs of tasks!
        for task1_id in self.task_runtimes:
            for task2_id in self.task_runtimes:
                if task1_id == task2_id:
                    continue
                task1_output_files = self.task_output_files[task1_id]
                task2_input_files = self.task_input_files[task2_id]
                has_intersection = bool(set(task1_output_files) & set(task2_input_files))
                if has_intersection:
                    self.workflow.add_dependency(task_map[task1_id].name, task_map[task2_id].name)
                    # sys.stderr.write(f"Added dependency {task_map[task1_id].name} -> {task_map[task2_id].name}\n")