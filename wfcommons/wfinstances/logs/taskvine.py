#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2026 The WfCommons Team.
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
from ...common.file import File
from ...common.machine import Machine
from ...common.task import Task, TaskType
from ...common.workflow import Workflow

####
## Example logs
##
## https://github.com/cooperative-computing-lab/taskvine-example-logs?tab=readme-ov-file
###

class TaskVineLogsParser(LogsParser):
    """
    Parse TaskVine logs to generate workflow instance. This parser has some limitations in that
    it may miss task input/output items due to them not being files or URLs. More importantly,
    because in TaskVine different tasks can have the same file names as input/output but those file names
    actually my correspond to different data sources, the parser will assume that these tasks have
    the exact same input/output files.  There is likely a way to address this, but it hasn't been done yet.
    For instance, the Gutenberg TaskVine example isn't parsed correctly by this parser due to the above feature.

    :param vine_logs_dir: TaskVine's vine-logs directory
    :type vine_logs_dir: pathlib.Path
    :param filenames_to_ignore: TaskVine sometimes considers that executables and package files
                                are input to tasks. This argument is the list of names of files that should be
                                ignored in the reconstructed instances, which typically do not include such
                                files at task input. For instance, if reconstructing a workflow from an execution
                                of a WfBench-generated benchmark, one could pass ["wfbench", "cpu-benchmark", "stress-ng"]
    :type filenames_to_ignore: List[str]
    :param description: Workflow instance description.
    :type description: Optional[str]
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Optional[Logger]
    """
    def __init__(self,
                 vine_logs_dir: pathlib.Path,
                 filenames_to_ignore: Optional[List[str]] = None,
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the makeflow log parser."""
        super().__init__('TaskVine', 'http://https://ccl.cse.nd.edu/software/taskvine/', description, logger)

        # Sanity check
        if not vine_logs_dir.is_dir():
            raise OSError(f'The provided path does not exist or is not a folder: {vine_logs_dir}')

        debug_file:  pathlib.Path = vine_logs_dir / "debug"
        if not debug_file.is_file():
            raise OSError(f'Cannot find file: {debug_file}')
        taskgraph_file: pathlib.Path = vine_logs_dir / "taskgraph"
        if not taskgraph_file.is_file():
            raise OSError(f'Cannot find file: {taskgraph_file}')
        transactions_file: pathlib.Path = vine_logs_dir / "transactions"
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
        self.known_task_ids = []

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
        # sys.stderr.write(str(self.task_command_lines.keys()) + "\n")

        # At this point, the ONLY tasks we care about are the ones for which we have a command-line
        self.known_task_ids = sorted(self.task_command_lines.keys())

        # Construct file map
        self._construct_file_map()
        # sys.stderr.write("FILEMAP: " + str(self.files_map) + "\n")
        for file_key in self.files_map.keys():
            if not "size" in self.files_map[file_key]:
                sys.stderr.write(f"Warning: Could not determine size for file with key {file_key}: assuming zero bytes.\n")
                self.files_map[file_key]["size"] = 0
        sys.stderr.write(f"Identified {len(self.files_map)} valid files\n")

        # Construct the task runtimes
        self._construct_task_runtimes()
        # sys.stderr.write("TASK RUN TIMES: " + str(self.task_runtimes) + "\n")

        # Check whether every known task has a runtime, and if not forget it :(
        to_remove = []
        for task_id in self.known_task_ids:
            if task_id not in self.task_runtimes.keys():
                sys.stderr.write(f"Warning: Ignoring task {task_id} because runtime could not be determined.\n")
                to_remove.append(task_id)
        for victim in to_remove:
            self.known_task_ids.remove(victim)

        sys.stderr.write(f"Identified {len(self.known_task_ids)} valid tasks\n")

        # Construct the input and output file for each task
        self._construct_task_input_output_files()

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
                    # May not be full-proof in case of commands like "export A=b; executable ..." but
                    # may help.....
                    self.filenames_to_ignore.add(executable)
                previous_line = line


    def _construct_file_map(self) -> None:
        filename_to_key_map = {}

        with open(self.taskgraph_file) as f:
            for line in f:
                if line[0] == "#":
                    continue
                
                if line.startswith("FILE"):
                    # FILE file-meta-df1e8b0d0e056c4aedb917abe198a2ff "taskvine_poncho.tar.gz" 0
                    [ignore, file_key, filename, ignore] = line.split()[:4]
                    if len(filename) > 1 and filename[0] == "\"" and filename[-1] == "\"":
                        filename = filename[1:-1]
                    if self.filenames_to_ignore and any(
                        ignore_string in filename for ignore_string in self.filenames_to_ignore
                    ):
                        continue
                    if file_key not in self.files_map:
                        self.files_map[file_key] = {"filename": filename}
        
        with open(self.transactions_file) as f:
            for line in f:
                if line[0] == "#":
                    continue
                if "TRANSFER INPUT" in line:
                    #1769907492293772 245 WORKER worker-8ac3adb4bf1e86b025d0df194b115b8c TRANSFER INPUT file-meta-04276d901d8d096cf981f7ab55f6d1a5 147667979 72409 1769907492221291
                    [ignore, ignore, ignore, ignore, ignore, ignore, file_key, file_size_in_mb] = line.split()[0:8]
                    if file_key not in self.files_map:
                        continue
                    self.files_map[file_key]["size"] = int(file_size_in_mb)
                elif "TRANSFER OUTPUT" in line:
                    #1769907502293773 245 WORKER worker-8ac3adb4bf1e86b025d0df194b115b8c TRANSFER OUTPUT file-meta-04276d901d8d096cf981f7ab55f6d1a5 147667979 72409 1769907502221292
                    [ignore, ignore, ignore, ignore, ignore, ignore, file_key, file_size_in_mb] = line.split()[0:8]
                    if file_key not in self.files_map:
                        continue
                    self.files_map[file_key]["size"] = int(file_size_in_mb)

    def _construct_task_runtimes(self) -> None:
        task_start_times = {}
        task_end_times = {}

        # This method consists only of

        with open(self.transactions_file) as f:
            for line in f:
                if line[0] == "#":
                    continue
                if "RUNNING" in line:
                    [start_date, ignore, ignore, task_index] = line.split()[0:4]
                    if int(task_index) in self.known_task_ids:
                        task_start_times[int(task_index)] = int(start_date)
                elif "DONE" in line:
                    [end_date, ignore, ignore, task_index] = line.split()[0:4]
                    if int(task_index) in self.known_task_ids:
                        task_end_times[int(task_index)] = int(end_date)

        for task_index in task_start_times:
            if task_index in task_end_times:
                self.task_runtimes[task_index] = (
                    float(task_end_times[task_index] - task_start_times[task_index]) / 1_000_000.0)

    def _construct_task_input_output_files(self) -> None:
        # Initialize all entries
        for task_id in self.known_task_ids:
            self.task_input_files[task_id] = []
            self.task_output_files[task_id] = []

        with open(self.taskgraph_file) as f:
            for line in f:               
                if line.startswith("TASK"):
                    # TASK T23 "__vine_env_task-rnd-twtxpejwzsyiebf/bin/run_in_env" INPUTS task-rnd-twtxpejwzsyiebf file-meta-d7504c061a7afd9401c612b4ac7d6be6 file-meta-baab4e4516c4d93a8fcdcbba1a680af7 file-meta-693cb61fedd032b4ddec444b8cce6c89 file-rnd-fnrudlxsaqmlpqq OUTPUTS file-rnd-pdtdqayfmmxyxyp
                    parts = line.split()
                    task_key = parts[1]
                    if not task_key.startswith("T"):
                        continue
                    task_id = int(task_key[1:])  # Remove the T
                    if task_id not in self.known_task_ids:
                        continue
                    input_section = False
                    output_section = False
                    for part in parts[2:]:
                        if part == "INPUTS":
                            input_section = True
                            output_section = False
                            continue
                        elif part == "OUTPUTS":
                            input_section = False
                            output_section = True
                            continue
                        else:
                            if input_section:
                                file_key = part
                                if file_key not in self.files_map:
                                    continue
                                self.task_input_files[task_id].append(self.files_map[file_key]["filename"])
                            elif output_section:
                                file_key = part
                                if file_key not in self.files_map:
                                    continue
                                self.task_output_files[task_id].append(self.files_map[file_key]["filename"])

    def _construct_workflow(self) -> None:
        # Create files and put them in a map
        file_object_map = {}
        for file_key in self.files_map:
            filename = self.files_map[file_key]["filename"]
            file_size = self.files_map[file_key]["size"]
            file_object_map[filename] = File(file_id=filename,
                                    size=file_size,
                                    logger=self.logger)

        # Create all tasks
        task_map = {}
        for task_id in self.known_task_ids:
            task_name = "Task_%d" % task_id
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
