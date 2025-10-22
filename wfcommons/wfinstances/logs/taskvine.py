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
                    # May not be full-proof in case of commands like "export A=b; executable ..." but
                    # may help.....
                    self.filenames_to_ignore.add(executable)
                previous_line = line


    def _construct_file_map(self) -> None:

        filename_to_key_map = {}
        # One pass through the debug file to create the initial file key -> filename mapping
        with open(self.debug_file) as f:
            for line in f:
                if "__vine_env_task" in line: # Ignore that weird task/file
                    continue
                if "infile " in line :
                    # 2025/09/09 21:12:48.02 vine_manager[239]vine: tx to dab178765b01 (127.0.0.1:34382): infile file-rnd-fmtpwpiobiumeze blastall_00000016_outfile_0016 0
                    [file_key, filename] = line[line.find("infile ") + len("infile "):].split()[:2]
                elif "outfile " in line and "completed with outfile " not in line and "outfile =" not in line:
                    # 2025/09/30 18:37:19.74 vine_manager[1849324]vine: tx to d64cepyc028.crc.nd.edu (10.32.94.18:47558): outfile temp-rnd-pidiwheippcwbeu fde2b5eb-9713-423a-8bc6-f4f9263ad20b.pkl 0 3
                    [file_key, filename] = line[line.find("outfile ") + len("outfile "):].split()[:2]
                else:
                    continue
                if filename in self.filenames_to_ignore:
                    continue
                # NOTE THAT THE FILENAME MAY NOT BE UNIQUE IN TASKVINE WORKFLOWS, SO
                # WE ADD THE KEY
                self.files_map[file_key] = {"filename": filename + "." + file_key}
                filename_to_key_map[filename] = file_key

        # Pass through the transactions file to get the file sizes
        with open(self.debug_file) as f:
            for line in f:
                if "): file " in line:
                    [file_key, file_size] = line[line.find("): file ") + len("): file "):].split()[0:2]
                else:
                    continue
                if file_key in self.files_map:
                    self.files_map[file_key]["size"] = int(file_size)


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
                if "->" not in line:
                    continue
                if "file-task" in line:  # Ignoring what I think are taskvine internal/specific things
                    continue
                line = line[:-1]
                # print(f"LINE: {line}")
                [source, ignore, destination] = line.split()
                # Remove quotes
                source = source [1:-1]
                destination = destination [1:-2]
                # Remove weird file- prefix
                source = source.replace("--", "-")  # Sometimes there is an unexpected "--"!!
                destination = destination.replace("--", "-")  # Sometimes there is an unexpected "--"!!
                # print(f"source: {source}  destination: {destination}")
                if source.startswith("file-"):
                    source = source[len("file-"):]
                if destination.startswith("file-"):
                    destination = destination[len("file-"):]

                if "task-" in source and "file-" not in source:
                    try:
                        task_id = int(source.split("-")[1])
                    except ValueError as e:
                        raise Exception(f"The source was {source} and the split around '-' failed!")

                    if task_id not in self.task_runtimes:
                        continue
                    file_key = destination
                    if file_key not in self.files_map:
                        continue
                    output_file = self.files_map[file_key]["filename"]
                    self.task_output_files[task_id].append(output_file)
                elif "task" in destination and "file" not in destination:
                    try:
                        task_id = int(destination.split("-")[1])
                    except ValueError as e:
                        raise Exception(f"The destination was {destination} and the split around '-' failed!")
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
            file_object_map[filename] = File(file_id=filename,
                                    size=file_size,
                                    logger=self.logger)

        # Create all tasks
        task_map = {}
        # print(self.task_runtimes[16])
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
