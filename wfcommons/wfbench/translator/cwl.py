#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import shutil
import logging
import pathlib
import shlex
from typing import Union, Optional
from collections import defaultdict, deque

from .abstract_translator import Translator
from ...common import Workflow

this_dir = pathlib.Path(__file__).resolve().parent

class CWLTranslator(Translator):
    """
    A WfFormat parser for creating CWL workflow benchmarks.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path],
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """
    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 logger: Optional[logging.Logger] = None) -> None:
        super().__init__(workflow, logger)
        self.cwl_script = ["cwlVersion: v1.2",
                           "class: Workflow",
                           "requirements:",
                           "  MultipleInputFeatureRequirement: {}",
                           "  StepInputExpressionRequirement: {}",
                           "  InlineJavascriptRequirement: {}"]
        self.yml_script = []
        self.parsed_tasks = []
        self.task_level_map = defaultdict(lambda: [])

        queue = deque(self.root_task_names)
        visited = set()
        top_sort = []

        while queue:
            task_name = queue.popleft()

            if task_name not in visited:
                top_sort.append(task_name)
                visited.add(task_name)

            for child in self.task_children[task_name]:
                queue.append(child)

        assert (len(top_sort) == len(self.tasks))

        levels = {task_name: 0 for task_name in top_sort}

        for task_name in top_sort:
            for child in self.task_children[task_name]:
                levels[child] = max(levels[child], levels[task_name] + 1)

        for task_name, level in levels.items():
            self.task_level_map[level].append(task_name)

    def translate(self, output_folder: pathlib.Path) -> None:
        # Create output folder
        output_folder.mkdir(parents=True)

        # Parsing the inputs and outputs of the workflow
        self._parse_inputs_outputs()

        # Parsing the steos
        self._parse_steps()

        # additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)

        # Writing the CWL files to the output folder
        self._write_cwl_files(output_folder)

        return 0

    def _parse_steps(self) -> None:
        steps_folder_source = []
        self.cwl_script.append("\nsteps:")
        # Parsing each steps by Workflow levels
        for level in sorted(self.task_level_map.keys()):
            # Parsing each task within a Workflow level
            for task_name in self.task_level_map[level]:

                # Getting the task object
                task = self.tasks[task_name]

                # Parsing the arguments of the task
                args_array = []
                benchmark_name = False

                for item in task.args:
                    # Split elements that contain both an option and a value
                    if item.startswith("--"):
                        item = item.replace("\'", "\"")
                        item = item.split(" ", 1)
                        args_array.append(item[0])
                        args_array.append(item[1])
                    elif not benchmark_name:
                        args_array.append(item)
                        benchmark_name = True

                output_files = [
                    f.file_id for f in task.output_files]

                # Adding the step to the cwl script

                self.cwl_script.append(f"  {task.task_id}:")
                # TODO: change so that it doesn't only run wfbench programs
                if not task.program.startswith("wfbench"):
                    raise ValueError("Only wfbench programs are supported")
                self.cwl_script.append("    run: clt/wfbench.cwl")

                self.cwl_script.append("    in:")
                if level == 0:
                    self.cwl_script.append(
                        f"      input_files: {task.task_id}_input")
                else:
                    self.cwl_script.append(
                        "      input_files:")
                    self.cwl_script.append(
                        "        linkMerge: merge_flattened")
                    self.cwl_script.append(
                        "        source:")
                    for parent in self.task_parents[task_name]:
                        self.cwl_script.append(
                            f"          - {parent}/output_files")
                self.cwl_script.append(
                    f"      input_params: {{ default: {args_array} }}")
                self.cwl_script.append("      step_name:")
                self.cwl_script.append(f"        valueFrom: {task.task_id}")
                self.cwl_script.append(f"      output_filenames: {{ default: {output_files} }}")
                self.cwl_script.append(
                    "    out: [out, err, output_files]\n")

                # Adding a step to create a directory with the output files
                self.cwl_script.append(f"  {task.task_id}_folder:")
                self.cwl_script.append("    run: clt/folder.cwl")
                self.cwl_script.append("    in:")
                self.cwl_script.append("      - id: item")
                self.cwl_script.append("        linkMerge: merge_flattened")
                self.cwl_script.append("        source:")
                self.cwl_script.append(f"          - {task.task_id}/out")
                self.cwl_script.append(f"          - {task.task_id}/err")
                self.cwl_script.append(f"          - {task.task_id}/output_files")
                self.cwl_script.append("      - id: name")
                self.cwl_script.append(f"        valueFrom: \"{level}_{task.task_id}\"")
                self.cwl_script.append("    out: [out]\n")

                # adding the folder id to grand list of step folders
                steps_folder_source.append(f"{task.task_id}_folder")

        self.cwl_script.append("  final_folder:")
        self.cwl_script.append("    run: clt/folder.cwl")
        self.cwl_script.append("    in:")
        self.cwl_script.append("      - id: item")
        self.cwl_script.append("        linkMerge: merge_flattened")
        self.cwl_script.append("        source:")
        for folder in steps_folder_source:
            self.cwl_script.append(f"          - {folder}/out")
        self.cwl_script.append("      - id: name")
        self.cwl_script.append("        valueFrom: \"final_output\"")
        self.cwl_script.append("    out: [out]")

    def _parse_inputs_outputs(self) -> None:
        # Parsing the inputs of all root tasks
        self.cwl_script.append("\ninputs:")
        for task_name in self.root_task_names:
            task = self.tasks[task_name]
            cwl_written = False
            yml_written = False
            for f in task.input_files:
                if not cwl_written:
                    self.cwl_script.append(f"  {task.task_id}_input:")
                    self.cwl_script.append("    type: File[]")
                    cwl_written = True
                if not yml_written:
                    self.yml_script.append(f"{task.task_id}_input:")
                    yml_written = True

                self.yml_script.append("  - class: File")
                self.yml_script.append(f"    path: data/{f.file_id}")

        # Appending the output to the cwl script
        self.cwl_script.append("\noutputs:")
        self.cwl_script.append("  final_output_folder:")
        self.cwl_script.append("    type: Directory")
        self.cwl_script.append("    outputSource: final_folder/out")

    def _write_cwl_files(self, output_folder: pathlib.Path) -> None:
        cwl_folder = output_folder

        clt_folder = cwl_folder.joinpath("clt")
        clt_folder.mkdir(exist_ok=True)
        shutil.copy(this_dir.joinpath("templates/cwl_templates/wfbench.cwl"), clt_folder)
        shutil.copy(this_dir.joinpath("templates/cwl_templates/folder.cwl"), clt_folder)

        with open(cwl_folder.joinpath("main.cwl"), "w", encoding="utf-8") as f:
            f.write("\n".join(self.cwl_script))

        with (open(cwl_folder.joinpath("config.yml"), "w", encoding="utf-8")) as f:
            f.write("\n".join(self.yml_script))
