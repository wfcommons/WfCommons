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
import ast
import json
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
                           "  InlineJavascriptRequirement: {}\n"]
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

    def _parse_steps(self) -> None:
        self.cwl_script.append("steps:")

        output_files_sources = []
        log_files_sources = []

        # Parsing each steps by Workflow levels
        for level in sorted(self.task_level_map.keys()):
            # Parsing each task within a Workflow level
            for task_name in self.task_level_map[level]:

                # Getting the task object
                task = self.tasks[task_name]

                # Parsing the arguments of the task
                args_array = []

                for a in task.args:
                    if a.startswith("--output-files"):
                        flag, output_files_dict = a.split(" ", 1)
                        output_files_dict = {f"{key}": value for key, value in ast.literal_eval(output_files_dict).items()}
                        a = f"{flag} '{json.dumps(output_files_dict).replace('"', '\\"')}'"
                    if a.startswith("--input-files"):
                        flag, input_files_arr = a.split(" ", 1)
                        input_files_arr = [f"{file}" for file in ast.literal_eval(input_files_arr)]
                        a = f"{flag} '{json.dumps(input_files_arr).replace('"', '\\"')}'"
                    args_array.append(a)


                cmd = f"{task.program} {' '.join(args_array)}"

                output_files = [
                    f.file_id for f in task.output_files]

                if level == 0:
                    input_files = [f"      input_files: {task.task_id}_input"]
                else:
                    input_files = [  "      input_files:",
                                     "        linkMerge: merge_flattened",
                                     "        source:"]
                    input_files += [f"          - {p}/output_files" for p in self.task_parents[task.task_id]]

                code = [
                    f"  {task.task_id}:",
                    "    run: clt/bash.cwl",
                    "    in:",
                ]

                code += input_files

                code += [
                    f"      command: {{default: \"{cmd}\"}}",
                     "      step_name:",
                    f"        valueFrom: \"{task.task_id}\"",
                    f"      output_filenames: {{default: {output_files}}}",
                     "    out: [out, err, output_files]\n"
                ]

                self.cwl_script.extend(code)
                output_files_sources.append(f"          - {task.task_id}/output_files")
                log_files_sources.append(f"          - {task.task_id}/out")
                log_files_sources.append(f"          - {task.task_id}/err")

        code = [
             "  compile_output_files:",
             "    run: clt/folder.cwl",
             "    in:",
             "      - id: name",
             "        valueFrom: \"output\"",
             "      - id: item",
             "        linkMerge: merge_flattened",
             "        source:",
        ]

        code += output_files_sources

        code += [
             "    out: [out]\n"
        ]

        code += [
             "  compile_log_files:",
             "    run: clt/folder.cwl",
             "    in:",
             "      - id: name",
             "        valueFrom: \"logs\"",
             "      - id: item",
             "        linkMerge: merge_flattened",
             "        source:",
        ]

        code += log_files_sources

        code += [
             "    out: [out]\n"
        ]

        self.cwl_script.extend(code)

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
        code = ["\noutputs:",
                "  data_folder:",
                "    type: Directory",
                "    outputSource: compile_output_files/out",
                "  log_folder:",
                "    type: Directory",
                "    outputSource: compile_log_files/out\n"]

        self.cwl_script.extend(code)

    def _write_cwl_files(self, output_folder: pathlib.Path) -> None:
        cwl_folder = output_folder

        clt_folder = cwl_folder.joinpath("clt")
        clt_folder.mkdir(exist_ok=True)
        shutil.copy(this_dir.joinpath("templates/cwl_templates/wfbench.cwl"), clt_folder)
        shutil.copy(this_dir.joinpath("templates/cwl_templates/folder.cwl"), clt_folder)
        shutil.copy(this_dir.joinpath("templates/cwl_templates/bash.cwl"), clt_folder)

        with open(cwl_folder.joinpath("main.cwl"), "w", encoding="utf-8") as f:
            f.write("\n".join(self.cwl_script))

        with (open(cwl_folder.joinpath("config.yml"), "w", encoding="utf-8")) as f:
            f.write("\n".join(self.yml_script))
