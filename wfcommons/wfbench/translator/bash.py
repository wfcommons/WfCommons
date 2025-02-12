#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging
from typing import Union, Optional
from collections import defaultdict, deque
import pathlib
import ast
import json
from ...common import Workflow
from .abstract_translator import Translator

this_dir = pathlib.Path(__file__).resolve().parent


class BashTranslator(Translator):
    """
    A WfFormat parser for creating a sequential bash workflow benchmarks.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path],
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """
    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 logger: Optional[logging.Logger] = None) -> None:
        super().__init__(workflow, logger)
        self.parsl_script = []
        self.task_level_map = defaultdict(lambda: [])

        indegree = {}

        for task in self.tasks.values():
            indegree[task.task_id] = len(self.task_parents[task.task_id])
        queue = deque(self.root_task_names)
        top_sort = []

        while queue:
            task_name = queue.popleft()
            top_sort.append(task_name)

            for child in self.task_children[task_name]:
                indegree[child] -= 1
                if indegree[child] == 0:
                    queue.append(child)

        assert (len(top_sort) == len(self.tasks)), "Error: The workflow contains a cycle"

        levels = {task_name: 0 for task_name in top_sort}

        for task_name in top_sort:
            for child in self.task_children[task_name]:
                levels[child] = max(levels[child], levels[task_name] + 1)

        for task_name, level in levels.items():
            self.task_level_map[level].append(task_name)

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """

        # Parsing each of the WfFormat Tasks as bash apps in Parsl
        codelines = self._bash_wftasks_codelines()

        wf_codelines = "\n".join(codelines)

        # Generate an output folder
        output_folder.mkdir(parents=True)
        with open(output_folder.joinpath("run_workflow.sh"), "w", encoding="utf-8") as fp:
            fp.write(wf_codelines)

        # Additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)

    def _bash_wftasks_codelines(self) -> None:
        codelines = []

        # Parsing each steps by Workflow levels
        for level in sorted(self.task_level_map.keys()):
            # Parsing each task within a Workflow level
            for task_name in self.task_level_map[level]:
                # Getting the task object
                task = self.tasks[task_name]

                args = []
                for a in task.args:
                    if a.startswith("--output-files"):
                        flag, output_files_dict = a.split(" ", 1)
                        output_files_dict = {f"data/{key}": value for key, value in ast.literal_eval(output_files_dict).items()}
                        a = f"{flag} '{json.dumps(output_files_dict).replace('"', '\\"')}'"

                    if a.startswith("--input-files"):
                        flag, input_files_arr = a.split(" ", 1)
                        input_files_arr = [f"data/{file}" for file in ast.literal_eval(input_files_arr)]
                        a = f"{flag} '{json.dumps(input_files_arr).replace('"', '\\"')}'"

                    args.append(a)

                code = f"bin/{task.program} {' '.join(args)}"

                codelines.append(code)

        return codelines
