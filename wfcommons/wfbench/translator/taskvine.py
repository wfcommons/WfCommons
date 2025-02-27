#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import shutil

from logging import Logger
from typing import Optional, Union

from .abstract_translator import Translator
from ...common import Workflow

this_dir = pathlib.Path(__file__).resolve().parent


class TaskVineTranslator(Translator):
    """
    A WfFormat parser for creating TaskVine workflow applications.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path],
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """
    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow, logger)
        self.parsed_tasks = []
        self.task_counter = 1
        self.output_files_map = {}

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        self.script = "# workflow tasks\n"

        # add tasks per level
        self.next_level = self.root_task_names.copy()
        while self.next_level:
            self.next_level = self._add_level_tasks(self.next_level)
            self.script += "wait_for_tasks_completion()\n\n"

        # generate code
        run_workflow_code = self._merge_codelines("templates/taskvine_template.py", self.script)
    
        # write benchmark files
        output_folder.mkdir(parents=True)
        with open(output_folder.joinpath("taskvine_workflow.py"), "w") as fp:
            fp.write(run_workflow_code)

        # additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)
        shutil.copy(this_dir.joinpath("templates/taskvine_poncho.json"), output_folder)
        
    def _add_level_tasks(self, tasks_list: list[str]) -> list[str]:
        """
        Add all tasks from a level in the workflow.

        :param tasks_list: list of tasks in the level
        :type tasks_list: list[str]

        :return: List of next level tasks
        :rtype: list[str]
        """
        next_level = set()
        level_parsed_tasks = set()
        for task_name in tasks_list:
            if set(self.task_parents[task_name]).issubset(self.parsed_tasks):
                next_level.update(self._add_task(task_name))
                level_parsed_tasks.add(task_name)
            else:
                next_level.add(task_name)
        
        self.parsed_tasks.extend(list(level_parsed_tasks))
        return list(next_level)

    def _add_task(self, task_name: str, parent_task: Optional[str] = None) -> list[str]:
        """
        Add a task and its dependencies to the workflow.

        :param task_name: name of the task
        :type task_name: str
        :param parent_task: name of the parent task
        :type parent_task: Optional[str]

        :return: List of children tasks
        :rtype: list[str]
        """
        if task_name not in self.parsed_tasks:
            task = self.tasks[task_name]

            # input files
            f_counter = 1
            task_script = f"t_{self.task_counter}.add_poncho_package(poncho_pkg)\n" \
                           f"t_{self.task_counter}.add_input(wfbench, 'wfbench')\n" \
                           f"t_{self.task_counter}.add_input(cpu_bench, 'cpu-benchmark')\n" \
                           f"t_{self.task_counter}.add_input(stress_ng, 'stress-ng')\n"
            input_spec = "["
            for file in task.input_files:
                if file.file_id in self.output_files_map.keys():
                    task_script += f"t_{self.task_counter}.add_input({self.output_files_map[file.file_id]}, '{file}')\n"
                else:
                    task_script += f"in_{self.task_counter}_f_{f_counter} = m.declare_file('data/{file}')\n" \
                                    f"t_{self.task_counter}.add_input(in_{self.task_counter}_f_{f_counter}, '{file}')\n"
                f_counter += 1
                input_spec += f"\"{file.file_id}\","
            input_spec = input_spec[:-1] + "]"

            # output files
            f_counter = 1
            output_spec = "\"{"
            for file in task.output_files:
                output_spec += f"\\\\\"{file.file_id}\\\\\":{str(file.size)},"
                task_script += f"out_{self.task_counter}_f_{f_counter} = m.declare_file('outputs/{file}')\n" \
                                f"t_{self.task_counter}.add_output(out_{self.task_counter}_f_{f_counter}, '{file}')\n"
                self.output_files_map[file.file_id] = f"out_{self.task_counter}_f_{f_counter}"
                f_counter += 1
            output_spec = output_spec[:-1] + "}\""

            task_script += f"m.submit(t_{self.task_counter})\n" \
                            f"print(f'submitted task {{t_{self.task_counter}.id}}: {{t_{self.task_counter}.command}}')\n\n"            
            self.task_counter += 1

            # arguments
            args = []
            for a in task.args:
                if "--output-files" in a:
                    args.append(f"--output-files {output_spec}")
                elif "--input-files" in a:
                    args.append(f"--input-files {input_spec}")
                else:
                    args.append(a)
            args = " ".join(f"{a}" for a in args)

            # write task
            self.script += f"t_{self.task_counter} = vine.Task('{task.program} {args}')\n" \
                f"t_{self.task_counter}.set_cores(1)\n{task_script}"
            
            return self.task_children[task_name]
        
        return []
