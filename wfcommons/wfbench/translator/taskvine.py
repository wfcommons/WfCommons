#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import shutil
import textwrap
from logging import Logger
from typing import Optional, Union

from .abstract_translator import Translator
from ...common import Workflow

this_dir = pathlib.Path(__file__).resolve().parent


def get_flowcept_init(workflow_id, workflow_name):
    code = textwrap.dedent(f"""
    from flowcept.flowcept_api.flowcept_controller import Flowcept
    f = Flowcept(workflow_id="{workflow_id}", workflow_name="{workflow_name}", bundle_exec_id="{workflow_id}")
    f.start()
    """)
    return code


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
                 with_flowcept: Optional[bool] = False,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow, logger)
        self.parsed_tasks = []
        self.task_counter = 1
        self.with_flowcept = with_flowcept
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
        with open(this_dir.joinpath("templates/taskvine_template.py")) as fp:
            run_workflow_code = fp.read()
        run_workflow_code = run_workflow_code.replace("# Generated code goes here", self.script)

        if self.with_flowcept:
            run_workflow_code = run_workflow_code.replace("# FLOWCEPT_INIT", get_flowcept_init(self.workflow.workflow_id, self.workflow.name))
            run_workflow_code = run_workflow_code.replace("# FLOWCEPT_END", "f.stop()")

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
            # arguments
            args = []
            for a in task.args:
                a = a.replace("'", "\"") if "--out" not in a else a.replace("{", "\"{").replace("}", "}\"").replace("'", "\\\\\"").replace(": ", ":")
                args.append(a)
            args = " ".join(f"{a}" for a in args)

            self.script += f"t_{self.task_counter} = vine.Task('{task.program} {args}')\n" \
                            f"t_{self.task_counter}.set_cores(1)\n"

            # input files
            f_counter = 1
            self.script += f"t_{self.task_counter}.add_poncho_package(poncho_pkg)\n" \
                            f"t_{self.task_counter}.add_input(wfbench, 'wfbench')\n" \
                            f"t_{self.task_counter}.add_input(cpu_bench, 'cpu-benchmark')\n" \
                            f"t_{self.task_counter}.add_input(stress_ng, 'stress-ng')\n"
            for in_file in task.input_files:
                if in_file.file_id in self.output_files_map.keys():
                    self.script += f"t_{self.task_counter}.add_input({self.output_files_map[in_file.file_id]}, '{in_file}')\n"
                else:
                    self.script += f"in_{self.task_counter}_f_{f_counter} = m.declare_file('data/{in_file}')\n" \
                                    f"t_{self.task_counter}.add_input(in_{self.task_counter}_f_{f_counter}, '{in_file}')\n"

            # output files
            f_counter = 1
            for out_file in task.output_files:
                self.script += f"out_{self.task_counter}_f_{f_counter} = m.declare_file('outputs/{out_file}')\n" \
                                f"t_{self.task_counter}.add_output(out_{self.task_counter}_f_{f_counter}, '{out_file}')\n"
                self.output_files_map[out_file.file_id] = f"out_{self.task_counter}_f_{f_counter}"
                f_counter += 1

            self.script += f"m.submit(t_{self.task_counter})\n" \
                            f"print(f'submitted task {{t_{self.task_counter}.id}}: {{t_{self.task_counter}.command}}')\n\n"

            self.task_counter += 1
            # self.parsed_tasks.append(task_name)
            
            return self.task_children[task_name]
        
        return []
