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

from logging import Logger
from typing import Dict, Optional, Union

from .abstract_translator import Translator
from ...common import Workflow

this_dir = pathlib.Path(__file__).resolve().parent


class TaskVineTranslator(Translator):
    """
    A WfFormat parser for creating Pegasus workflow applications.

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

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        self.script = ""
        # adding tasks
        for task_name in self.root_task_names:
            # input file
            task = self.tasks[task_name]
            for file in task.input_files:
                self.script += f"if_{self.task_counter} = m.declare_file('data/{file.file_id}')\n"
            self._add_task(task_name)

        # generate code
        with open(this_dir.joinpath("templates/task_vine_template.py")) as fp:
            run_workflow_code = fp.read()
        run_workflow_code = run_workflow_code.replace("# Generated code goes here", self.script)
    
        # write benchmark files
        output_folder.mkdir(parents=True)
        with open(output_folder.joinpath("task_vine_workflow.py"), "w") as fp:
            fp.write(run_workflow_code)

        # additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)
        
    def _add_task(self, task_name: str, parent_task: Optional[str] = None) -> None:
        """
        Add a task and its dependencies to the workflow.

        :param task_name: name of the task
        :type task_name: str
        :param parent_task: name of the parent task
        :type parent_task: Optional[str]
        """
        if task_name not in self.parsed_tasks:
            task = self.tasks[task_name]
            # arguments
            args = []
            for a in task.args:
                a = a.replace("'", "\"") if "--out" not in a else a.replace("{", "\"{").replace("}", "}\"").replace("'", "\\\\\"").replace(": ", ":")
                args.append(a)
            args = " ".join(f"{a}" for a in args)

            self.script += f"t_{self.task_counter} = vine.Task('{args}')\n"

            self.task_counter += 1
