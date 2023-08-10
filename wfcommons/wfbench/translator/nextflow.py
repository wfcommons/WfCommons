#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2022 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib

from logging import Logger
from typing import Dict, List, Optional, Union

from .abstract_translator import Translator
from ...common import FileLink, Workflow
from ...common.task import Task

class NextflowTranslator(Translator):
    """
    A WfFormat parser for creating Nextflow workflow applications.

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

        self.script = "nextflow.enable.dsl=2\n\n"

    def translate(self, output_file_path: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into a Nextflow workflow application.

        :param output_file_path: The name of the output file (e.g., workflow.py).
        :type output_file_path: pathlib.Path
        """
        
        self.all_inputs = []
        self.all_outputs = []
        self.task_inputs: Dict[str, List[str]] = {}
        self.task_outputs: Dict[str, List[str]] = {}
        
        task_written: Dict[str, bool] = {}
        for task in self.tasks.values():
            task_written[task.name] = False
            self._add_task(task)

        self.workflow_inputs = [i for i in self.all_inputs if i not in self.all_outputs]
        self.workflow_outputs = [o for o in self.all_outputs if o not in self.all_outputs]
        
        self.script += "workflow {\n"
        self.script += "    data = channel.fromPath('./inputs/*')\n"

        def _write_task(task_name: str, inputs: List[str]):
            task_params = []
            for parent in self.task_parents[task_name]:
                if not task_written[parent]:
                    _write_task(parent, self.task_inputs[parent])
            i = 0
            for input_name in self.task_inputs[task_name]:
                for parent in self.task_parents[task_name]:
                    if input_name in self.task_outputs[parent]:
                        task_params.append(f"{parent}_out.get({self.task_outputs[parent].index(input_name)})")
                        break
                else:
                    if input_name in self.workflow_inputs:
                        continue
                    varname = f"{task_name}_in{i}"
                    i += 1
                    self.script += f"    file('{input_name.replace('-', '_')}').text = \"=== FILE {input_name.replace('-', '_')} ===\"\n"
                    self.script += f"    {varname.replace('-', '_')} = channel.fromPath(file('{input_name.replace('-', '_')}'))\n"
                    task_params.append(varname)
            if not set(self.workflow_inputs).isdisjoint(inputs):
                task_params.append("data")
            task_params = {param.replace('-', '_') for param in task_params}
            self.script += f"    {task_name.replace('-', '_')}_out = {task_name.replace('-', '_')}({', '.join(task_params)})\n"
            task_written[task_name] = True

        for task_name, inputs in self.task_inputs.items():
            if not task_written[task_name]:
                _write_task(task_name, inputs)

        self.script += "}\n"

        # write script to file
        self._write_output_file(self.script, output_file_path)

    def _add_task(self, task: Task) -> None:
        """
        Add a task to the workflow.

        :param task: the task
        :type task: Task
        """

        input_files = [file.name for file in task.files if file.link == FileLink.INPUT]
        output_files = [file.name for file in task.files if file.link == FileLink.OUTPUT]
        self.all_inputs.extend(input_files)
        self.all_outputs.extend(output_files)
    
        self.script += f"process {task.name.replace('-', '_')}" + " {\n"
        if task.cores:
            self.script += f"  cpus {task.cores+1}\n"
        if task.memory:
            self.script += f"  memory {task.memory}\n"
        else:
            self.script += f"  memory 10.GB\n"
        self.script += "  errorStrategy { task.exitStatus=143 ? 'ignore' : 'terminate' }\n\n"
        self.script += "  input:\n"
        for input_file in input_files:
            self.script += f"    path \"{input_file}\"\n"
        self.script += "  output:\n"
        for output_file in output_files:
            self.script += f"    path \"{output_file}\"\n"
        self.script += "\n"
        self.script += "  script:\n"
        self.script += "  \"\"\"\n"
        self.script += f"  {pathlib.Path(task.program).name}"
        for a in task.args:
            self.script += ' '
            a = a.replace("'", "\"") if "--out" not in a else a.replace("{", "\"{").replace("}", "}\"").replace("'", "\\\\\"").replace(": ", ":")
            self.script += a
        self.script += '\n'
        self.script += "  \"\"\"\n"
        self.script += "}\n"

        self.task_inputs[task.name] = input_files
        self.task_outputs[task.name] = output_files

