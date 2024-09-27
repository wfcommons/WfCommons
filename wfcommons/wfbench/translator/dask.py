#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib

from logging import Logger
from typing import Optional, Union

from .abstract_translator import Translator
from ...common import FileLink, Workflow

this_dir = pathlib.Path(__file__).resolve().parent


class DaskTranslator(Translator):
    """
    A WfFormat parser for creating Dask workflow applications.

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
        self.tasks_futures = {}
        self.task_id = 0

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        noindent_python_codelines = self._dask_wftasks_codelines("randomizer", output_folder)
        
        for task_name in self.root_task_names:
            noindent_python_codelines.extend(self._parse_tasks(task_name))
        
        # generate results
        while self.task_id > 0:
            self.task_id -= 1
            noindent_python_codelines.append(f"TASKS['{self.parsed_tasks[self.task_id]}'] = fut_dv_{self.task_id}.result()")

        # generate code
        INDENT = "    "
        wf_codelines = "\n".join(["%s%s" % (INDENT, codeline) for codeline in noindent_python_codelines])
        with open(this_dir.joinpath("templates/dask_template.py")) as fp:
            run_workflow_code = fp.read()
        run_workflow_code = run_workflow_code.replace("# Generated code goes here", wf_codelines)

        # write benchmark files
        output_folder.mkdir(parents=True)
        with open(output_folder.joinpath("dask_workflow.py"), "w") as fp:
            fp.write(run_workflow_code)

        # additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)
        
    def _dask_wftasks_codelines(self, 
                                randomizer_varname: str, 
                                output_folder: pathlib.Path,
                                simulate_minimum_execution_time: float = 0.1,
                                simulate_maximum_execution_time: float = 1.1) -> list[str]:
        """
        Build the code definining all tasks in the workflow, i.e. WorkflowTask instances.
        
        :param randomizer_varname: The name of the randomizer.
        :type randomizer_varname: str
        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path

        :return: The non-indented Python lines of code used to instantiate the WorkflowTask instances.
        :rtype: list[str]
        """
        codelines = ["randomizer = random.Random(seed)",
                     "TASKS = {}"]
        for task in self.tasks.values():
            input_files = [str(output_folder.joinpath(f"data/{f.file_id}")) for f in task.files if f.link == FileLink.INPUT]
            output_files = [str(output_folder.joinpath(f"data/{f.file_id}")) for f in task.files if f.link == FileLink.OUTPUT]
            program = output_folder.joinpath(f'bin/{task.program}')
            args = []
            print(task.args)
            for a in task.args:
                if "--out" in a:
                    a = a.replace("{", "\"{").replace("}", "}\"").replace(".txt'", ".txt\\\\\"").replace("'", "\\\\\"" + str(output_folder.joinpath("data")) + "/").replace(": ", ":")
                elif "--" not in a: 
                    a = str(output_folder.joinpath("data", a))
                else: 
                    a = a.replace("'", "\"") 
                args.append(a)
            print(args)
            print("")
            code = [f"WorkflowTask(dag_id = '{task.task_id}',",
                    f"             name = '{task.task_id}',",
                    f"             command_arguments = {[str(program)] + args},",
                    f"             inputs = {input_files},",
                    f"             outputs = {output_files},",
                    "             simulate = simulate,",
                    f"             randomizer = {randomizer_varname},",
                    f"             simulate_minimum_execution_time = {simulate_minimum_execution_time},",
                    f"             simulate_maximum_execution_time = {simulate_maximum_execution_time},",
                    "             )"]
            codelines.append(f"TASKS['{task.task_id}'] = {code[0]}")
            codelines.extend([codeline for codeline in code[1:]])
        # exit(1)
        return codelines

    def _parse_tasks(self, task_name: str) -> list[str]:
        """
        Recursively iterates over workflow tasks to generate submit command.
        
        :param task_name: The name of a task.
        :type task_name: str

        :return: The 
        :rtype: list[str]
        """
        if task_name not in self.parsed_tasks:
            # check for dependencies
            for parent in self.task_parents[task_name]:
                if parent not in self.parsed_tasks:
                    return []
            
            self.parsed_tasks.append(task_name)
            self.tasks_futures[task_name] = f"fut_dv_{self.task_id}"
            self.task_id += 1
            noindent_python_codelines = [f"{self.tasks_futures[task_name]} = client.submit(execute_task, TASKS['{task_name}'], {self.task_parents[task_name]})"]
            
            # parse children
            for child in self.task_children[task_name]:
                noindent_python_codelines.extend(self._parse_tasks(child))
        
        return noindent_python_codelines
