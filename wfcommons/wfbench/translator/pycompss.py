#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib

from logging import Logger
from typing import Optional, Union

from .abstract_translator import Translator
from ...common import Workflow

this_dir = pathlib.Path(__file__).resolve().parent


class PyCompssTranslator(Translator):
    """
    A WfFormat parser for creating PyCOMPSs workflow applications.

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
        self.output_folder = pathlib.Path("")

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        self.output_folder = output_folder
        self.script = "\n# workflow tasks\n"
        # PyCOMPSs translator
        self._pycompss_code()

        # Generates pycompss workflow file: template + script
        with open(this_dir.joinpath("templates/pycompss_template.py")) as fp:
            run_workflow_code = fp.read()
        run_workflow_code = run_workflow_code.replace("# Generated code goes here", self.script)
        # write benchmark files
        output_folder.mkdir(parents=True)
        with open(output_folder.joinpath("pycompss_workflow.py"), "w") as fp:
            fp.write(run_workflow_code)
        # additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)


    def _pycompss_code(self) -> None:
        # GENERATES PYCOMPSS TASKS (functions)
        bin_path = self.output_folder.joinpath("bin")
        all_pycompss_tasks_as_functions = []
        for task in self.tasks.values():
            ############################
            # CREATE FUNCTION DECORATOR: @binary parameters
            ############################
            self.script += f"@binary(binary='${bin_path}/{task.program}'"
            if len(task.args) > 0:
                task_args = " ".join(task.args)
                self.script += f", args='{task_args}'"
            self.script += f")"
            # @binary(binary="date", args= "-d {{param_1}}")
            ############################
            # CREATE FUNCTION DECORATOR: @task parameters
            ############################
            task_parameter_names_file_in = ""
            for i in range(len(task.input_files)):
                if len(task.input_files) == 1:
                    task_parameter_names_file_in += f"file_in_{i}=FILE_IN"
                else:
                    if i == 0:
                        task_parameter_names_file_in += f"file_in_{i}=FILE_IN"
                    else:
                        task_parameter_names_file_in += f", file_in_{i}=FILE_IN"
            task_parameter_names_file_out = ""
            for i in range(len(task.output_files)):
                if len(task.output_files) == 1:
                    task_parameter_names_file_out += f"file_out_{i}=FILE_OUT"
                else:
                    if i == 0:
                        task_parameter_names_file_out += f"file_out_{i}=FILE_OUT"
                    else:
                        task_parameter_names_file_out += f", file_out_{i}=FILE_OUT"
            self.script += f"@task"
            self.script += "("
            self.script += task_parameter_names_file_in
            self.script += ", " if len(task_parameter_names_file_in) > 0 else ""
            self.script += task_parameter_names_file_out
            self.script += ")\n"
            ############################
            # CREATE FUNCTION DEFINITION
            ############################
            function_name = task.name
            # function parameters
            function_parameter_names_file_in = ""
            for i in range(len(task.input_files)):
                if len(task.input_files) == 1:
                    function_parameter_names_file_in += f"file_in_{i}"
                else:
                    if i == 0:
                        function_parameter_names_file_in += f"file_in_{i}"
                    else:
                        function_parameter_names_file_in += f", file_in_{i}"
            function_parameter_names_file_out = ""
            for i in range(len(task.output_files)):
                if len(task.output_files) == 1:
                    function_parameter_names_file_out += f"file_out_{i}"
                else:
                    if i == 0:
                        function_parameter_names_file_out += f"file_out_{i}"
                    else:
                        function_parameter_names_file_out += f", file_out_{i}"
            self.script += f"def {function_name}"
            self.script += "("
            self.script += function_parameter_names_file_in
            self.script += ", " if len(function_parameter_names_file_in) > 0 else ""
            self.script += function_parameter_names_file_out
            self.script += "):\n"
            ############################
            # CREATE FUNCTION BODY
            ############################
            for outfile in function_parameter_names_file_out.replace(' ', '').split(','):
                # this method is in the template file 'pycompss_template.py'
                self.script += f"\t_create_out_file({outfile})\n"
            self.script += f"\tpass\n\n"

            ############################
            # STORE FUNCTION CALL
            ############################
            is_root_task = True if task.name in self.root_task_names else False
            data_folder = self.output_folder.joinpath("data")
            function_parameters_in = ""
            for i in range(len(task.input_files)):
                if len(task.input_files) == 1:
                    if is_root_task:
                        function_parameters_in += f"\'{data_folder}/{task.input_files[i].file_id}\'"
                    else:
                        function_parameters_in += f"\'{task.input_files[i].file_id}\'"
                else:
                    if i == 0:
                        if is_root_task:
                            function_parameters_in += f"\'{data_folder}/{task.input_files[i].file_id}\'"
                        else:
                            function_parameters_in += f"\'{task.input_files[i].file_id}\'"
                    else:
                        if is_root_task:
                            function_parameters_in += f", \'{data_folder}/{task.input_files[i].file_id}\'"
                        else:
                            function_parameters_in += f", \'{task.input_files[i].file_id}\'"
            function_parameters_out = ""
            for i in range(len(task.output_files)):
                if len(task.output_files) == 1:
                    function_parameters_out += f"\'{task.output_files[i].file_id}\'"
                else:
                    if i == 0:
                        function_parameters_out += f"\'{task.output_files[i].file_id}\'"
                    else:
                        function_parameters_out += f", \'{task.output_files[i].file_id}\'"
            function_parameters_in_out = ""
            function_parameters_in_out += function_parameters_in
            function_parameters_in_out += ", " if len(function_parameters_in) > 0 else ""
            function_parameters_in_out += function_parameters_out
            all_pycompss_tasks_as_functions.append(f"{function_name}({function_parameters_in_out})")


        # INVOKE PYCOMPSS TASKS (functions)
        self.script += f"\n\ndef main_program():\n"
        # tasks
        for func in all_pycompss_tasks_as_functions:
            self.script += f"\t{func}\n"

        # CALL TO MAIN METHOD
        self.script += f"\n\nif __name__ == \"__main__\":\n"
        self.script += f"\tmain_program()\n"
