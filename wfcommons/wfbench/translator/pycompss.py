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
import ast
from logging import Logger
from typing import Optional, Union
import copy
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
        self.script = ""

        # PyCOMPSs translator
        self.script += "\n# workflow tasks\n"
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
        bin_path = "${WFBENCH_BIN}"
        data_folder = "os.getenv('WFBENCH_DATA')"
        all_pycompss_tasks_as_functions = {}
        task_number = 1
        for task in self.tasks.values():
            function_name = f"{task.name}{task_number}"
            task_number += 1
            is_root_task = True if task.task_id in self.root_task_names else False
            all_input_files_name = []
            task_parameter_names_file_in = ""
            function_parameter_names_file_in = ""
            function_parameters_in = ""
            task_parameter_names_file_out = ""
            function_parameter_names_file_out = ""
            function_parameters_out = ""

            for i in range(len(task.input_files)):
                all_input_files_name.append(task.input_files[i].file_id)
                if len(task.input_files) == 1:
                    task_parameter_names_file_in += f"file_in_{i}=FILE_IN"
                    function_parameter_names_file_in += f"file_in_{i}"
                    if is_root_task:
                        function_parameters_in += "f\"{" + data_folder + "}" + f"/{task.input_files[i].file_id}" + "\""
                    else:
                        function_parameters_in += f"\'{task.input_files[i].file_id}\'"
                else:
                    if i == 0:
                        task_parameter_names_file_in += f"file_in_{i}=FILE_IN"
                        function_parameter_names_file_in += f"file_in_{i}"
                        if is_root_task:
                            function_parameters_in += "f\"{" + data_folder + "}" + f"/{task.input_files[i].file_id}" + "\""
                        else:
                            function_parameters_in += f"\'{task.input_files[i].file_id}\'"
                    else:
                        task_parameter_names_file_in += f", file_in_{i}=FILE_IN"
                        function_parameter_names_file_in += f", file_in_{i}"
                        if is_root_task:
                            function_parameters_in += ", f\"{" + data_folder + "}" + f"/{task.input_files[i].file_id}" + "\""
                        else:
                            function_parameters_in += f", \'{task.input_files[i].file_id}\'"
            for i in range(len(task.output_files)):
                if len(task.output_files) == 1:
                    task_parameter_names_file_out += f"file_out_{i}=FILE_OUT"
                    function_parameter_names_file_out += f"file_out_{i}"
                    function_parameters_out += f"\'{task.output_files[i].file_id}\'"
                else:
                    if i == 0:
                        task_parameter_names_file_out += f"file_out_{i}=FILE_OUT"
                        function_parameter_names_file_out += f"file_out_{i}"
                        function_parameters_out += f"\'{task.output_files[i].file_id}\'"
                    else:
                        task_parameter_names_file_out += f", file_out_{i}=FILE_OUT"
                        function_parameter_names_file_out += f", file_out_{i}"
                        function_parameters_out += f", \'{task.output_files[i].file_id}\'"
            ############################
            # STORE FUNCTION CALL
            ############################
            function_parameters_in_out = ""
            function_parameters_in_out += function_parameters_in
            function_parameters_in_out += ", " if len(function_parameters_in) > 0 else ""
            function_parameters_in_out += function_parameters_out
            all_pycompss_tasks_as_functions[task.task_id] = f"{function_name}({function_parameters_in_out})"
            ############################
            # create function decorator: @binary parameters
            ############################
            self.script += f"@binary(binary='{bin_path}/{task.program}'"
            if len(task.args) > 0:
                all_task_args = ""
                for task_arg in task.args:
                    if task_arg.startswith("--output-files"):
                        all_out_params = function_parameter_names_file_out.replace(' ', '').split(',')
                        i = 0
                        all_task_args += "--output-files "
                        all_task_args += "{"
                        if task_arg.find('"{') != -1:
                            task_arg = task_arg.replace('"{', '{')
                            task_arg = task_arg.replace('}"', '}')
                        if task_arg.find('\\\"') != -1:
                            task_arg = task_arg.replace('\\\"', '"')
                        for file_out_name, file_out_size in ast.literal_eval(task_arg.split('--output-files ')[1]).items():
                            if i == 0:
                                all_task_args += "\\\\\\\"" + "{{"+ all_out_params[i] +"}}" + "\\\\\\\":" + str(file_out_size)
                            else:
                                all_task_args += ", \\\\\\\"" + "{{" + all_out_params[i] + "}}" + "\\\\\\\":" + str(file_out_size)
                            i += 1
                        all_task_args += "} "
                    elif task_arg.startswith("--input-files"):
                        all_task_args += "--input-files "
                        all_task_args += "["
                        i = 0
                        for input_param in function_parameter_names_file_in.replace(' ', '').split(','):
                            if i == 0:
                                all_task_args += "\\\\\\\"" + "{{"+ input_param +"}}" + "\\\\\\\""
                            else:
                                all_task_args += ", \\\\\\\"" + "{{" + input_param + "}}" + "\\\\\\\""
                            i += 1
                        all_task_args += "] "
                    else:
                        all_task_args += f"{task_arg} "
                self.script += f", args='{all_task_args}'"
            self.script += f")\n"


            ############################
            # CREATE FUNCTION DECORATOR: @task parameters
            ############################
            self.script += f"@task"
            self.script += "("
            self.script += task_parameter_names_file_in
            self.script += ", " if len(task_parameter_names_file_in) > 0 else ""
            self.script += task_parameter_names_file_out
            self.script += ")\n"
            ############################
            # CREATE FUNCTION DEFINITION: function parameters
            ############################
            self.script += f"def {function_name}"
            self.script += "("
            self.script += function_parameter_names_file_in
            self.script += ", " if len(function_parameter_names_file_in) > 0 else ""
            self.script += function_parameter_names_file_out
            self.script += "):\n"
            ############################
            # CREATE FUNCTION BODY
            ############################
            self.script += f"\tpass\n\n"



        # INVOKE PYCOMPSS TASKS (functions)
        self.script += f"\n\ndef main_program():\n"
        added_tasks = []
        copy_of_task_parents = copy.deepcopy(self.task_parents)
        # call root tasks first (no parents)
        for task_id, my_parents in copy_of_task_parents.items():
            if len(my_parents) == 0:
                self.script += f"\t{all_pycompss_tasks_as_functions.pop(task_id)}\n"
                added_tasks.append(task_id)
        # remove added tasks from dict
        for task_id in added_tasks:
            copy_of_task_parents.pop(task_id) if task_id in copy_of_task_parents else None

        # call tasks with parents
        while len(all_pycompss_tasks_as_functions) > 0:
            for task_id, my_parents in copy_of_task_parents.items():
                is_ready = False
                for parent in my_parents:
                    if parent not in added_tasks:
                        is_ready = False
                        break
                    else:
                        is_ready = True
                if is_ready:
                    self.script += f"\t{all_pycompss_tasks_as_functions.pop(task_id)}\n"
                    added_tasks.append(task_id)
            # remove added tasks from dict
            for task_id in added_tasks:
                copy_of_task_parents.pop(task_id) if task_id in copy_of_task_parents else None

        # CALL TO MAIN METHOD
        self.script += f"\n\nif __name__ == \"__main__\":\n"
        # START Flowcept
        if self.workflow.workflow_id is not None:
            flowcept_init_code = self._flowcept_init_python(self.workflow.workflow_id, self.workflow.name)
            self.script += "".join("\t" + line + "\n" for line in flowcept_init_code.splitlines())
        # main
        self.script += f"\tmain_program()\n"
        # STOP Flowcept
        if self.workflow.workflow_id is not None:
            self.script += f"\t{self._flowcept_stop_python()}\n"