#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import shutil

from logging import Logger
from typing import List, Optional, Union

from .abstract_translator import Translator
from ...common import Workflow
from ...common.task import Task

this_dir = pathlib.Path(__file__).resolve().parent


class NextflowTranslator(Translator):
    """
    A WfFormat parser for creating Nextflow workflow applications.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path]
    :param slurm: Whether to generate a Slurm template script for workflow submission using :code:`sbatch`.
    :type slurm: bool
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """
    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 slurm: Optional[bool] = False,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow, logger)
        self.slurm = slurm
        self.script = ""
        self.out_files = set()

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description(WfFormat) into a Nextflow workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        # Create the output folder
        self.output_folder = output_folder
        self.output_folder.mkdir(parents=True)

        # Create benchmark files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)
        
        if self.slurm:
            shutil.copy(this_dir.joinpath("templates/nextflow/nextflow_hyperqueue_job.sh"), output_folder)

        if self.workflow.workflow_id:
            shutil.copy(this_dir.joinpath("templates/flowcept_agent.py"), output_folder.joinpath("bin"))
            self.logger.info(f"Workflow ID: {self.workflow.workflow_id}")

        # Create a topological order of the tasks
        sorted_tasks = self._get_tasks_in_topological_order()

        # Create the bash script for each task
        for task in sorted_tasks:
            self._create_task_script(task)

        # Create the Nextflow workflow script and file
        self._create_workflow_script(sorted_tasks)
        run_workflow_code = self._merge_codelines("templates/nextflow/workflow.nf", self.script)
        self._write_output_file(run_workflow_code, output_folder.joinpath("workflow.nf"))

        # Create the README file
        self._write_readme_file(output_folder)

        return

    def _create_workflow_script(self, tasks: list[Task]):
        """
        Create the Nextflow script.

        :param tasks: The (sorted) list of tasks.
        :type tasks: list[Task]
        """       
        # Add Flowcept code if enabled
        if self.workflow.workflow_id:
            self.script += self._generate_flowcept_code()

        # Output the code for each task
        for task in tasks:
            self.script += self._generate_task_code(task)

        # Output the code for the workflow
        self.script += self._generate_workflow_code(tasks)

        return

    def _generate_flowcept_code(self) -> str:
        """

        :return: The code.
        :rtype: str
        """
        out_files = ", ".join(f"\"{item}\"" for item in self.out_files)
        return "process flowcept(){\n" \
               "    input:\n" \
               "    output:\n" \
               "    script:\n" \
               "        \"\"\"\n" \
		       "        ${pwd}/bin/flowcept_agent.py " \
               f"{self.workflow.name} {self.workflow.workflow_id} '[{out_files}]' \n" \
		       "        \"\"\"\n" \
               "}\n\n"                     

    def _get_tasks_in_topological_order(self) -> List[Task]:
        """
        Sort the workflow tasks in topological order.

        :param output_folder: The path to the output folder.
        :type output_folder: pathlib.Path
        :return: A sorted list of tasks.
        :rtype: List[Task]
        """
        levels = {0: self._find_root_tasks()}
        sorted_tasks: List[Task] = levels[0]
        current_level = 1
        while True:
            tasks_in_current_level = []
            all_children = [self._find_children(p.task_id) for p in levels[current_level-1]]
            all_children = [item for sublist in all_children for item in sublist]
            all_children = list(set(all_children))
            if not all_children:
                break
            for potential_task in all_children:
                num_children = len(self.task_children[potential_task.task_id])
                if not num_children:
                    self.out_files.add(f"{self.output_folder.absolute()}/data/{potential_task.output_files[0]}")
                if all(parent in sorted_tasks for parent in self._find_parents(potential_task.task_id)):
                    tasks_in_current_level.append(potential_task)
            levels[current_level] = tasks_in_current_level
            sorted_tasks += tasks_in_current_level
            current_level += 1
        return sorted_tasks

    def _create_task_script(self, task: Task):
        """
        Generate the bash script for invoking a task.

        :param task: The task.
        :type task: Task
        :return: The code.
        :rtype: str
        """

        code = "#!/bin/bash\n\n"

        # Generate input spec
        input_spec = "'\\["
        for f in task.input_files:
            input_spec += f"\"{self.output_folder.resolve()}/data/{f.file_id}\","
        input_spec = input_spec[:-1] + "\\]'"

        # Generate output spec
        output_spec = "'\\{"
        for f in task.output_files:
            output_spec += f"\"{self.output_folder.resolve()}/data/{f.file_id}\":{str(f.size)},"
        output_spec = output_spec[:-1] + "\\}'"

        code += f"{self.output_folder.resolve()}/bin/{task.program} "

        for a in task.args:
            if "--output-files" in a:
                code += f"--output-files {output_spec} "
            elif "--input-files" in a:
                code += f"--input-files {input_spec} "
            else:
                code += f"{a} "
        code += "\n"

        script_file_path = self.output_folder.joinpath(f"bin/script_{task.task_id}.sh")
        with open(script_file_path, "w") as out:
            out.write(code)

    def _generate_task_code(self, task: Task) -> str:
        """
        Generate the code for a task, as a Nextflow process.

        :param task: The task.
        :type task: Task
        :return: The code.
        :rtype: str
        """
        function_name = task.task_id.replace(".", "_")
        code = f"process {function_name}()" + "{\n"

        # File variables
        if self._find_children(task.task_id):
            for f in task.output_files:
                code += f"\tdef {f.file_id} = " + "\"${pwd}/data/" + f.file_id + "\"\n"

        # Input declaration
        code += f"\tinput:\n"
        if self._find_parents(task.task_id):
            code += f"\t\tval {task.task_id}_ifs\n"

        # Output declaration
        code += f"\toutput:\n"
        if self._find_children(task.task_id):
            code += "\t\tval([\n"
            for f in task.output_files:
                code += f"\t\t\t{f.file_id},\n"
            code = code[:-2]
            code += "\n\t\t])\n"

        # Script
        code += "\tscript:\n"
        # Input convenience vars
        if self._find_parents(task.task_id):
            counter = 0
            for f in task.input_files:
                code += f"\t\tdef {f.file_id} = {task.task_id}_ifs[{counter}]\n"
                counter += 1
            code += "\n"

        # Generate output variables
        if self._find_children(task.task_id):
            for f in task.output_files:
                code += "\t\t" + f.file_id + " = \"${pwd}/data/" + f.file_id + "\"\n"


        code += "\t\t\"\"\"\n"
        code += "\t\t${params.simulate ? 'sleep 1' : \"bash ${pwd}/bin/script_" + task.task_id + ".sh\"}\n"
        code += "\t\t\"\"\"\n"
        code += "}\n\n"
        return code

    def _generate_workflow_code(self, sorted_tasks: List[Task]) -> str:
        """
        Generate the code for the workflow.

        :param sorted_tasks: The task.
        :type sorted_tasks: List[Task]
        :return: The code.
        :rtype: str
        """

        code = ""

        # Generate bootstrap function
        code += ("// Bootstrap function\n"
                 "def bootstrap() {\n"
                 "\tdef outputs = [:]\n"
                 "\treturn outputs\n"
                 "}\n\n")

        # Generate all task functions
        for task in sorted_tasks:
            code += self._generate_task_function(task)

        # Generate workflow function
        code += "workflow {\n"
        if self.workflow.workflow_id:
            code += "\tflowcept()\n"
        code += "\tresults = bootstrap()\n"
        for task in sorted_tasks:
            function_name = task.task_id.replace(".", "_")
            code += f"\tresults = function_{function_name}(results)\n"
        code += "}\n"
        return code

    def _generate_task_function(self, task: Task) -> str:
        """
        Generate the code for a starting a task's Nextflow process, as a Nextflow function.

        :param task: The task.
        :type task: Task
        :return: The code.
        :rtype: str
        """
        code = f"// Function to call task {task.task_id}\n"
        function_name = task.task_id.replace(".", "_")
        code += f"def function_{function_name}(Map inputs) " + "{\n"

        if self._find_parents(task.task_id):
            # Input channel mixing and then call
            code += f"\tdef {function_name}_necessary_input = Channel.empty()\n"
            for f in task.input_files:
                code += f"\t{function_name}_necessary_input = {function_name}_necessary_input.mix(inputs.{f.file_id})\n"
            code += f"\tdef {function_name}_necessary_input_future = {function_name}_necessary_input.collect()\n"
            code += f"\tdef {function_name}_produced_output = {function_name}({function_name}_necessary_input_future)\n"
        else:
            # Simple call
            code += f"\tdef {function_name}_produced_output = {function_name}()\n"

        # Pass on the outputs
        code += "\n"
        code += "\tdef outputs = inputs.clone()\n"
        if self._find_children(task.task_id):
            counter = 0
            for f in task.output_files:
                code += f"\toutputs.{f.file_id} = {function_name}_produced_output.map" + "{it[" + str(counter) + "]}\n"
                counter += 1
        code += "\treturn outputs\n"
        code += "}\n\n"

        return code


    def _write_readme_file(self, output_folder: pathlib.Path) -> None:
        """
        Write the README  file.

        :param output_folder: The path of the output folder.
        :type output_folder: pathlib.Path
        """
        readme_file_path = output_folder.joinpath("README")
        with open(readme_file_path, "w") as out:
            out.write(f"Run the workflow in directory {str(output_folder)} using the following command:\n")

            out.write(f"\tnextflow run ./workflow.nf --pwd `pwd`\n")
