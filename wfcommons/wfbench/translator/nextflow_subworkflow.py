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
from typing import Dict, List, Optional, Set, Union

from .abstract_translator import Translator
from ...common import Workflow
from ...common.task import Task

this_dir = pathlib.Path(__file__).resolve().parent


class NextflowSubworkflowTranslator(Translator):
    """
    A WfFormat parser for creating Nextflow workflow applications with code split across multiple files.

    This translator splits large workflows into multiple .nf files (max tasks per file),
    enabling better scalability for workflows with thousands of tasks. The processes and
    helper functions are split across files, while orchestration remains in the main workflow.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path]
    :param max_tasks_per_subworkflow: Maximum number of tasks per subworkflow file (default: 1000).
    :type max_tasks_per_subworkflow: int
    :param slurm: Whether to generate a Slurm template script for workflow submission using :code:`sbatch`.
    :type slurm: bool
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """
    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 max_tasks_per_subworkflow: Optional[int] = 1000,
                 slurm: Optional[bool] = False,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow, logger)
        self.max_tasks_per_subworkflow = max_tasks_per_subworkflow
        self.slurm = slurm
        self.script = ""
        self.out_files = set()
        self.subworkflows: List[List[Task]] = []

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description(WfFormat) into a Nextflow workflow application
        with code split across multiple files.

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

        # Split tasks into chunks for separate files
        self._split_into_subworkflows(sorted_tasks)

        self.logger.info(f"Split workflow into {len(self.subworkflows)} file(s)")

        # Create the bash script for each task
        for task in sorted_tasks:
            self._create_task_script(task)

        # Create modules directory
        modules_dir = output_folder.joinpath("modules")
        modules_dir.mkdir(exist_ok=True)

        # Create each module file (processes + helper functions)
        for idx, subworkflow_tasks in enumerate(self.subworkflows):
            module_code = self._generate_module_file(idx, subworkflow_tasks)
            self._write_output_file(
                module_code,
                modules_dir.joinpath(f"tasks_{idx}.nf")
            )

        # Create the main Nextflow workflow script
        main_workflow_code = self._generate_main_workflow(sorted_tasks)
        self._write_output_file(main_workflow_code, output_folder.joinpath("workflow.nf"))

        # Create the README file
        self._write_readme_file(output_folder)

        return

    def _get_tasks_in_topological_order(self) -> List[Task]:
        """
        Sort the workflow tasks in topological order.

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

    def _split_into_subworkflows(self, sorted_tasks: List[Task], max_parents_threshold: int = 100) -> None:
        """
        Split the sorted tasks into chunks for separate files.

        Tasks with more than max_parents_threshold parents are placed in their own
        subworkflow to avoid generating overly long channel mixing code.

        :param sorted_tasks: The topologically sorted list of tasks.
        :type sorted_tasks: List[Task]
        :param max_parents_threshold: Tasks with more parents than this get their own subworkflow.
        :type max_parents_threshold: int
        """
        self.subworkflows = []
        current_chunk = []

        for task in sorted_tasks:
            num_parents = len(self._find_parents(task.task_id))

            # If task has many parents, put it in its own subworkflow
            if num_parents > max_parents_threshold:
                # First, save the current chunk if not empty
                if current_chunk:
                    self.subworkflows.append(current_chunk)
                    current_chunk = []
                # Add this task as its own subworkflow
                self.subworkflows.append([task])
            else:
                current_chunk.append(task)
                # If chunk is full, start a new one
                if len(current_chunk) >= self.max_tasks_per_subworkflow:
                    self.subworkflows.append(current_chunk)
                    current_chunk = []

        # Don't forget the last chunk
        if current_chunk:
            self.subworkflows.append(current_chunk)

    def _create_task_script(self, task: Task):
        """
        Generate the bash script for invoking a task.

        :param task: The task.
        :type task: Task
        """
        code = "#!/bin/bash\n\n"

        # Generate input spec
        input_spec = "'\\["
        for f in task.input_files:
            input_spec += f"\"{self.output_folder}/data/{f.file_id}\","
        input_spec = input_spec[:-1] + "\\]'"

        # Generate output spec
        output_spec = "'\\{"
        for f in task.output_files:
            output_spec += f"\"{self.output_folder.absolute()}/data/{f.file_id}\":{str(f.size)},"
        output_spec = output_spec[:-1] + "\\}'"

        code += f"{self.output_folder.absolute()}/bin/{task.program} "

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

    def _generate_module_file(self, module_idx: int, tasks: List[Task]) -> str:
        """
        Generate a module .nf file containing processes and an orchestrating function.

        Each module file contains:
        - Process definitions (private to the file)
        - A single exported function that orchestrates those processes

        Using a function instead of a workflow block allows passing Maps of channels
        between modules without Nextflow channel type casting issues.

        :param module_idx: The index of this module.
        :type module_idx: int
        :param tasks: The tasks in this module.
        :type tasks: List[Task]
        :return: The complete module file content.
        :rtype: str
        """
        code = f"// Module {module_idx} - Tasks file\n"
        code += "// Auto-generated by WfCommons NextflowSubworkflowTranslator\n\n"

        # Define pwd variable from params (modules can't access main workflow variables)
        code += "// Resolve working directory from params\n"
        code += "def pwd = params.pwd ? file(params.pwd).toAbsolutePath().toString() : null\n\n"

        # Generate process definitions for each task (private to this file)
        for task in tasks:
            code += self._generate_task_process(task)

        # Generate the function that orchestrates these tasks
        code += self._generate_module_function(module_idx, tasks)

        return code

    def _generate_module_function(self, module_idx: int, tasks: List[Task]) -> str:
        """
        Generate a Groovy function that orchestrates a set of tasks.

        Using a function instead of a workflow block allows passing Maps of channels
        between modules without Nextflow channel type casting issues.

        :param module_idx: The index of this module.
        :type module_idx: int
        :param tasks: The tasks in this module.
        :type tasks: List[Task]
        :return: The function code.
        :rtype: str
        """
        code = f"// Function to execute tasks in module {module_idx}\n"
        code += f"def run_module_{module_idx}(Map inputs) {{\n"
        code += "\tdef results = inputs.clone()\n\n"

        # Call each task's process
        for task in tasks:
            function_name = task.task_id.replace(".", "_")
            code += self._generate_task_call_in_function(task, function_name)

        code += "\treturn results\n"
        code += "}\n\n"

        return code

    def _generate_task_call_in_function(self, task: Task, function_name: str) -> str:
        """
        Generate the code to call a task's process within a module function.

        :param task: The task.
        :type task: Task
        :param function_name: The sanitized function name.
        :type function_name: str
        :return: The code.
        :rtype: str
        """
        code = ""

        if self._find_parents(task.task_id):
            # Mix input channels and call process
            code += f"\t// Task: {task.task_id}\n"
            code += f"\tdef {function_name}_input = Channel.empty()\n"
            for f in task.input_files:
                code += f"\t{function_name}_input = {function_name}_input.mix(results.{f.file_id})\n"
            code += f"\tdef {function_name}_output = {function_name}({function_name}_input.collect())\n"
        else:
            # Root task - no inputs needed
            code += f"\t// Task: {task.task_id} (root)\n"
            code += f"\tdef {function_name}_output = {function_name}()\n"

        # Update results map with outputs
        if self._find_children(task.task_id):
            counter = 0
            for f in task.output_files:
                code += f"\tresults.{f.file_id} = {function_name}_output.map{{it[{counter}]}}\n"
                counter += 1
        code += "\n"

        return code

    def _generate_task_process(self, task: Task) -> str:
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

    def _generate_main_workflow(self, sorted_tasks: List[Task]) -> str:
        """
        Generate the main workflow file that orchestrates all module functions.

        :param sorted_tasks: All tasks in topological order (used for Flowcept).
        :type sorted_tasks: List[Task]
        :return: The main workflow file content.
        :rtype: str
        """
        # Start with the template header
        code = """params.simulate = false
params.pwd = null
params.help = null
pwd = null

def printUsage(error_msg, exit_code) {
    def usage_string = \"\"\"
Usage: nextflow run workflow.nf --pwd /path/to/directory [--simulate] [--help]

    Required parameters:
      --pwd         Working directory (where the workflow.nf file is located)

    Optional parameters:
      --help        Show this message and exit.
      --simulate    Use a "sleep 1" for all tasks instead of the WfBench benchmark.
\"\"\"
    if (error_msg) {
        def RED = '\\u001B[31m'
        def RESET = '\\u001B[0m'
        System.err.println \"${RED}Error: ${RESET}\" + error_msg
    }
    System.err.println usage_string
    exit exit_code
}

def validateParams() {
    if (params.help) {
        printUsage(msg = \"\", exit_code=0)
    }
    if (params.pwd == null) {
        printUsage(msg = \"Missing required parameter: --pwd\", exit_code=1)
    }
    pwd = file(params.pwd).toAbsolutePath().toString()
    if (!file(pwd).exists()) {
        printUsage(msg = \"Directory not found: ${pwd}\", exit_code=1)
    }
}

// Call validation at the start
validateParams()

"""

        # Include module functions (one per line to avoid long strings)
        code += "// Include module functions\n"
        for idx in range(len(self.subworkflows)):
            code += f"include {{ run_module_{idx} }} from './modules/tasks_{idx}.nf'\n"
        code += "\n"

        # Add Flowcept process if enabled
        if self.workflow.workflow_id:
            code += self._generate_flowcept_code()

        # Generate the main workflow
        code += "workflow {\n"

        if self.workflow.workflow_id:
            code += "\tflowcept()\n"

        # Initialize empty channel map
        code += "\t// Initialize empty channel map\n"
        code += "\tdef ch_data = [:]\n\n"

        # Call each module function in sequence, passing the channel map
        code += "\t// Execute modules in sequence\n"
        for idx in range(len(self.subworkflows)):
            code += f"\tch_data = run_module_{idx}(ch_data)\n"

        code += "}\n"

        return code

    def _generate_flowcept_code(self) -> str:
        """
        Generate the Flowcept process code.

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

    def _write_readme_file(self, output_folder: pathlib.Path) -> None:
        """
        Write the README file.

        :param output_folder: The path of the output folder.
        :type output_folder: pathlib.Path
        """
        readme_file_path = output_folder.joinpath("README")
        with open(readme_file_path, "w") as out:
            out.write(f"Run the workflow in directory {str(output_folder)} using the following command:\n")
            out.write(f"\tnextflow run ./workflow.nf --pwd `pwd`\n\n")
            out.write(f"This workflow has been split into {len(self.subworkflows)} module file(s), ")
            out.write(f"each containing a maximum of {self.max_tasks_per_subworkflow} tasks.\n")
            out.write(f"\nModule files are located in the 'modules/' directory.\n")
