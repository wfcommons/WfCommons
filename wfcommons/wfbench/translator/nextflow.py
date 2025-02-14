#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import json
import ast

from collections import defaultdict
from math import ceil
from logging import Logger
from typing import Dict, List, Optional, Union, MutableSet

from pyparsing import empty

from .abstract_translator import Translator
from ...common import File, FileLink, Workflow
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

        self.script = ""

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description(WfFormat) into a Nextflow workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """

        # Create the output folder
        output_folder.mkdir(parents=True)

        # Create a topological order of the tasks
        sorted_tasks = self._get_tasks_in_topological_order()

        # Output the code for each task
        for task in sorted_tasks:
            self.script += self._generate_task_code(task)

        # Output the code for the workflow
        self.script += self._generate_workflow_code(sorted_tasks)

        # Output the code to the workflow file
        self._write_output_file(self.script, output_folder.joinpath("workflow.nf"))

        # Create the README file
        self._write_readme_file(output_folder)

        # Create additional files
        self._copy_binary_files(output_folder)
        self._generate_input_files(output_folder)

        return


        # # determine the abstract tasks and their abstract parents and children
        # self.abstract_tasks = defaultdict(list)
        # self.abstract_parents: Dict[str, MutableSet[str]] = defaultdict(set)
        # self._determine_abstract_relations()
        #
        # self.task_inputs: Dict[str, List[File]] = {}
        # self.task_outputs: Dict[str, List[File]] = {}
        # self.task_input_amounts: Dict[str, int] = {}
        # self._determine_input_output()
        #
        #
        #
        # self.workflow_inputs = set().union(*self.task_inputs.values()).difference(*self.task_outputs.values())
        #
        # self._create_file_task_mappings(output_folder)
        #
        # for abstract_task_name, physical_tasks in self.abstract_tasks.items():
        #     self._create_task_args_map(output_folder, abstract_task_name, physical_tasks)
        # self.script += "\n\n"
        #
        # self.task_written: Dict[str, bool] = {}
        # for abstract_task_name, physical_tasks in self.abstract_tasks.items():
        #     self.task_written[abstract_task_name] = False
        #     self._add_abstract_task_definition(abstract_task_name, physical_tasks)
        #
        # self.script += "workflow {\n"
        # self.script += "  workflow_inputs = Channel.fromPath(\"${params.indir}/*\")\n"
        # self.script += "\n"
        #
        # for abstract_task_name, physical_tasks in self.abstract_tasks.items():
        #     if not self.task_written[abstract_task_name]:
        #         self._add_call_to_abstract_task(abstract_task_name, physical_tasks)
        # self.script += "}\n"
        #
        # self._write_output_file(self.script, output_folder.joinpath("workflow.nf"))
        #
        # self._write_readme_file(output_folder)
        #
        # # additional files
        # self._copy_binary_files(output_folder)
        # self._generate_input_files(output_folder)

    def _get_tasks_in_topological_order(self) -> List[Task]:
        levels = {0: self._find_root_tasks()}
        sorted_tasks: List[Task] = levels[0]
        current_level = 1
        while (True):
            # print(f"Dealing with level {current_level}")
            tasks_in_current_level = []
            all_children = [self._find_children(p.task_id) for p in levels[current_level-1]]
            all_children = [item for sublist in all_children for item in sublist]
            print(all_children)
            if not all_children:
                # print("NO MORE CHILDREN - DONE")
                break
            for potential_task in all_children:
                # print(f"Looking at child {potential_task.task_id}")
                if all(parent in sorted_tasks for parent in self._find_parents(potential_task.task_id)):
                    tasks_in_current_level.append(potential_task)
            levels[current_level] = tasks_in_current_level
            sorted_tasks.extend(tasks_in_current_level)
            current_level += 1
        # for level in levels:
            # print(f"Level {level}: {[t.task_id for t in levels[level]]}")
        return sorted_tasks


    def _generate_task_code(self, task: Task) -> str:
        code = f"process {task.task_id}()" + "{\n"
        code += f"\tinput:\n"
        if self._find_parents(task.task_id):
            for f in task.input_files:
                code += f"\t\tval {f.file_id}\n"
        code += "\n"

        code += f"\toutput:\n"
        for f in task.output_files:
            code += f"\t\tval {f.file_id}\n"
        code += "\n"

        code += "\tscript:\n"
        code += "\t\t\"\"\"\n"

        code += "\t\t\"${params.pwd}/bin/" + task.program + "\""
        for a in task.args:
            if "--output-files" in a:
                flag, output_files_dict = a.split(" ", 1)
                output_files_dict = {str("\"${params.pwd}/data/" + key): value for key, value in
                                     ast.literal_eval(output_files_dict).items()}
                a = f"{flag} '{json.dumps(output_files_dict)}'"
            elif "--input-files" in a:
                flag, input_files_arr = a.split(" ", 1)
                input_files_arr = [str("\"${params.pwd}/data/" + file) for file in
                                   ast.literal_eval(input_files_arr)]
                a = f"{flag} '{json.dumps(input_files_arr)}'"

            code += " " + a
        code += " ".join(task.args) + "\n"

        code += "\t\t\"\"\"\n"
        code += "}\n\n"
        print(code)
        return code

    def _generate_workflow_code(self, sorted_tasks: List[Task]) -> str:
        code = "workflow {\n"
        for task in sorted_tasks:
            code += self._generate_task_invocation_code(task)
        code += "}\n"
        return code

    def _generate_task_invocation_code(self, task: Task) -> str:

        # Figure out task output values
        if task.output_files and self._find_children(task.task_id):
            output_values = "(" + ",".join(["SSS"+f.file_id for f in task.output_files]) + ")"
        else:
            output_values = "_"

        # Figure out task input values
        input_values = ",".join([f.file_id for f in task.input_files])

        code = output_values + " = " + task.task_id + "(" + input_values + ")\n\n"

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
            out.write(f"nextflow run ./workflow.nf --indir .\n")


    def _is_resource_arg(self, arg: str) -> bool:
        return arg.startswith("--percent-cpu") or arg.startswith("--mem") \
            or arg.startswith("--cpu-work") or arg.startswith("--gpu-work")

    def _determine_abstract_relations(self) -> None:
        """
        Determines the abstract tasks that will be used for the nextflow definition of the workflow.
        """

        for task in self.tasks.values():
            abstract_task: str = task.name
            self.abstract_tasks[abstract_task].append(task)

            for parent in self.task_parents[task.task_id]:
                abstract_parent: str = self.tasks[parent].name
                self.abstract_parents[abstract_task].add(abstract_parent)

        tasks_with_iterations = set()
        for abstract_task in self.abstract_tasks:
            for abstract_parent in self.abstract_parents[abstract_task]:
                if abstract_parent == abstract_task:
                    tasks_with_iterations.add(abstract_task)
        if tasks_with_iterations:
            error_msg: str = "Iterations are not supported by Nextflow. "\
                "Thus, this workflow has no Nextflow translation since the following "\
                f"{'task uses' if len(tasks_with_iterations) == 1 else 'tasks use'} "\
                "iterations:\n    "
            error_msg += ", ".join(tasks_with_iterations)
            raise RuntimeError(error_msg)

    def _determine_input_output(self) -> None:
        """
        Determines the inputs and outputs for the physical and abstract tasks.
        """
        for task in self.tasks.values():
            self.task_inputs[task.task_id] = [file for file in task.input_files]
            self.task_outputs[task.task_id] = [file for file in task.output_files]

        self.script += "// define amount of input files for abstracts tasks where the amount is not constant\n"
        for abstract_task_name, physical_tasks in self.abstract_tasks.items():
            input_amounts = {task.task_id: len(self.task_inputs[task.task_id]) for task in physical_tasks}
            if (max(input_amounts.values()) == min(input_amounts.values())):
                # all physical tasks have the same amount of inputs
                self.task_input_amounts[abstract_task_name] = max(input_amounts.values())
            else:
                # define amount of inputs for each physical task
                self.script += f"def {self.valid_task_name(abstract_task_name)}_input_amounts = [\n"
                for k, v in input_amounts.items():
                    self.script += f"  \"{k}\": {v},\n"
                self.script += "]\n"
                self.task_input_amounts[abstract_task_name] = None
        self.script += "\n"

    def _create_file_task_mappings(self, output_folder: pathlib.Path) -> None:
        file_task_map = defaultdict(list)
        for task_name, task_input_files in self.task_inputs.items():
            task = self.tasks[task_name]
            for file in task_input_files:
                file_task_map[str(output_folder) + "/data/" + file.file_id].append([task.name, task.task_id])
        self._write_map_file(file_task_map, "file_inputs", output_folder)

    def _write_map_file(self, map_dict: Dict, map_name: str, output_folder: pathlib.Path) -> None:
        path = output_folder.joinpath(f"{map_name}.json")
        with open(path, "w") as f:
            f.write(json.dumps(map_dict, indent=4))
        self.script += f"{map_name} = jsonSlurper.parseText(file(\"${{projectDir}}/{map_name}.json\").text)\n"

    def _create_task_args_map(self, output_folder: pathlib.Path, abstract_task_name: str, physical_tasks: List[Task]) -> None:
        map_name = f"{self.valid_task_name(abstract_task_name)}_args"
        task_args_map = {}
        for ptask in physical_tasks:
            out_file_sizes = {str(output_folder) + "/data/" + file.file_id: file.size for file in self.task_outputs[ptask.task_id]}
            out_arg = str(out_file_sizes).replace("{", "").replace("}", "").replace("'", "\\\"").replace(": ", ":")
            task_args_map[ptask.task_id] = {
                "out": out_arg,
                "resources": " ".join((arg for arg in ptask.args if self._is_resource_arg(arg)))
            }
        self._write_map_file(task_args_map, map_name, output_folder)

    def valid_task_name(self, original_task_name: str) -> str:
        return original_task_name.replace('-', '_')

    def _add_abstract_task_definition(self, abstract_task_name: str, physical_tasks: List[Task]) -> None:
        """
        Add an abstract task to the workflow considering its physical tasks.

        : param abstract_task_name: the name of the abstract task
        : type abstract_task_name: str
        : param physical_tasks: a list of physical tasks for this abstract tasks
        : type physical_tasks: List[Task]
        """

        cores_values = [task.cores for task in physical_tasks if task.cores is not None]
        if len(cores_values) == 0:
            cores = None
        else:
            cores = int(max(cores_values))
        memory_values = [task.memory for task in physical_tasks if task.memory is not None]
        if len(memory_values) == 0:
            memory = None
        else:
            memory = max(memory_values) * 1.05

        # creating the command for the abstract task using the first physical task as a template
        example_task = physical_tasks[0]
        print(example_task)
        cmd = pathlib.Path(example_task.program).name
        resource_args_done = False
        for a in example_task.args:
            print("ARG = " , a)
            print(f"{abstract_task_name}_{example_task.task_id}")
            if a == f"{abstract_task_name}_{example_task.task_id}":
                print("GOT THE ABSTRACT TASK NAME")
                cmd += f" {abstract_task_name}_${{id}}"
                continue
            if self._is_resource_arg(a):
                if resource_args_done:
                    continue
                cmd += " ${" + self.valid_task_name(abstract_task_name) + "_args.get(id).get(\"resources\")}"
                resource_args_done = True
            elif a.startswith("--out"):
                cmd += " --out \"{${" + self.valid_task_name(abstract_task_name) + "_args.get(id).get(\"out\")}}\""
                cmd += " \\$inputs"
                break
            else:
                print("===> ", a)
                a = a.replace(f"{abstract_task_name}_{example_task.task_id}", f"{abstract_task_name}_${{id}}")
                cmd += " " + a.replace("'", "\"")

        # creating the abstract task
        self.script += f"process task_{self.valid_task_name(abstract_task_name)}" + " {\n"
        if cores:
            self.script += f"  cpus {cores}\n"
        if memory:
            self.script += f"  memory '{human_readable_memory(memory)}'\n"
        self.script += "  input:\n"
        self.script += "    tuple val( id ), path( \"*\" )\n"
        self.script += f"  output:\n    path( \"{self.valid_task_name(abstract_task_name)}_????????_outfile_????*\" )\n"
        self.script += "  script:\n"
        self.script += "  \"\"\"\n"
        self.script += "  inputs=\\$(find . -maxdepth 1 -name \\\"workflow_infile_*\\\" -or -name \\\"*_outfile_0*\\\")\n"
        self.script += f"  {cmd}\n"
        self.script += "  \"\"\"\n"
        self.script += "}\n"

    def _add_call_to_abstract_task(self, abstract_task_name: str, physical_tasks: List[Task]) -> None:
        parents = self.abstract_parents[abstract_task_name]
        for parent in parents:
            if parent == abstract_task_name:
                raise RuntimeError("Iterations are not supported by Nextflow.")
            if not self.task_written[parent]:
                self._add_call_to_abstract_task(parent, self.abstract_tasks[parent])

        # determining the channel of all raw inputs (outputs of other tasks)
        input_channels = [
            f"{self.valid_task_name(parent)}_out" for parent in parents]
        example_task = physical_tasks[0]

        for input_file in self.task_inputs[example_task.task_id]:
            if input_file in self.workflow_inputs:
                input_channels.append("workflow_inputs")
                break

        if len(input_channels) == 1:
            inputs_channel = input_channels[0]
        elif len(input_channels) != 0:
            # concatenating all the input channels into one big channel
            one_channel = input_channels.pop()
            inputs_channel = f"concatenated_FOR_{self.valid_task_name(abstract_task_name)}"
            self.script += f"  {inputs_channel} = {one_channel}.concat({', '.join(input_channels)})\n"
        else:
            raise RuntimeError(f"The abstract task {abstract_task_name} has no inputs.")

        # creating the input channel for this abstract task by grouping the outputs from the parents by id
        self.script += f"  {self.valid_task_name(abstract_task_name)}_in = {inputs_channel}.flatten().flatMap{{\n"
        self.script += f"    List<String> ids = extractTaskIDforFile(it, \"{abstract_task_name}\")\n"
        self.script += "    def pairs = new ArrayList()\n"
        self.script += "    for (id : ids) pairs.add([id, it])\n"
        self.script += "    return pairs\n"
        self.script += "  }"

        if self.task_input_amounts[abstract_task_name]:
            self.script += f".groupTuple(size: {self.task_input_amounts[abstract_task_name]})\n"
        else:
            self.script += f".map {{ id, file -> tuple( groupKey(id, {self.valid_task_name(abstract_task_name)}_input_amounts[id]), file ) }}\n"
            self.script += "  .groupTuple()\n"

        self.script += f"  {self.valid_task_name(abstract_task_name)}_out = task_"
        self.script += f"{self.valid_task_name(abstract_task_name)}({self.valid_task_name(abstract_task_name)}_in)\n\n"

        self.task_written[abstract_task_name] = True


def human_readable_memory(mem_bytes: int) -> str:
    idx = 0
    memory = mem_bytes
    memory_units = ["B", "KB", "MB", "GB", "TB"]
    while memory > 4096 and idx < len(memory_units) - 1:
        memory /= 1024
        idx += 1
    memory = ceil(memory * 100) / 100  # ensure that it is an upper bound
    return f"{memory:.2f} {memory_units[idx]}"
