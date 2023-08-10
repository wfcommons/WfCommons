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
import json

from collections import defaultdict
from math import ceil
from logging import Logger
from typing import Dict, List, Optional, Union, MutableSet

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

        self.script = """
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

List<String> extractTaskIDforFile(Path filepath, String task_name) {
  String filename = filepath as String
  filename = filename[filename.lastIndexOf('/')+1..-1]

  List<String> ids_for_file = new ArrayList<String>()
  for (destination : file_inputs[filename]) {
    def destination_task_name = destination[0]
    def destination_task_id = destination[1]
    if (destination_task_name == task_name)
      ids_for_file.add(destination_task_id)
  }
  return ids_for_file
}

"""

    def translate(self, output_file_path: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description(WfFormat) into a Nextflow workflow application.

        : param output_file_path: The name of the output file(e.g., workflow.py).
        : type output_file_path: pathlib.Path
        """

        # determine the abstract tasks and their abstract parents and children
        self.abstract_tasks = defaultdict(list)
        self.abstract_parents: Dict[str, MutableSet[str]] = defaultdict(set)
        self._determine_abstract_relations()

        self.task_inputs: Dict[str, List[File]] = {}
        self.task_outputs: Dict[str, List[File]] = {}
        self.task_input_amounts: Dict[str, int] = {}
        self._determine_input_output()

        self.workflow_inputs = set().union(*self.task_inputs.values()).difference(*self.task_outputs.values())

        self._create_file_task_mappings(output_file_path)

        for abstract_task_name, physical_tasks in self.abstract_tasks.items():
            self._create_task_output_map(output_file_path, abstract_task_name, physical_tasks)
        self.script += "\n\n"

        self.task_written: Dict[str, bool] = {}
        for abstract_task_name, physical_tasks in self.abstract_tasks.items():
            self.task_written[abstract_task_name] = False
            self._add_abstract_task_definition(abstract_task_name, physical_tasks)

        self.script += "workflow {\n"
        self.script += "  workflow_inputs = Channel.fromPath(\"${params.indir}/*\")\n"
        self.script += "\n"

        for abstract_task_name, physical_tasks in self.abstract_tasks.items():
            if not self.task_written[abstract_task_name]:
                self._add_call_to_abstract_task(abstract_task_name, physical_tasks)
        self.script += "}\n"

        self._write_output_file(self.script, output_file_path)

    def _determine_abstract_relations(self) -> None:
        """
        Determines the abstract tasks that will be used for the nextflow definition of the workflow.
        """

        for task in self.tasks.values():
            abstract_task: str = task.category
            self.abstract_tasks[abstract_task].append(task)

            for parent in self.task_parents[task.name]:
                abstract_parent: str = self.tasks[parent].category
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
            self.task_inputs[task.name] = [file for file in task.files if file.link == FileLink.INPUT]
            self.task_outputs[task.name] = [file for file in task.files if file.link == FileLink.OUTPUT]

        self.script += "// define amount of input files for abstracts tasks where the amount is not constant\n"
        for abstract_task_name, physical_tasks in self.abstract_tasks.items():
            input_amounts = {task.task_id: len(self.task_inputs[task.name]) for task in physical_tasks}
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

    def _create_file_task_mappings(self, output_file_path: pathlib.Path) -> None:
        file_task_map = defaultdict(list)
        for task_name, task_input_files in self.task_inputs.items():
            task = self.tasks[task_name]
            for file in task_input_files:
                file_task_map[file.name].append([task.category, task.task_id])
        self._write_map_file(file_task_map, "file_inputs", output_file_path)

    def _write_map_file(self, map_dict: Dict, map_name: str, output_file_path: pathlib.Path) -> None:
        path = output_file_path.parent.joinpath(f"{map_name}.json")
        with open(path, "w") as f:
            f.write(json.dumps(map_dict, indent=4))
        self.script += f"{map_name} = jsonSlurper.parseText(file(\"${{params.indir}}/{map_name}.json\").text)\n"

    def _create_task_output_map(self, output_file_path: pathlib.Path, abstract_task_name: str, physical_tasks: List[Task]) -> None:
        map_name = f"{self.valid_task_name(abstract_task_name)}_out_args"
        task_output_map = {}
        for ptask in physical_tasks:
            out_file_sizes = {file.name: file.size for file in self.task_outputs[ptask.name]}
            out_arg = str(out_file_sizes).replace("{", "").replace("}", "").replace("'", "\\\"").replace(": ", ":")
            task_output_map[ptask.task_id] = out_arg
        self._write_map_file(task_output_map, map_name, output_file_path)

    def valid_task_name(self, original_task_name: str) -> str:
        return original_task_name.replace('-', '_')

    def _add_abstract_task_definition(self, abstract_task_name: str, physical_tasks: List[Task]) -> None:
        """
        Add an abstract task to the workflow considering it's physical tasks.

        : param abstract_task_name: the name of the abstract task
        : type abstract_task_name: str
        : param physical_tasks: a list of physical tasks for this abstract tasks
        : type physical_tasks: List[Task]
        """

        cores_values = [task.cores for task in physical_tasks if task.cores]
        cores = max(cores_values) if cores_values else None
        memory_values = [task.memory for task in physical_tasks if task.memory]
        memory = max(memory_values)*1024 if memory_values else 8 * 1024**3  # default to 8.GB

        example_task = physical_tasks[0]

        cmd = pathlib.Path(example_task.program).name
        for a in example_task.args:
            if a == f"{abstract_task_name}_{example_task.task_id}":
                cmd += f" {abstract_task_name}_${{id}}"
                continue
            cmd += ' '
            if a.startswith("--out"):
                cmd += "--out \"{${" + self.valid_task_name(abstract_task_name) + "_out_args.get(id)}}\""
                cmd += " \\$inputs"
                break
            else:
                a = a.replace(f"{abstract_task_name}_{example_task.task_id}", f"{abstract_task_name}_${{id}}")
                cmd += a.replace("'", "\"")

        self.script += f"process task_{self.valid_task_name(abstract_task_name)}" + " {\n"
        if cores:
            self.script += f"  cpus {cores}\n"
        if memory:
            self.script += f"  memory {human_readable_memory(memory)}\n"
        self.script += "  input:\n"
        self.script += f"    tuple val( id ), path( \"*\" )\n"

        self.script += f"  output:\n    path( \"{self.valid_task_name(abstract_task_name)}_????????_outfile_????*\" )\n"
        self.script += "  script:\n"
        self.script += "  \"\"\"\n"
        self.script += f"  inputs=\\$(find . -maxdepth 1 -name \\\"workflow_infile_*\\\" -or -name \\\"*_outfile_0*\\\")\n"
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

        for input_file in self.task_inputs[example_task.name]:
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
        self.script += f"    def pairs = new ArrayList()\n"
        self.script += f"    for (id : ids) pairs.add([id, it])\n"
        self.script += f"    return pairs\n"
        self.script += "  }"

        if self.task_input_amounts[abstract_task_name]:
            self.script += f".groupTuple(size: {self.task_input_amounts[abstract_task_name]})\n"
        else:
            self.script += f".map {{ id, file -> tuple( groupKey(id, {self.valid_task_name(abstract_task_name)}_input_amounts[id]), file ) }}\n"
            self.script += f"  .groupTuple()\n"

        self.script += f"  {self.valid_task_name(abstract_task_name)}_out = task_"
        self.script += f"{self.valid_task_name(abstract_task_name)}({self.valid_task_name(abstract_task_name)}_in)\n\n"

        self.task_written[abstract_task_name] = True


def human_readable_memory(mem_bytes: int) -> str:
    idx = 0
    memory = mem_bytes
    l = ["B", "KB", "MB", "GB", "TB"]
    while memory > 4096 and idx < len(l)-1:
        memory /= 1024
        idx += 1
    memory = ceil(memory * 100) / 100  # ensure that it is an upper bound
    return f"{memory:.2f}{l[idx]}"
