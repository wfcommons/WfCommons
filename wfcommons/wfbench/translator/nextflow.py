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
import re

from collections import defaultdict
from logging import Logger
from typing import Dict, List, Optional, Union, MutableSet

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

        self.script = """
String extractIDfromFile(Path filepath, String task_name) {
  String filename = filepath as String
  filename = filename[filename.lastIndexOf('/')+1..-1]
  
  int first_id_start = filename.indexOf('0')
  if (first_id_start < 0) return
  String first_task_name = filename.substring(0, first_id_start - 1)
  if (first_task_name == task_name) return filename.substring(first_id_start, first_id_start + 8)

  if (first_id_start + 9 >= filename.length()) return

  int second_id_start = filename.indexOf('0', first_id_start + 8)
  if (second_id_start < 0) return
  String second_task_name = filename.substring(first_id_start + 9, second_id_start - 1)
  if (second_task_name == task_name) return filename.substring(second_id_start, second_id_start + 8)
}

"""

    def translate(self, output_file_path: pathlib.Path, input_file_directory: pathlib.Path = None) -> None:
        """
        Translate a workflow benchmark description(WfFormat) into a Nextflow workflow application.

        : param output_file_path: The name of the output file(e.g., workflow.py).
        : type output_file_path: pathlib.Path
        """

        self.task_inputs: Dict[str, List[str]] = {}
        self.task_outputs: Dict[str, List[str]] = {}

        self.missing_files = []

        # determine the abstract tasks and their abstract parents and children
        self.abstract_tasks = defaultdict(list)
        self.abstract_children: Dict[str, MutableSet[str]] = defaultdict(set)
        self.abstract_parents: Dict[str, MutableSet[str]] = defaultdict(set)
        self._determine_abstract_tasks()

        all_inputs = []
        all_outputs = []
        for inputs in self.task_inputs.values():
            all_inputs.extend(inputs)
        for outputs in self.task_outputs.values():
            all_outputs.extend(outputs)
        self.workflow_inputs = [i for i in all_inputs if i not in all_outputs]

        self.task_outstreams: Dict = {}
        for abstract_task_name, physical_tasks in self.abstract_tasks.items():
            self._add_output_map(abstract_task_name, physical_tasks)
        self.task_written: Dict[str, bool] = {}
        for abstract_task_name, physical_tasks in self.abstract_tasks.items():
            self.task_written[abstract_task_name] = False
            self._add_abstract_task(abstract_task_name, physical_tasks)

        self.script += "workflow {\n"
        workflow_inputs_path: str = input_file_directory.absolute() if input_file_directory else pathlib.Path("./")
        for missing_file in self.missing_files:
            self.script += f"  file(\"{workflow_inputs_path.joinpath(missing_file)}\").text = \"?\"\n"
        self.script += "\n"

        self.script += f"  workflow_inputs = Channel.fromPath(\"{workflow_inputs_path.joinpath('*')}\")\n"
        self.script += "\n"

        for abstract_task_name, physical_tasks in self.abstract_tasks.items():
            if not self.task_written[abstract_task_name]:
                self._write_abstract_task(abstract_task_name, physical_tasks)

        self.script += "}\n"

        # write script to file
        self._write_output_file(self.script, output_file_path)

    def _add_output_map(self, abstract_task_name: str, physical_tasks: List[Task]) -> None:
        self.script += f"def {self.valid_task_name(abstract_task_name)}_out_args = [\n"
        for ptask in physical_tasks:
            out_args = [arg for arg in ptask.args if arg.startswith("--out ")]
            assert (len(out_args) == 1)
            out_arg = out_args[0].replace("--out {", "").replace("}", "").replace("'", "\\\\\\\"").replace(": ", ":")
            self.script += f"  \"{ptask.task_id}\": \"{out_arg}\",\n"
        self.script += "]\n\n"

    def valid_task_name(self, original_task_name: str) -> str:
        return original_task_name.replace('-', '_')

    def _write_abstract_task(self, abstract_task_name: str, physical_tasks: List[Task]) -> None:
        parents = self.abstract_parents[abstract_task_name]
        for parent in parents:
            if parent == abstract_task_name:
                continue
            if not self.task_written[parent]:
                self._write_abstract_task(parent, self.abstract_tasks[parent])

        # determining the channel of all raw inputs (outputs of other tasks)
        input_channels = [
            f"{self.valid_task_name(parent)}_FOR_{self.valid_task_name(abstract_task_name)}" for parent in parents]
        example_task = physical_tasks[0]
        for input_file in self.task_inputs[example_task.name]:
            if input_file in self.workflow_inputs or input_file in self.missing_files:
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
        self.script += f"  {self.valid_task_name(abstract_task_name)}_in = {inputs_channel}.flatten().map{{\n"
        self.script += f"    String id = extractIDfromFile(it, \"{abstract_task_name}\")\n"
        self.script += f"    if (id) return [id, it]\n"
        self.script += "  }.groupTuple()\n"

        # creating one output channel for each child
        out_channels = [f"{self.valid_task_name(abstract_task_name)}_FOR_{self.valid_task_name(name)}" for name,
                        _ in self.task_outstreams[abstract_task_name]]
        self.script += "  "
        if len(out_channels) == 1:
            self.script += f"{out_channels[0]} = "
        elif len(out_channels) > 1:
            self.script += f"({', '.join(out_channels)}) = "
        self.script += f"task_{self.valid_task_name(abstract_task_name)}({self.valid_task_name(abstract_task_name)}_in)\n\n"

        self.task_written[abstract_task_name] = True

    def _determine_abstract_tasks(self) -> None:
        """
        Determines the abstract tasks that will be used for the nextflow definition of the workflow.
        """

        abstract_input_amounts = defaultdict(list)
        abstract_output_amounts = defaultdict(list)

        tasks_with_iterations = set()
        for task in self.tasks.values():
            abstract_task: str = task.category
            self.abstract_tasks[abstract_task].append(task)
            self.task_inputs[task.name] = [file.name for file in task.files if file.link == FileLink.INPUT]
            self.task_outputs[task.name] = [file.name for file in task.files if file.link == FileLink.OUTPUT]
            self.missing_files.extend([filename for filename in self.task_inputs[task.name]
                                       if filename in self.task_outputs[task.name]])
            abstract_input_amounts[abstract_task].append(len(self.task_inputs[task.name]))
            abstract_output_amounts[abstract_task].append(len(self.task_outputs[task.name]))

            for child in self.task_children[task.name]:
                abstract_child: str = self.tasks[child].category
                if abstract_child == abstract_task:
                    tasks_with_iterations.add(abstract_task)
                self.abstract_children[abstract_task].add(abstract_child)
            for parent in self.task_parents[task.name]:
                abstract_parent: str = self.tasks[parent].category
                if abstract_parent == abstract_task:
                    tasks_with_iterations.add(abstract_task)
                self.abstract_parents[abstract_task].add(abstract_parent)

        if tasks_with_iterations:
            error_msg: str = "Iterations are not supported by Nextflow. "\
                "Thus, this workflow has no Nextflow translation since the following "\
                f"{'task uses' if len(tasks_with_iterations) == 1 else 'tasks use'} "\
                "iterations:\n    "
            error_msg += ", ".join(tasks_with_iterations)
            raise RuntimeError(error_msg)

    def _add_abstract_task(self, abstract_task_name: str, physical_tasks: List[Task]) -> None:
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
        memory = f"{max(memory_values)}.KB" if memory_values else "8.GB"

        self.script += f"process task_{self.valid_task_name(abstract_task_name)}" + " {\n"
        if cores:
            self.script += f"  cpus {cores}\n"
        if memory:
            self.script += f"  memory {memory}\n"
        self.script += "  input:\n"
        self.script += f"    tuple val( id ), path( \"*\" )\n"

        example_task = physical_tasks[0]

        self.script += "  output:\n"
        outstream_contains_infiles = defaultdict(bool)
        for ptask in physical_tasks:
            outfile_regex = f"{abstract_task_name}_{ptask.task_id}_(.*)_(\d{{8}})(.*)"
            for outfile in self.task_outputs[ptask.name]:
                outfile_search = re.search(outfile_regex, outfile)
                if outfile_search:
                    is_infile_aswell = outfile in self.task_inputs[ptask.name]
                    outstream = (outfile_search.group(1), outfile_search.group(3))
                    outstream_contains_infiles[outstream] |= is_infile_aswell
        outstreams = list(outstream_contains_infiles.keys())
        self.task_outstreams[abstract_task_name] = outstreams
        for outstream in outstreams:
            self.script += f"    path( \"{abstract_task_name}_${{id}}_{outstream[0]}_????????{outstream[1]}\""
            if outstream_contains_infiles[outstream]:
                self.script += ", includeInputs: true"
            self.script += " )\n"

        cmd = pathlib.Path(example_task.program).name
        for a in example_task.args:
            if a == f"{abstract_task_name}_{example_task.task_id}":
                cmd += f" {abstract_task_name}_${{id}}"
                continue
            cmd += ' '
            if a.startswith("--out"):
                cmd += "--out \"{${" + self.valid_task_name(abstract_task_name) + "_out_args.get(id)}}\" \\$inputs"
                break
            else:
                a = a.replace(f"{abstract_task_name}_{example_task.task_id}", f"{abstract_task_name}_${{id}}")
                cmd += a.replace("'", "\"")

        self.script += "\n"
        self.script += "  script:\n"
        self.script += "  \"\"\"\n"
        self.script += f"  inputs=*{abstract_task_name}_${{id}}*\n"
        self.script += f"  {cmd}"
        self.script += '\n'
        self.script += "  \"\"\"\n"
        self.script += "}\n"
