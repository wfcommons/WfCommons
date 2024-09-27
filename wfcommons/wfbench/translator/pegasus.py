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

from logging import Logger
from typing import Dict, Optional, Union

from .abstract_translator import Translator
from ...common import FileLink, Workflow


class PegasusTranslator(Translator):
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

        self.script = "import os\n" \
                      "from Pegasus.api import *\n\n\n" \
                      "def which(file):\n" \
                      "    for path in os.environ['PATH'].split(os.pathsep):\n" \
                      "        if os.path.exists(os.path.join(path, file)):\n" \
                      "            return os.path.join(path, file)\n" \
                      "    return None\n\n\n"
        self.parsed_tasks = []
        self.tasks_map = {}
        self.task_counter = 1

    def translate(self, output_file_path: pathlib.Path, tasks_priorities: Optional[Dict[str, int]] = None) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into a Pegasus workflow application.

        :param output_file_path: The path of the output file (e.g., workflow.py).
        :type output_file_path: pathlib.Path
        :param tasks_priorities: Priorities to be assigned to tasks.
        :type tasks_priorities: Optional[Dict[str, int]]
        """
        # overall workflow
        self.script += f"wf = Workflow('{self.workflow.name}', infer_dependencies=True)\n" \
                       "tc = TransformationCatalog()\n" \
                       "rc = ReplicaCatalog()\n\n"
        self.script += "task_output_files = {}\n\n"

        # transformation catalog
        transformations = []

        # cpu-benchmark
        self.script += "t_cpu_benchmark = Transformation('cpu-benchmark', site='local',\n" \
                       "pfn=os.getcwd() + '/cpu-benchmark', is_stageable=True)\n" \
                       "tc.add_transformations(t_cpu_benchmark)\n\n"

        # tasks' programs
        for task in self.tasks.values():
            if task.name not in transformations:
                transformations.append(task.name)
                self.script += f"transformation_path = which('{task.program}')\n" \
                               "if transformation_path is None:\n" \
                               f"    raise RuntimeError('Unable to find {task.program}')\n" \
                               f"transformation = Transformation('{task.name}', site='local',\n" \
                               f"                                pfn=transformation_path,\n" \
                               "                                is_stageable=True)\n" \
                               "transformation.add_env(PATH='/usr/bin:/bin:.')\n" \
                               "transformation.add_profiles(Namespace.CONDOR, 'request_disk', '10')\n" \
                               "transformation.add_requirement(t_cpu_benchmark)\n" \
                               "tc.add_transformations(transformation)\n\n"

        # adding tasks
        for task_name in self.root_task_names:
            self._add_task(task_name, tasks_priorities=tasks_priorities)
            # input file
            task = self.tasks[task_name]
            for file in task.files:
                if file.link == FileLink.INPUT:
                    self.script += f"in_file_{self.task_counter} = File('{file.file_id}')\n"
                    self.script += f"rc.add_replica('local', '{file.file_id}', 'file://' + os.getcwd() + " \
                                   f"'/data/{file.file_id}')\n"
                    self.script += f"{self.tasks_map[task_name]}.add_inputs(in_file_{self.task_counter})\n" \
                                   f"print('Using input data: ' + os.getcwd() + '/data/{file.file_id}')\n"

        self.script += "\n"

        # write out the workflow
        self.script += "wf.add_replica_catalog(rc)\n" \
                       "wf.add_transformation_catalog(tc)\n" \
                       f"wf.write('{self.workflow.name}-benchmark-workflow.yml')\n"

        # write script to file
        self._write_output_file(self.script, output_file_path)

    def _add_task(self, task_name: str, parent_task: Optional[str] = None, tasks_priorities: Optional[Dict[str, int]] = None) -> None:
        """
        Add a task and its dependencies to the workflow.

        :param task_name: name of the task
        :type task_name: str
        :param parent_task: name of the parent task
        :type parent_task: Optional[str]
        :param tasks_priorities: Priorities to be assigned to tasks.
        :type tasks_priorities: Optional[Dict[str, int]]
        """
        if task_name not in self.parsed_tasks:
            task = self.tasks[task_name]
            job_name = f"job_{self.task_counter}"
            self.script += f"{job_name} = Job('{task.name}', _id='{task_name}')\n" \
                f"task_output_files.setdefault('{job_name}', [])\n"

            # task priority
            if tasks_priorities and task.name in tasks_priorities:
                self.script += f"{job_name}.add_condor_profile(priority='{tasks_priorities[task.name]}')\n"

            # find children
            children = self.task_children[task_name]

            # output file
            for file in task.files:
                if file.link == FileLink.OUTPUT:
                    out_file = file.file_id
                    # task.args.append(f"--out={out_file}")
                    stage_out = "True" if len(children) == 0 else "False"
                    self.script += f"out_file_{self.task_counter} = File('{out_file}')\n" \
                        f"task_output_files['{job_name}'].append(out_file_{self.task_counter})\n" \
                        f"{job_name}.add_outputs(out_file_{self.task_counter}, " \
                        f"stage_out={stage_out}, register_replica={stage_out})\n"

            # arguments
            args = []
            for a in task.args:
                a = a.replace("'", "\"") if "--out" not in a else a.replace("{", "\"{").replace("}", "}\"").replace("'", "\\\\\"").replace(": ", ":")
                args.append(a)
            args = ", ".join(f"'{a}'" for a in args)
            self.script += f"{job_name}.add_args({args})\n"

            self.script += f"wf.add_jobs({job_name})\n\n"
            self.task_counter += 1
            self.parsed_tasks.append(task_name)
            self.tasks_map[task_name] = job_name

            for child_task_name in children:
                self._add_task(child_task_name, job_name, tasks_priorities)

        if parent_task:
            self.script += f"if '{parent_task}' in task_output_files:\n" \
                        f"    for f in task_output_files['{parent_task}']:\n" \
                        f"        {self.tasks_map[task_name]}.add_inputs(f)\n" \
                        f"wf.add_dependency({self.tasks_map[task_name]}, parents=[{parent_task}])\n\n"
