#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import uuid

from logging import Logger
from typing import Optional

from .abstract_translator import Translator


class PegasusTranslator(Translator):
    """A WfFormat parser for creating Pegasus workflow applications.

    :param workflow_json_file:
    :type workflow_json_file: str
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 workflow_json_file: str,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow_json_file, logger)

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

    def translate(self, output_file: str) -> None:
        """
        Translates a workflow description (WfFormat) into a Pegasus workflow application.

        :param output_file: The name of the output file (e.g., workflow.py).
        :type output_file: str
        """
        # overall workflow
        self.script += "wf = Workflow('{}', infer_dependencies=True)\n" \
                       "tc = TransformationCatalog()\n" \
                       "rc = ReplicaCatalog()\n\n".format(self.instance.name)
        self.script += "task_output_files = {}\n\n"

        # transformation catalog
        self.script += "full_path = which('sys_test.py')\n" \
                       "if full_path is None:\n" \
                       "    raise RuntimeError('sys_test.py is not in the $PATH')\n" \
                       "base_dir = os.path.dirname(full_path)\n" \
                       "transformation = Transformation('wfperf', site='local',\n" \
                       "                                pfn=os.path.join(base_dir, 'sys_test.py'),\n" \
                       "                                is_stageable=True)\n" \
                       "transformation.add_env(PATH='/usr/bin:/bin:.')\n" \
                       "transformation.add_profiles(Namespace.CONDOR, 'request_disk', '10')\n" \
                       "tc.add_transformations(transformation)\n\n"
        # adding tasks
        for task_name in self.parent_task_names:
            self._add_task(task_name)
            # input file
            input_file = str(uuid.uuid4())
            self.script += "in_file_{} = File('{}')\n".format(self.task_counter, input_file)
            self.script += "rc.add_replica('local', '{}', 'file://' + os.getcwd() + '/data/{}')\n".format(input_file,
                                                                                                          input_file)
            self.script += "{}.add_inputs(in_file_{})\n".format(self.tasks_map[task_name], self.task_counter)
            self.script += "print('{}')\n".format(input_file)

        self.script += "\n"

        # write out the workflow
        self.script += "wf.add_replica_catalog(rc)\n" \
                       "wf.add_transformation_catalog(tc)\n" \
                       "wf.write('{}-benchmark-workflow.yml')\n".format(self.instance.name)

        # write script to file
        with open(output_file, 'w') as out:
            out.write(self.script)

    def _add_task(self, task_name: str, parent_task: Optional[str] = None) -> None:
        """
        Add a task and its dependencies to the workflow.

        :param task_name: name of the task
        :type task_name: str
        :param parent_task: name of the parent task
        :type parent_task: Optional[str]
        """
        if task_name not in self.parsed_tasks:
            task = self.tasks[task_name]
            job_name = "job_{}".format(self.task_counter)
            out_file = str(uuid.uuid4())
            self.script += "{} = Job(transformation=transformation, _id='{}')\n".format(job_name, task_name)
            task.args.insert(0, task_name.split("_")[0])
            task.args.append("--out={}".format(out_file))
            self.script += "{}.add_args({})\n".format(job_name, ", ".join(f"'{a}'" for a in task.args))

            self.script += "out_file_{} = File('{}')\n".format(self.task_counter, out_file)
            self.script += "task_output_files['{}'] = out_file_{}\n".format(job_name, self.task_counter)
            self.script += "{}.add_outputs(out_file_{}, stage_out=False, register_replica=False)\n".format(
                job_name, self.task_counter)

            self.script += "wf.add_jobs({})\n\n".format(job_name)
            self.task_counter += 1
            self.parsed_tasks.append(task_name)
            self.tasks_map[task_name] = job_name

            for node in self.instance.instance['workflow']['jobs']:
                if node['name'] == task_name:
                    for child_task_name in node['children']:
                        self._add_task(child_task_name, job_name)

        if parent_task:
            self.script += "{}.add_inputs(task_output_files['{}'])\n".format(self.tasks_map[task_name], parent_task)
            self.script += "wf.add_dependency({}, parents=[{}])\n\n".format(self.tasks_map[task_name], parent_task)
