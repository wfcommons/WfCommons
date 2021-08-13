#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import getpass
import json
import networkx as nx

from datetime import datetime
from typing import Optional
from ..common.task import Task
from ..version import __version__

from ..wfchef.utils import create_graph
import tempfile


class Workflow(nx.DiGraph):
    """
    Representation of a workflow. The workflow representation is an extension of the
    `NetworkX DiGraph class <https://networkx.github.io/documentation/stable/reference/classes/digraph.html>`_.

    :param name: Workflow name.
    :type name: str
    :param description: Workflow instance description.
    :type description: str
    :param wms_name: WMS name.
    :type wms_name: str
    :param wms_version: WMS version.
    :type wms_version: str
    :param wms_url: URL for the WMS website.
    :type wms_url: str
    :param executed_at: Workflow start timestamp in the ISO 8601 format.
    :type executed_at: str
    :param makespan: Workflow makespan in seconds.
    :type makespan: int
    """

    def __init__(self,
                 name: str,
                 description: Optional[str] = None,
                 wms_name: Optional[str] = None,
                 wms_version: Optional[str] = None,
                 wms_url: Optional[str] = None,
                 executed_at: Optional[str] = None,
                 makespan: Optional[int] = 0.0
                 ) -> None:
        """Create an object of a workflow representation."""
        self.description = description if description else 'Instance generated with WfCommons - https://wfcommons.org'
        self.created_at = str(datetime.utcnow().isoformat())
        self.schema_version = '1.2'
        self.wms_name = 'WfCommons' if not wms_name else wms_name
        self.wms_version = str(__version__) if not wms_version else wms_version
        self.wms_url = 'https://docs.wfcommons.org/en/v{}/'.format(__version__) if not wms_url else wms_url
        self.executed_at = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z") if not executed_at else executed_at
        self.makespan = makespan
        super().__init__(name=name, makespan=self.makespan, executedat=self.executed_at)

    def write_json(self, json_filename: Optional[str] = None) -> None:
        """Write a JSON file of the workflow instance.

        :param json_filename: JSON output file name.
        :type json_filename: str
        """
        workflow_machines = []
        machines_list = []
        workflow_tasks = []

        workflow_json = {
            'name': self.name,
            'description': self.description,
            'createdAt': self.created_at,
            'schemaVersion': self.schema_version,
            'author': {
                'name': str(getpass.getuser()),
                'email': 'support@wfcommons.org'
            },
            'wms': {
                'name': self.wms_name,
                'version': self.wms_version,
                'url': self.wms_url
            },
            'workflow': {
                'executedAt': self.executed_at,
                'makespan': self.makespan,
                'jobs': workflow_tasks,
                'machines': workflow_machines
            }
        }

        # generate tasks parents and children
        tasks_dependencies = {}
        for edge in self.edges:
            for task_name in edge:
                if task_name not in tasks_dependencies:
                    tasks_dependencies[task_name] = {'parents': [], 'children': []}
            tasks_dependencies[edge[0]]['children'].append(edge[1])
            tasks_dependencies[edge[1]]['parents'].append(edge[0])

        # add tasks to the workflow json object
        for node in self.nodes:
            task: Task = self.nodes[node]['task']
            task_obj = task.as_dict()

            # manage task dependencies
            task_obj['parents'] = tasks_dependencies[task.name]['parents']
            task_obj['children'] = tasks_dependencies[task.name]['children']

            workflow_tasks.append(task_obj)

            # add machines to the workflow json object
            if task.machine and task.machine.name not in machines_list:
                machines_list.append(task.machine.name)
                workflow_machines.append(task.machine.as_dict())

        # write to file
        if not json_filename:
            json_filename = '{}.json'.format(self.name.lower())
        with open(json_filename, 'w') as outfile:
            outfile.write(json.dumps(workflow_json, indent=4))

    def write_dot(self, dot_filename: str = None) -> None:
        """Write a dot file of the workflow instance.

        :param dot_filename: DOT output file name.
        :type dot_filename: str
        """
        if not dot_filename:
            dot_filename = "{}.dot".format(self.name.lower())
        nx.nx_agraph.write_dot(self, dot_filename)

    def to_nx_digraph(self) -> nx.DiGraph:
        with tempfile.NamedTemporaryFile() as temp:
            self.write_json(temp.name)
            return create_graph(temp.name)
