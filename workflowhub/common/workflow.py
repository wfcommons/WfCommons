#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import datetime
import getpass
import json
import networkx as nx

from typing import Optional
from ..common.task import Task
from ..version import __version__

class Workflow(nx.DiGraph):
    """
    Representation of a workflow. The workflow representation is an extension of the
    `NetworkX DiGraph class <https://networkx.github.io/documentation/stable/reference/classes/digraph.html>`_.

    :param name: Workflow name.
    :type name: str
    :param makespan: Workflow makespan in seconds.
    :type makespan: int
    """

    def __init__(self, name: str, makespan: Optional[int]) -> None:
        """Create an object of a workflow representation."""
        self.makespan = makespan
        self.executedat = datetime.datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
        super().__init__(name=name, makespan=self.makespan, executedat=self.executedat)

    def write_json(self, json_filename: Optional[str] = None) -> None:
        """Write a JSON file of the workflow trace.

        :param json_filename: JSON output file name.
        :type json_filename: str
        """
        workflow_json = {
            'name': self.name,
            'description': 'Trace generated with WorkflowHub - https://workflowhub.org',
            'createdAt': str(datetime.datetime.utcnow().isoformat()),
            'schemaVersion': '1.0',
            'author': {
                'name': str(getpass.getuser()),
                'email': 'support@workflowhub.org'
            },
            'wms': {
                'name': 'WorkflowHub',
                'version': str(__version__),
                'url': 'https://workflowhub.readthedocs.io/en/v{}/'.format(__version__)
            },
            'workflow': {
                'executedAt': self.executedat,
                'makespan': 0.0 if not self.makespan else self.makespan,
                'machines': [
                    {
                        'nodeName': 'fake-1',
                        'system': 'linux', 
                        'cpu': {
                            'count': 1, 
                            'speed': 1
                        }
                    }
                ],
                'jobs': []
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
            task_files = []
            for f in task.files:
                task_files.append({'link': f.link.value,
                                  'name': f.name,
                                  'size': f.size})

            workflow_json['workflow']['jobs'].append({'name': task.name,
                                                      'type': task.type.value,
                                                      'runtime': task.runtime,
                                                      'parents': tasks_dependencies[task.name]['parents'],
                                                      'children': tasks_dependencies[task.name]['children'],
                                                      'files': task_files
                                                      })

        # write to file
        if not json_filename:
            json_filename = '{}.json'.format(self.name.lower())
        with open(json_filename, 'w') as outfile:
            outfile.write(json.dumps(workflow_json, indent=4))

    def write_dot(self, dot_filename: str = None) -> None:
        """Write a dot file of the workflow trace.

        :param dot_filename: DOT output file name.
        :type dot_filename: str
        """
        if not dot_filename:
            dot_filename = "{}.dot".format(self.name.lower())
        nx.nx_agraph.write_dot(self, dot_filename)
