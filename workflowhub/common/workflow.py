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
from ..common.job import Job


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
        super().__init__(name=name, makespan=makespan)

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
            'workflow': {
                'jobs': []
            }
        }

        # generate jobs parents and children
        jobs_dependencies = {}
        for edge in self.edges:
            for job_name in edge:
                if job_name not in jobs_dependencies:
                    jobs_dependencies[job_name] = {'parents': [], 'children': []}
            jobs_dependencies[edge[0]]['children'].append(edge[1])
            jobs_dependencies[edge[1]]['parents'].append(edge[0])

        # add jobs to the workflow json object
        for node in self.nodes:
            job: Job = self.nodes[node]['job']
            job_files = []
            for f in job.files:
                job_files.append({'link': f.link.value,
                                  'name': f.name,
                                  'size': f.size})

            workflow_json['workflow']['jobs'].append({'name': job.name,
                                                      'type': job.type.value,
                                                      'runtime': job.runtime,
                                                      'parents': jobs_dependencies[job.name]['parents'],
                                                      'children': jobs_dependencies[job.name]['children'],
                                                      'files': job_files
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
