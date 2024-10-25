#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import getpass
import importlib
import json
import networkx as nx
import pathlib

from datetime import datetime
from typing import Optional, List
from ..common.task import Task, TaskType
from ..version import __version__, __schema_version__

from ..wfchef.utils import create_graph
import tempfile


class Workflow(nx.DiGraph):
    """
    Representation of a workflow. The workflow representation is an extension of the
    `NetworkX DiGraph class <https://networkx.github.io/documentation/stable/reference/classes/digraph.html>`_.

    :param name: Workflow name.
    :type name: str
    :param description: Workflow instance description.
    :type description: Optional[str]
    :param runtime_system_name: WMS name.
    :type runtime_system_name: Optional[str]
    :param runtime_system_version: WMS version.
    :type runtime_system_version: Optional[str]
    :param runtime_system_url: URL for the WMS website.
    :type runtime_system_url: Optional[str]
    :param executed_at: Workflow start timestamp in the ISO 8601 format.
    :type executed_at: Optional[str]
    :param makespan: Workflow makespan in seconds.
    :type makespan: Optional[int]
    :param author_name: Author name.
    :type author_name: Optional[str]
    :param author_email: Author email.
    :type author_email: Optional[str]
    :param author_institution: Author institution.
    :type author_institution: Optional[str]
    :param author_country: Author country (preferably country code, ISO ALPHA-2 code).
    :type author_country: Optional[str]
    """

    def __init__(self,
                 name: Optional[str] = "workflow",
                 description: Optional[str] = None,
                 runtime_system_name: Optional[str] = None,
                 runtime_system_version: Optional[str] = None,
                 runtime_system_url: Optional[str] = None,
                 executed_at: Optional[str] = None,
                 makespan: Optional[int] = 0.0,
                 author_name: Optional[str] = None,
                 author_email: Optional[str] = None,
                 author_institution: Optional[str] = None,
                 author_country: Optional[str] = None
                 ) -> None:
        """Create an object of a workflow representation."""
        self.description: Optional[
            str] = description if description else "Instance generated with WfCommons - https://wfcommons.org"
        self.created_at: str = str(datetime.now().astimezone().isoformat())
        self.schema_version: str = f"{__schema_version__}"
        self.runtime_system_name: Optional[str] = "WfCommons" if not runtime_system_name else runtime_system_name
        self.runtime_system_version: Optional[str] = str(__version__) if not runtime_system_version else runtime_system_version
        self.runtime_system_url: Optional[str] = f"https://docs.wfcommons.org/en/v{__version__}/" if not runtime_system_url else runtime_system_url
        self.executed_at: Optional[str] = str(datetime.now().astimezone().isoformat()) if not executed_at else executed_at
        self.makespan: Optional[int] = makespan
        self.author_name: Optional[str] = author_name if author_name else str(getpass.getuser())
        self.author_email: Optional[str] = author_email if author_email else "support@wfcommons.org"
        self.author_institution: Optional[str] = None
        self.author_country: Optional[str] = None
        self.tasks: Task = {}
        self.tasks_parents = {}
        self.tasks_children = {}
        super().__init__(name=name, makespan=self.makespan, executedat=self.executed_at)

    def add_task(self, task: Task) -> None:
        """
        Add a Task to the workflow.

        :param task: A Task object.
        :type task: Task
        """
        self.tasks[task.task_id] = task
        self.tasks_parents.setdefault(task.task_id, set())
        self.tasks_children.setdefault(task.task_id, set())
        self.add_node(task.task_id, task=task, label=task.task_id)

    def add_dependency(self, parent: str, child: str) -> None:
        """
        Add a dependency between tasks.

        :param parent: Parent task id.
        :type parent: str
        :param child: Child task id.
        :type child: str
        """
        self.tasks_parents[child].add(parent)
        self.tasks_children[parent].add(child)
        self.add_edge(parent, child, weight=0)

    def write_json(self, json_file_path: Optional[pathlib.Path] = None) -> None:
        """
        Write a JSON file of the workflow instance.

        :param json_file_path: JSON output file name.
        :type json_file_path: Optional[pathlib.Path]
        """
        workflow_machines = []
        machines_list = []
        specification_tasks = []
        execution_tasks = []
        files = set()

        workflow_json = {
            "name": self.name,
            "description": self.description,
            "createdAt": self.created_at,
            "schemaVersion": self.schema_version,
            "author": {
                "name": self.author_name,
                "email": self.author_email
            },
            "workflow": {
                "specification": {
                    "tasks": specification_tasks,
                    "files": []
                },
                "execution": {
                    "makespanInSeconds": self.makespan,
                    "executedAt": self.executed_at,
                    "tasks": execution_tasks
                }
            },
            "runtimeSystem": {
                "name": self.runtime_system_name,
                "version": self.runtime_system_version,
                "url": self.runtime_system_url
            }
        }

        # generate tasks parents and children
        tasks_dependencies = {}
        for edge in self.edges:
            for task_id in edge:
                if task_id not in tasks_dependencies:
                    tasks_dependencies[task_id] = {"parents": [], "children": []}
            tasks_dependencies[edge[0]]["children"].append(edge[1])
            tasks_dependencies[edge[1]]["parents"].append(edge[0])

        # add tasks to the workflow json object
        for node in self.nodes:
            task: Task = self.nodes[node]["task"]
            task_spec = task.specification_as_dict()
            execution_tasks.append(task.execution_as_dict())

            # manage task dependencies
            if task.task_id in tasks_dependencies:
                task_spec["parents"] = tasks_dependencies[task.task_id]["parents"]
                task_spec["children"] = tasks_dependencies[task.task_id]["children"]

            specification_tasks.append(task_spec)

            # add machines to the workflow json object
            if task.machines:
                for machine in task.machines:
                    if machine.name not in machines_list:
                        machines_list.append(machine.name)
                        workflow_machines.append(machine.as_dict())

            # add files to the workflow json object (input and output)
            for file in task.input_files:
                files.add(file)
            for file in task.output_files:
                files.add(file)

        if workflow_machines:
            workflow_json["workflow"]["execution"]["machines"] = workflow_machines

        if files and len(files) > 0:
            workflow_json["workflow"]["specification"]["files"] = list(file.as_dict() for file in files)

        # write to file
        if not json_file_path:
            json_file_path = pathlib.Path(f"{self.name.lower()}.json")
        with open(json_file_path, "w") as outfile:
            outfile.write(json.dumps(workflow_json, indent=4))
        
        self.workflow_json = workflow_json

    def write_dot(self, dot_file_path: Optional[pathlib.Path] = None) -> None:
        """
        Write a dot file of the workflow instance.

        :param dot_file_path: DOT output file name.
        :type dot_file_path: Optional[pathlib.Path]
        """
        if not dot_file_path:
            dot_file_path = pathlib.Path(f"{self.name.lower()}.dot")
        nx.nx_agraph.write_dot(self, dot_file_path)

    def read_dot(self, dot_file_path: Optional[pathlib.Path] = None) -> None:
        """
        Read a dot file of the workflow instance.

        :param dot_file_path: DOT input file name.
        :type dot_file_path: Optional[pathlib.Path]
        """
        if not dot_file_path:
            raise FileNotFoundError(f"Not able to find the dot file: {dot_file_path}.")
        
        graphviz_found = importlib.util.find_spec('pydot')
        if graphviz_found is None:
            raise ModuleNotFoundError(
                f"\'pydot\' package not found: call to {type(self).__name__}.read_dot() failed.")
        
        graph = nx.drawing.nx_pydot.read_dot(dot_file_path)

        tasks_map = {}
        for node in graph.nodes(data=True):
            task_id = f"{node[1]['label']}_ID{node[0]}"
            self.add_task(Task(name=node[1]['label'], task_id=task_id, runtime=0))
            tasks_map[node[0]] = task_id

        for edge in graph.edges:
            self.add_dependency(tasks_map[edge[0]], tasks_map[edge[1]])

    def to_nx_digraph(self) -> nx.DiGraph:
        with tempfile.NamedTemporaryFile() as temp:
            self.write_json(pathlib.Path(temp.name))
            return create_graph(pathlib.Path(temp.name))

    def roots(self) -> List[Task]:
        return [n for n,d in self.in_degree() if d==0]

    def leaves(self) -> List[Task]:
        return [n for n,d in self.out_degree() if d==0]
