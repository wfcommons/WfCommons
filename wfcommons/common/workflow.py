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
import logging
import networkx as nx
import pathlib

from datetime import datetime
from logging import Logger
from typing import Any, Dict, Optional, List, Tuple
from ..common.machine import Machine
from ..common.task import Task, TaskType
from ..utils import to_rfc3339
from ..version import __version__, __schema_version__

from ..wfchef.utils import create_graph
import tempfile


def _normalize_timestamp(timestamp: Any, logger: Logger) -> Any:
    """
    Convert a timestamp into the RFC 3339 date-time format required by WfFormat,
    leaving it untouched if it cannot be parsed.
    """
    try:
        return to_rfc3339(timestamp)
    except (ValueError, OverflowError, TypeError):
        logger.warning(f"Unable to convert timestamp '{timestamp}' into the RFC 3339 format.")
        return timestamp


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
    :type makespan: Optional[float]
    :param author_name: Author name.
    :type author_name: Optional[str]
    :param author_email: Author email.
    :type author_email: Optional[str]
    :param author_institution: Author institution.
    :type author_institution: Optional[str]
    :param author_country: Author country (preferably country code, ISO ALPHA-2 code).
    :type author_country: Optional[str]
    :param created_at: Instance creation timestamp in the ISO 8601 format (defaults to now).
    :type created_at: Optional[str]
    :param include_execution: Whether the JSON instance includes the execution section
                              (a WfFormat instance with only a specification is valid).
    :type include_execution: Optional[bool]
    :param logger: The logger where to log information/warning or errors.
    :type logger: Optional[Logger]
    """

    def __init__(self,
                 name: Optional[str] = "workflow",
                 description: Optional[str] = None,
                 runtime_system_name: Optional[str] = None,
                 runtime_system_version: Optional[str] = None,
                 runtime_system_url: Optional[str] = None,
                 executed_at: Optional[str] = None,
                 makespan: Optional[float] = 0.0,
                 author_name: Optional[str] = None,
                 author_email: Optional[str] = None,
                 author_institution: Optional[str] = None,
                 author_country: Optional[str] = None,
                 created_at: Optional[str] = None,
                 include_execution: Optional[bool] = True,
                 logger: Optional[Logger] = None
                 ) -> None:
        """Create an object of a workflow representation."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        name = name if name else "workflow"
        self.description: Optional[
            str] = description if description else "Instance generated with WfCommons - https://wfcommons.org"
        self.created_at: str = created_at if created_at else str(datetime.now().astimezone().isoformat())
        self.schema_version: str = f"{__schema_version__}"
        self.runtime_system_name: Optional[str] = "WfCommons" if not runtime_system_name else runtime_system_name
        self.runtime_system_version: Optional[str] = "unknown" if not runtime_system_version else runtime_system_version
        self.runtime_system_url: Optional[str] = f"https://docs.wfcommons.org/en/v{__version__}/" if not runtime_system_url else runtime_system_url
        self.executed_at: Optional[str] = str(datetime.now().astimezone().isoformat()) if not executed_at else executed_at
        self.makespan: Optional[float] = makespan if makespan is not None else 0.0
        self.author_name: Optional[str] = author_name if author_name else str(getpass.getuser())
        self.author_email: Optional[str] = author_email if author_email else "support@wfcommons.org"
        self.author_institution: Optional[str] = author_institution
        self.author_country: Optional[str] = author_country
        self.include_execution: bool = include_execution
        # machines declared for the workflow, including the ones not used by any task
        self.machines: Dict[str, Machine] = {}
        self.tasks: dict[str, Task] = {}
        self.tasks_parents = {}
        self.tasks_children = {}
        self.workflow_id: str|None = None
        self.workflow_json = {}
        super().__init__(name=name, makespan=self.makespan, executed_at=self.executed_at)

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

        self.generate_json()

         # write to file
        if not json_file_path:
            json_file_path = pathlib.Path(f"{self.name.lower()}.json")
        with open(json_file_path, "w") as outfile:
            outfile.write(json.dumps(self.workflow_json, indent=4))

    def generate_json(self) -> None:
        """
        Generate a JSON representation of the workflow instance.

        :return: A JSON representation of the workflow instance.
        :rtype: dict
        """
        workflow_machines = [machine.as_dict() for machine in self.machines.values()]
        machines_list = list(self.machines.keys())
        specification_tasks = []
        execution_tasks = []
        files = {}

        # optional instance properties are omitted when not set
        workflow_json = {"name": self.name}
        if self.description:
            workflow_json["description"] = self.description
        if self.created_at:
            workflow_json["createdAt"] = _normalize_timestamp(self.created_at, self.logger)
        workflow_json["schemaVersion"] = self.schema_version
        if self.author_name:
            workflow_json["author"] = {
                "name": self.author_name,
                "email": self.author_email
            }
            if self.author_institution:
                workflow_json["author"]["institution"] = self.author_institution
            if self.author_country:
                workflow_json["author"]["country"] = self.author_country
        workflow_json["workflow"] = {
            "specification": {
                "tasks": specification_tasks,
                "files": []
            }
        }
        if self.runtime_system_name:
            workflow_json["runtimeSystem"] = {
                "name": self.runtime_system_name,
                "version": self.runtime_system_version
            }
            if self.runtime_system_url:
                workflow_json["runtimeSystem"]["url"] = self.runtime_system_url
        if self.include_execution:
            workflow_json["workflow"]["execution"] = {
                "makespanInSeconds": self.makespan if self.makespan is not None else 0.0,
                "executedAt": _normalize_timestamp(self.executed_at, self.logger),
                "tasks": execution_tasks
            }

        if len(self.nodes) == 0:
            self.logger.warning("The workflow has no tasks: WfFormat requires at least one task, "
                                "so the generated instance will not validate.")

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

            # add files to the workflow json object (input and output); a file ID must
            # be declared only once, even if tasks report different sizes for it
            for file in task.input_files + task.output_files:
                if file.file_id not in files:
                    files[file.file_id] = file
                elif files[file.file_id].size != file.size:
                    self.logger.warning(f"File '{file.file_id}' has different sizes ({files[file.file_id].size} "
                                        f"and {file.size}): keeping the first one.")

        if workflow_machines and self.include_execution:
            workflow_json["workflow"]["execution"]["machines"] = workflow_machines

        if files:
            workflow_json["workflow"]["specification"]["files"] = [file.as_dict() for file in files.values()]

        self.workflow_json = workflow_json

        # Augment the workflow_json with the metrics
        specification_metrics, execution_metrics = self.compute_metrics(workflow_json)
        if specification_metrics != {}:
            self.workflow_json["workflow"]["specification"]["metrics"] = specification_metrics
        if execution_metrics != {}:
            self.workflow_json["workflow"]["execution"]["metrics"] = execution_metrics

    def compute_metrics(self, workflow_json: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Compute workflow specification and execution metrics, as defined by WfFormat.

        :param workflow_json: Workflow instance as a JSON object.
        :type workflow_json: Dict[str, Any]

        :return: specification and execution metrics.
        :rtype: Tuple[Dict[str, Any], Dict[str, Any]]
        """
        # Specification metrics
        specification = workflow_json["workflow"]["specification"]
        tasks = specification["tasks"]
        specification_metrics = {"numberOfTasks": len(tasks)}
        if "files" in specification:
            files = specification["files"]
            specification_metrics["numberOfFiles"] = len(files)
            if len(files) > 0:
                specification_metrics["sumOfFileSizesInBytes"] = sum(file.get("sizeInBytes", 0) for file in files)

        # A task with no parents is at level 0, and every other task is at one plus the
        # maximum level of any of its parents. The graph is keyed by task ID, since task
        # names are arbitrary and may be shared by many tasks.
        graph = nx.DiGraph()
        graph.add_nodes_from(task["id"] for task in tasks)
        for task in tasks:
            graph.add_edges_from((parent, task["id"]) for parent in task.get("parents", [])
                                 if parent in graph and parent != task["id"])
        try:
            widths = [len(level) for level in nx.topological_generations(graph)]
        except nx.NetworkXUnfeasible:
            self.logger.warning("The workflow graph contains a cycle: DAG level metrics are not computed.")
            widths = []

        if len(widths) > 0:
            specification_metrics["numberOfLevels"] = len(widths)
            specification_metrics["minimumWidth"] = min(widths)
            specification_metrics["maximumWidth"] = max(widths)

        # Execution metrics
        execution_metrics = {}
        if "execution" in workflow_json["workflow"]:
            execution_tasks = workflow_json["workflow"]["execution"]["tasks"]

            total_runtime_in_seconds = sum(task.get("runtimeInSeconds", 0) for task in execution_tasks)
            if total_runtime_in_seconds > 0:
                execution_metrics["sumTaskRuntimesInSeconds"] = total_runtime_in_seconds

            # byte totals are integers in WfFormat, while per-task byte counts may be fractional
            for task_key, metric_key in (("readBytes", "totalNumBytesRead"),
                                         ("writtenBytes", "totalNumBytesWritten")):
                values = [task[task_key] for task in execution_tasks if task_key in task]
                if values:
                    total = sum(values)
                    if float(total).is_integer() and total >= 0:
                        execution_metrics[metric_key] = int(total)
                    else:
                        self.logger.debug(f"Metric '{metric_key}' not computed: {total} is not a "
                                          f"non-negative integer.")

        return specification_metrics, execution_metrics

    def write_dot(self, dot_file_path: Optional[pathlib.Path] = None) -> None:
        """
        Write a dot file of the workflow instance.

        :param dot_file_path: DOT output file name.
        :type dot_file_path: Optional[pathlib.Path]
        """
        if not dot_file_path:
            dot_file_path = pathlib.Path(f"{self.name.lower()}.dot")
        nx.nx_agraph.write_dot(self, str(dot_file_path))

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
        
        # graph = nx.drawing.nx_pydot.read_dot(str(dot_file_path))
        graph = nx.nx_agraph.read_dot(str(dot_file_path))

        # clear everything
        self.tasks.clear()
        self.tasks_parents.clear()
        self.tasks_children.clear()
        self.clear()

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

    def roots(self) -> List[str]:
        return [n for n,d in self.in_degree() if d==0]

    def leaves(self) -> List[str]:
        return [n for n,d in self.out_degree() if d==0]





