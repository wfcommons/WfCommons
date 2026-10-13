#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import datetime
import dateutil.parser
import importlib.util
import logging
import pathlib

import networkx as nx
import matplotlib.pyplot as plt

from logging import Logger
from typing import Any, Dict, List, Optional

from .schema import SchemaValidator
from ..common.file import File
from ..common.machine import Machine, MachineSystem
from ..common.task import Task, TaskType
from ..common.workflow import Workflow
from ..utils import read_json


class Instance:
    """
    Representation of one execution of one workflow on a set of machines

    .. code-block:: python

        Instance(input_instance = 'instance.json')

    :param input_instance: The JSON instance.
    :type input_instance: pathlib.Path
    :param schema_file: The path to the JSON schema that defines the instance.
                        If no schema file is provided, it will use the schema of the
                        `WfFormat release <https://github.com/wfcommons/wfformat>`_
                        that matches the schema version written by WfCommons.
    :type schema_file: Optional[str]
    :param logger: The logger where to log information/warning or errors.
    :type logger: Optional[Logger]
    """

    def __init__(self, input_instance: pathlib.Path,
                 schema_file: Optional[str] = None,
                 logger: Optional[Logger] = None) -> None:
        """Create an object that represents a workflow execution instance."""
        self.logger: Logger = logging.getLogger(
            __name__) if logger is None else logger

        # Internal variables to be able to iterate directly on an instance
        self._n = 0
        self._order = None

        self.instance: Dict[str, Any] = read_json(input_instance)
        self.logger.info(f"Read a JSON instance: {input_instance}")

        # validate instance
        schema_validator = SchemaValidator(schema_file, logger=logger)
        schema_validator.validate_instance(self.instance)

        # Basic global properties (only name, schemaVersion, and workflow are required)
        self.name: str = self.instance["name"]
        self.desc: Optional[str] = self.instance.get("description", None)
        self.created_at: Optional[datetime.datetime] = dateutil.parser.parse(
            self.instance["createdAt"]) if "createdAt" in self.instance else None
        self.schema_version: str = self.instance["schemaVersion"]

        # Runtime system properties
        self.runtime_system: Dict[str, str] = dict(self.instance.get("runtimeSystem", {}))

        # Author properties
        self.author: Dict[str, str] = dict(self.instance.get("author", {}))

        # Workflow properties
        specification: Dict[str, Any] = self.instance["workflow"]["specification"]
        execution: Optional[Dict[str, Any]] = self.instance["workflow"].get("execution", None)

        # Global properties
        self.executed_at: Optional[datetime.datetime] = dateutil.parser.parse(
            execution["executedAt"]) if execution else None
        self.makespan: float = execution["makespanInSeconds"] if execution else 0.0

        # Machines
        self.machines: Dict[str, Machine] = {
            machine['nodeName']: Machine(
                name=machine['nodeName'],
                cpu=dict(machine.get('cpu', {})),
                system=MachineSystem(machine['system']) if machine.get('system', None) else None,
                architecture=machine.get('architecture', None),
                memory=machine.get('memoryInBytes', None),
                release=machine.get('release', None),
                hashcode=machine.get('machine_code', None),
                logger=self.logger
            ) for machine in (execution.get("machines", []) if execution else [])
        }

        # Files
        files_map = {file["id"]: file["sizeInBytes"] for file in specification.get("files", [])}
        if "files" not in specification and any(
                "inputFiles" in task or "outputFiles" in task for task in specification["tasks"]):
            self.logger.warning("The instance declares no list of files: file sizes are set to 0.")

        # Execution information, by task ID
        execution_tasks: Dict[str, Dict[str, Any]] = {
            t["id"]: t for t in (execution["tasks"] if execution else [])
        }

        # Workflow
        self.workflow: Workflow = Workflow(
            name=self.name,
            description=self.desc,
            created_at=self.instance.get("createdAt", None),
            executed_at=execution["executedAt"] if execution else None,
            makespan=self.makespan,
            runtime_system_name=self.runtime_system.get("name", None),
            runtime_system_url=self.runtime_system.get("url", None),
            runtime_system_version=self.runtime_system.get("version", None),
            author_name=self.author.get("name", None),
            author_email=self.author.get("email", None),
            author_institution=self.author.get("institution", None),
            author_country=self.author.get("country", None),
            include_execution=execution is not None,
            logger=self.logger
        )

        # do not make up the optional properties the instance does not declare
        if self.desc is None:
            self.workflow.description = None
        if self.created_at is None:
            self.workflow.created_at = None
        if not self.author:
            self.workflow.author_name = None
        if not self.runtime_system:
            self.workflow.runtime_system_name = None
        elif "url" not in self.runtime_system:
            self.workflow.runtime_system_url = None
        self.workflow.machines = self.machines

        # Tasks
        for task in specification["tasks"]:
            # Create the list of files associated to this task
            input_files = [File(
                file_id=f,
                size=files_map.get(f, 0),
                logger=self.logger
            ) for f in task.get('inputFiles', [])]
            output_files = [File(
                file_id=f,
                size=files_map.get(f, 0),
                logger=self.logger
            ) for f in task.get('outputFiles', [])]

            t = execution_tasks.get(task["id"], {})
            command = t.get("command", None)

            self.workflow.add_task(Task(
                name=task['name'],
                task_id=task['id'],
                runtime=t.get('runtimeInSeconds', 0),
                category=t.get('category', task.get('category', None)),
                input_files=input_files,
                output_files=output_files,
                cores=t.get('coreCount', None),
                avg_cpu=t.get('avgCPU', None),
                bytes_read=t.get('readBytes', None),
                bytes_written=t.get('writtenBytes', None),
                memory=t.get('memoryInBytes', None),
                energy=t.get('energyInKWh', None),
                avg_power=t.get('avgPowerInW', None),
                priority=t.get('priority', None),
                machines=[self._get_machine(machine) for machine in t.get("machines", [])],
                program=command.get('program', None) if command else None,
                args=command.get('arguments', None) if command else None,
                launch_dir=t.get('launchDir', None),
                executed_at=t.get('executedAt', None),
                logger=self.logger
            ))

        for task in specification["tasks"]:
            for parent in task['parents']:
                self.workflow.add_dependency(parent, task["id"])

        self.logger.info(
            f'Parsed an instance with {len(self.workflow.nodes)} tasks')

    def _get_machine(self, machine_name: str) -> Machine:
        """
        Get a machine by its node name, creating it if the instance does not declare it.

        :param machine_name: Machine node name.
        :type machine_name: str

        :return: The machine.
        :rtype: Machine
        """
        if machine_name not in self.machines:
            self.machines[machine_name] = Machine(name=machine_name, logger=self.logger)
        return self.machines[machine_name]

    def __iter__(self):
        """Produce an iterator based on a topological sort (e.g., scheduling order)"""
        self._n = 0
        self._order = list(nx.topological_sort(self.workflow))
        return self

    def __next__(self) -> str:
        """
        Return the next task from a topological sort.

        :return: task ID
        :rtype: str
        """
        if self._n < len(self.workflow):
            val = self._order[self._n]
            self._n += 1
            return val

        raise StopIteration

    def roots(self) -> List[str]:
        """
        Get the roots of the workflow (i.e., the tasks without any predecessors).

        :return: List of roots
        :rtype: List[str]
        """
        return [n for n, d in self.workflow.in_degree() if d == 0]

    def leaves(self) -> List[str]:
        """
        Get the leaves of the workflow (i.e., the tasks without any successors).

        :return: List of leaves
        :rtype: List[str]
        """
        return [n for n, d in self.workflow.out_degree() if d == 0]

    def write_dot(self, output_path: Optional[pathlib.Path] = None) -> None:
        """
        Write a dot file of the instance.

        :param output_path: The output ``dot`` file name (optional).
        :type output_path: Optional[pathlib.Path]
        """
        self.workflow.write_dot(output_path)

    # # TODO: improve drawing for large instances
    def draw(self, output_path: Optional[pathlib.Path] = None, extension: Optional[str] = "pdf") -> None:
        """
        Produce an image or a pdf file representing the instance.

        :param output_path: Name of the output file.
        :type output_path: Optional[pathlib.Path]
        :param extension: Type of the file extension (``pdf``, ``png``, or ``svg``).
        :type extension: Optional[str]
        """
        graphviz_found = importlib.util.find_spec('pygraphviz')
        if graphviz_found is None:
            self.logger.error(
                f"\'pygraphviz\' package not found: call to {type(self).__name__}.draw() ignored.")
            return

        pos = nx.nx_agraph.graphviz_layout(self.workflow, prog='dot')
        nx.draw(self.workflow, pos=pos, with_labels=False)
        if not output_path:
            output_path = pathlib.Path(f"{self.name.lower()}.{extension}")

        plt.savefig(output_path)
