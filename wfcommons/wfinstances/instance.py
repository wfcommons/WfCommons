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
from ..common.file import File, FileLink
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
                        If no schema file is provided, it will look for a local
                        copy of the WfFormat, and if not available it will fetch
                        the latest schema from the
                        `WfFormat schema GitHub <https://github.com/wfcommons/wfformat>`_
                        repository.
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

        # Basic global properties
        self.name: str = self.instance["name"]
        self.desc: str = self.instance["description"]
        self.created_at: datetime = dateutil.parser.parse(
            self.instance["createdAt"])
        self.schema_version: str = self.instance["schemaVersion"]

        # Runtime system properties
        self.runtime_system: Dict[str, str] = {
            u: v for u, v in self.instance["runtimeSystem"].items()
        }

        # Author properties
        self.author: Dict[str, str] = {
            u: v for u, v in self.instance["author"].items()
        }

        # Workflow properties
        # Global properties
        self.executed_at: datetime = dateutil.parser.parse(
            self.instance["workflow"]["execution"]["executedAt"])
        self.makespan: int = self.instance["workflow"]["execution"]["makespanInSeconds"]

        # Machines
        if "machines" in self.instance["workflow"]["execution"].keys():
            self.machines: Dict[str, Machine] = {
                machine['nodeName']: Machine(
                    name=machine['nodeName'],
                    cpu={k: v for k, v in machine['cpu'].items()},
                    system=MachineSystem(machine.get('system', None)) if machine.get(
                        'system', None) else None,
                    architecture=machine.get('architecture', None),
                    memory=machine.get('memoryInBytes', None),
                    release=machine.get('release', None),
                    hashcode=machine.get('machine_code', None),
                    logger=self.logger
                ) for machine in self.instance["workflow"]["execution"]["machines"]
            }

        # Files
        files_map = {}
        for file in self.instance["workflow"]["specification"]["files"]:
            files_map[file["id"]] = file["sizeInBytes"]

        # Tasks
        tasks_map = {}
        for task in self.instance["workflow"]["specification"]["tasks"]:
            # Required arguments are defined in the JSON scheme
            # Here name, type and runtime are required
            # By default the value is set to None if we do not find the value

            # Create the list of files associated to this task
            input_files = [File(
                file_id=f,
                size=files_map[f],
                link=FileLink.INPUT,
                logger=self.logger
            ) for f in task.get('inputFiles', [])]
            output_files = [File(
                file_id=f,
                size=files_map[f],
                link=FileLink.OUTPUT,
                logger=self.logger
            ) for f in task.get('outputFiles', [])]

            tasks_map[task['id']] = Task(
                name=task['name'],
                task_id=task['id'],
                runtime=0,
                category=task.get('category', None),
                input_files=input_files,
                output_files=output_files,
                logger=self.logger
            )

        # Workflow
        self.workflow: Workflow = Workflow(
            name=self.name, 
            makespan=self.makespan,
            runtime_system_name=self.runtime_system["name"],
            runtime_system_url=self.runtime_system["url"],
            runtime_system_version=self.runtime_system["version"],
            author_name=self.author["name"],
            author_email=self.author["email"]
        )

        for t in self.instance["workflow"]["execution"]["tasks"]:
            task = tasks_map[t["id"]]
            task.runtime=t['runtimeInSeconds'] if 'runtimeInSeconds' in t else 0
            task.cores=t.get('coreCount', None)
            task.avg_cpu=t.get('avgCPU', None)
            task.bytes_read=t.get('readBytes', None)
            task.bytes_written=t.get('writtenBytes', None)
            task.memory=t.get('memoryInBytes', None)
            task.energy=t.get('energyInMhz', None)
            task.avg_power=t.get('avgPowerInMhz', None)
            task.priority=t.get('priority', None)
            task.start_time=t.get('executedAt', None)

            # Fetch back the machine associated to this task
            machines_list = t["machines"] if "machines" in t else []
            machines = []
            for machine in machines_list:
                machines.append(self.machines[machine])
            task.machines = machines

            # Fetch the command associated to this task
            command = t.get("command", None)
            task.program=command.get('program', None) if command else None
            task.args=command.get('arguments', None) if command else None

            self.workflow.add_task(task)

        for task in self.instance["workflow"]["specification"]["tasks"]:
            for parent in task['parents']:
                self.workflow.add_dependency(parent, task["id"])

        self.logger.info(
            f'Parsed an instance with {len(self.workflow.nodes)} tasks')

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

        pos = nx.nx_pydot.graphviz_layout(self.workflow, prog='dot')
        nx.draw(self.workflow, pos=pos, with_labels=False)
        if not output_path:
            output_path = pathlib.Path(f"{self.name.lower()}.{extension}")

        plt.savefig(output_path)
