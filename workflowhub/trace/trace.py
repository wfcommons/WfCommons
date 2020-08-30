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
import dateutil.parser
import importlib.util
import logging

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


class Trace:
    """Representation of one execution of one workflow on a set of machines

    .. code-block:: python

        Trace(input_trace = 'trace.json')

    :param input_trace: The JSON trace.
    :type input_trace: str
    :param schema_file: The path to the JSON schema that defines the trace.
                        If no schema file is provided, it will look for a local
                        copy of the WorkflowHub schema, and if not available
                        it will fetch the latest schema from the
                        `WorkflowHub schema GitHub <https://github.com/workflowhub/workflow-schema>`_
                        repository.
    :type schema_file: str
    :param logger: The logger where to log information/warning or errors.
    :type logger: Logger
    """

    def __init__(self, input_trace: str, schema_file: Optional[str] = None, logger: Optional[Logger] = None) -> None:
        """Create an object that represents a workflow execution trace."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger

        # Internal variables to be able to iterate directly on a trace
        self._n = 0
        self._order = None

        self.trace: Dict[str, Any] = read_json(input_trace)
        self.logger.info("Read a JSON trace: " + input_trace)

        # validate trace
        schema_validator = SchemaValidator(schema_file)
        schema_validator.validate_trace(self.trace)

        # Basic global properties
        self.name: str = self.trace['name']
        self.desc: str = self.trace['description']
        self.created_at: datetime = dateutil.parser.parse(self.trace['createdAt'])
        self.schema_version: str = self.trace['schemaVersion']

        # WMS properties
        self.wms: Dict[str, str] = {
            u: v for u, v in self.trace['wms'].items()
        }

        # Author properties
        self.author: Dict[str, str] = {
            u: v for u, v in self.trace['author'].items()
        }

        # Workflow properties
        # Global properties
        self.executed_at: datetime = dateutil.parser.parse(self.trace['workflow']['executedAt'])
        self.makespan: int = self.trace['workflow']['makespan']

        # Machines
        self.machines: Dict[str, Machine] = {
            machine['nodeName']: Machine(
                name=machine['nodeName'],
                cpu={k: v for k, v in machine['cpu'].items()},
                system=MachineSystem(machine.get('system', None)),
                architecture=machine.get('architecture', None),
                memory=machine.get('memory', None),
                release=machine.get('release', None),
                hashcode=machine.get('machine_code', None),
                logger=self.logger
            ) for machine in self.trace['workflow']['machines']
        }

        # Tasks
        self.workflow: Workflow = Workflow(name=self.name, makespan=self.makespan)
        for task in self.trace['workflow']['jobs']:
            # Required arguments are defined in the JSON scheme
            # Here name, type and runtime are required
            # By default the value is set to None if we do not find the value

            # Create the list of files associated to this task
            list_files = task.get('files', [])
            list_files = [File(
                name=f['name'],
                size=f['size'],
                link=FileLink(f['link']),
                logger=self.logger
            ) for f in list_files]

            # Fetch back the machine associated to this task
            machine = task.get('machine', None)
            machine = None if machine is None else self.machines[machine]

            self.workflow.add_node(
                task['name'],
                task=Task(
                    name=task['name'],
                    task_type=TaskType(task['type']),
                    runtime=task['runtime'],
                    machine=machine,
                    args=task.get('arguments', None),
                    cores=task.get('cores', None),
                    avg_cpu=task.get('avgCPU', None),
                    bytes_read=task.get('bytesRead', None),
                    bytes_written=task.get('bytesWritten', None),
                    memory=task.get('memory', None),
                    energy=task.get('energy', None),
                    avg_power=task.get('avgPower', None),
                    priority=task.get('priority', None),
                    files=list_files,  # TODO: sum all files read/written by this task
                    logger=self.logger
                )
            )

        # TODO: handle the case of the output files of the leaves tasks (not taken into account yet)
        for task in self.trace['workflow']['jobs']:
            for parent in task['parents']:
                self.workflow.add_edge(parent, task['name'], weight=0)

        # TODO: instead of attaching files to tasks, attach them to edges based on the link direction.
        self.logger.info('Parsed a trace with {} tasks'.format(len(self.workflow.nodes)))

    def __iter__(self):
        """Produce an iterator based on a topological sort (e.g., scheduling order)"""
        self._n = 0
        self._order = list(nx.topological_sort(self.workflow))
        return self

    def __next__(self) -> str:
        """Return the next task from a topological sort.

        :return: task ID
        :rtype: str
        """
        if self._n < len(self.workflow):
            val = self._order[self._n]
            self._n += 1
            return val
        else:
            raise StopIteration

    def roots(self) -> List[str]:
        """Get the roots of the workflow (i.e., the tasks without any predecessors).

        :return: List of roots
        :rtype: List[str]
        """
        return [n for n, d in self.workflow.in_degree() if d == 0]

    def leaves(self) -> List[str]:
        """Get the leaves of the workflow (i.e., the tasks without any successors).

        :return: List of leaves
        :rtype: List[str]
        """
        return [n for n, d in self.workflow.out_degree() if d == 0]

    def write_dot(self, output: Optional[str] = None) -> None:
        """Write a dot file of the trace.

        :param output: The output ``dot`` file name (optional).
        :type output: str
        """
        self.workflow.write_dot(output)

    # # TODO: improve drawing for large traces
    def draw(self, output: Optional[str] = None, extension: str = "pdf") -> None:
        """Produce an image or a pdf file representing the trace.

        :param output: Name of the output file.
        :type output: str
        :param extension: Type of the file extension (``pdf``, ``png``, or ``svg``).
        :type output: str
        """
        graphviz_found = importlib.util.find_spec('pygraphviz')
        if graphviz_found is None:
            self.logger.error(
                "\'pygraphviz\' package not found: call to {0}.draw() ignored.".format(type(self).__name__))
            return

        pos = nx.nx_pydot.graphviz_layout(self.workflow, prog='dot')
        nx.draw(self.workflow, pos=pos, with_labels=False)
        if not output:
            output = "{0}.{1}".format(self.name.lower(), extension)

        plt.savefig(output)
