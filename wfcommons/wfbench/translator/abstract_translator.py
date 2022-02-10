#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2022 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging
import pathlib

from abc import ABC, abstractmethod
from typing import Optional, List

from ...common.task import Task
from ...wfinstances.instance import Instance


class Translator(ABC):
    """
    An abstract class of WfFormat parser for creating workflow benchmark applications.

    :param workflow_json_file_path: Path to the workflow benchmark JSON instance.
    :type workflow_json_file_path: pathlib.Path
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: logging.Logger
    """

    def __init__(self,
                 workflow_json_file_path: pathlib.Path,
                 logger: Optional[logging.Logger] = None) -> None:
        """Create an object of the translator."""
        self.logger = logging.getLogger(__name__) if logger is None else logger
        self.workflow_json_file_path = workflow_json_file_path
        self.instance = Instance(workflow_json_file_path, logger=logger)

        # find all tasks
        self.tasks = {}
        for node in self.instance.workflow.nodes.data():
            self.tasks[node[0]] = node[1]["task"]

        # find parent tasks
        self.parent_task_names = []
        for node in self.instance.instance["workflow"]["tasks"]:
            if len(node["parents"]) == 0 and node["name"] not in self.parent_task_names:
                self.parent_task_names.append(node["name"])

    @abstractmethod
    def translate(self, output_file_name: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_file_name: The name of the output file.
        :type output_file_name: pathlib.Path
        """

    def _write_output_file(self, contents: str, output_file_name: pathlib.Path) -> None:
        """
        Write the translated content to a file.

        :param contents: Contents to be written to the output file.
        :type contents: str
        :param output_file_name: The name of the output file.
        :type output_file_name: pathlib.Path
        """
        # file will be written to the same folder as for the original JSON instance.
        out_file = self.workflow_json_file_path.parent.joinpath(
            output_file_name)

        with open(out_file, "w") as out:
            out.write(contents)
        self.logger.info(f"Translated content written to '{out_file}'")

    def _find_children(self, task_name: str) -> List[Task]:
        """
        Find the children for a specific task.

        :param task_name: The task name.
        :type task_name: str

        :return: List of task's children.
        :rtype: List[Task]
        """
        children = None
        for node in self.instance.instance["workflow"]["tasks"]:
            if node["name"] == task_name:
                children = node["children"]

        return children

    def _find_parents(self, task_name: str) -> List[Task]:
        """
        Find the parents for a specific task.

        :param task_name: The task name.
        :type task_name: str

        :return: List of task's parents.
        :rtype: List[Task]
        """
        parents = None
        for node in self.instance.instance["workflow"]["tasks"]:
            if node["name"] == task_name:
                parents = node["parents"]

        return parents
