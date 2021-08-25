#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging

from abc import ABC, abstractmethod
from logging import Logger
from typing import Optional

from ...wfinstances.instance import Instance


class Translator(ABC):
    """An abstract class of WfFormat parser for creating workflow applications.

    :param workflow_json_file:
    :type workflow_json_file: str
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 workflow_json_file: str,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        self.logger = logging.getLogger(__name__) if logger is None else logger
        self.instance = Instance(workflow_json_file, logger=logger)

        # find all tasks
        self.tasks = {}
        for node in self.instance.workflow.nodes.data():
            self.tasks[node[0]] = node[1]['task']

        # find parent tasks
        self.parent_task_names = []
        for node in self.instance.instance['workflow']['jobs']:
            if len(node['parents']) == 0 and node['name'] not in self.parent_task_names:
                self.parent_task_names.append(node['name'])

    @abstractmethod
    def translate(self, output_file: str) -> None:
        """
        Translates a workflow description (WfFormat) into an actual workflow application.

        :param output_file: The name of the output file.
        :type output_file: str
        """
