#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging

from abc import ABC, abstractmethod
from logging import Logger
from typing import Optional

from ...common.workflow import Workflow


class LogsParser(ABC):
    """An abstract class of logs parser for creating workflow traces.

    :param description: Workflow trace description.
    :type description: str
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self, description: Optional[str] = None, logger: Optional[Logger] = None) -> None:
        """Create an object of the logs parser."""
        self.logger = logging.getLogger(__name__) if logger is None else logger
        self.description = description

    @abstractmethod
    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow trace based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: str

        :return: A workflow trace object.
        :rtype: Workflow
        """
