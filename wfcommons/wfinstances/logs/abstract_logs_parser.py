#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
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
    """An abstract class of logs parser for creating workflow instances.

    :param wms_name: Name of the workflow system.
    :type wms_name: str
    :param wms_url: URL for the workflow system.
    :type wms_url: str
    :param description: Workflow instance description.
    :type description: str
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 wms_name: str,
                 wms_url: Optional[str] = None,
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the logs parser."""
        self.logger = logging.getLogger(__name__) if logger is None else logger
        self.description = description
        self.wms_name = wms_name
        self.wms_url = wms_url
        self.workflow = None
        self.workflow_name = None
        self.schema_version = None
        self.executed_at = None
        self.makespan = None

    @abstractmethod
    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow instance based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: str

        :return: A workflow instance object.
        :rtype: Workflow
        """
