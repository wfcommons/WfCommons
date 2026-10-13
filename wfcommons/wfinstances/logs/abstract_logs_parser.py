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
from typing import Any, Dict, Hashable, Optional

from ...common.workflow import Workflow
from ...utils import sanitize_wfformat_id


class LogsParser(ABC):
    """An abstract class of logs parser for creating workflow instances.

    :param wms_name: Name of the workflow system.
    :type wms_name: str
    :param wms_version: Version of the workflow system.
    :type wms_version: Optional[str]
    :param wms_url: URL for the workflow system.
    :type wms_url: Optional[str]
    :param description: Workflow instance description.
    :type description: Optional[str]
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Optional[Logger]
    """

    def __init__(self,
                 wms_name: str,
                 wms_version: Optional[str] = None,
                 wms_url: Optional[str] = None,
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the logs parser."""
        self.logger: Optional[Logger] = logger if logger else logging.getLogger(__name__)
        self.description: Optional[str] = description
        self.wms_name: str = wms_name
        self.wms_version: Optional[str] = wms_version
        self.wms_url: Optional[str] = wms_url
        self.workflow = None
        self.workflow_name = None
        self.schema_version = None
        self.executed_at = None
        self.makespan = None
        # WfFormat file IDs, by file key (and the reverse mapping, to detect collisions)
        self._file_ids: Dict[Hashable, str] = {}
        self._file_keys: Dict[str, Hashable] = {}

    def _make_file_id(self, file_name: Any, key: Optional[Hashable] = None) -> str:
        """
        Build a WfFormat-compliant file ID from a file name found in the logs. The same file
        (identified by its key) always gets the same ID, while distinct files whose sanitized
        names collide get distinct IDs (suffixed with "_2", "_3", etc.).

        :param file_name: The file name found in the logs.
        :type file_name: Any
        :param key: A key that uniquely identifies the file in the logs (defaults to the file name).
        :type key: Optional[Hashable]

        :return: A file ID that matches the WfFormat ID pattern.
        :rtype: str
        """
        key = str(file_name) if key is None else key
        if key not in self._file_ids:
            base_id = sanitize_wfformat_id(file_name)
            file_id = base_id
            suffix = 1
            while file_id in self._file_keys:
                suffix += 1
                file_id = f"{base_id}_{suffix}"
            if file_id != str(file_name):
                self.logger.debug(f"File '{file_name}' is identified as '{file_id}'.")
            self._file_ids[key] = file_id
            self._file_keys[file_id] = key
        return self._file_ids[key]

    @abstractmethod
    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow instance based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: Optional[str]

        :return: A workflow instance object.
        :rtype: Workflow
        """
