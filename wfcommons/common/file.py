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

from logging import Logger
from typing import Dict, Optional

from ..utils import NoValue


class FileLink(NoValue):
    """Type of file link."""
    INPUT = "input"
    OUTPUT = "output"


class File:
    """Representation of a file.

    :param name: The name of the file.
    :type name: str
    :param size: File size in KB.
    :type size: int
    :param link: Type of file link.
    :type link: FileLink
    :param logger: The logger where to log information/warning or errors.
    :type logger: Logger
    """

    def __init__(self, name: str, size: int, link: FileLink, logger: Optional[Logger] = None) -> None:
        """A file used by tasks."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger

        self.name: str = name
        self.size: int = size
        self.link: FileLink = link

    def as_dict(self) -> Dict:
        """A JSON representation of the file.

        :return: A JSON object representation of the file.
        :rtype: Dict
        """
        return {
            'link': self.link.value,
            'name': self.name,
            'size': self.size
        }
