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


class Translator(ABC):
    """An abstract class of logs parser for creating workflow instances.
    """

    def __init__(self,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        self.logger = logging.getLogger(__name__) if logger is None else logger

    @abstractmethod
    def translate(self, output_file: str) -> None:
        """
        :param output_file: The name of the output file.
        :type output_file: str
        """
