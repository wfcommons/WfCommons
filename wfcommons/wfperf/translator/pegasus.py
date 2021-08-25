#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from logging import Logger
from typing import Optional

from .abstract_translator import Translator


class PegasusTranslator(Translator):
    """A WfFormat parser for creating Pegasus workflow applications.

    :param workflow_json_file:
    :type workflow_json_file: str
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 workflow_json_file: str,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the translator."""
        super().__init__(workflow_json_file, logger)

    def translate(self, output_file: str) -> None:
        """
        Translates a workflow description (WfFormat) into a Pegasus workflow application.

        :param output_file: The name of the output file.
        :type output_file: str
        """
