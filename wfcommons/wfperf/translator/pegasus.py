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

        :param output_file: The name of the output file (e.g., workflow.py).
        :type output_file: str
        """
        script = "import os\n" \
                 "from Pegasus.api import *\n\n\n" \
                 "def which(file):\n" \
                 "    for path in os.environ['PATH'].split(os.pathsep):\n" \
                 "        if os.path.exists(os.path.join(path, file)):\n" \
                 "            return os.path.join(path, file)\n" \
                 "    return None\n\n\n"

        # overall workflow
        script += "wf = Workflow('{}', infer_dependencies=True)\n\n".format(self.instance.instance['name'])

        # transformation catalog
        script += "tc = TransformationCatalog()\n" \
                  "full_path = which('sys_test.py')\n" \
                  "if full_path is None:\n" \
                  "    raise RuntimeError('sys_test.py is not in the $PATH')\n" \
                  "base_dir = os.path.dirname(full_path)\n" \
                  "transformation = Transformation('sys_test.py', site='local',\n" \
                  "                                pfn=os.path.join(base_dir, 'sys_test.py'), is_stageable=True)\n" \
                  "transformation.add_env(PATH='/usr/bin:/bin:.')\n" \
                  "tc.add_transformations(transformation)\n" \
                  "wf.add_transformation_catalog(tc)\n\n"

        # replica catalog
        script += "rc = ReplicaCatalog()\n" \
                  "wf.add_replica_catalog(rc)\n\n"

        # write out the workflow
        script += "wf.write('{}-benchmark-workflow.yml')\n".format(self.instance.instance['name'])

        # write script to file
        with open(output_file, 'w') as out:
            out.write(script)
