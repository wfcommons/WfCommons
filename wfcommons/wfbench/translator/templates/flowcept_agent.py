#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import ast
import logging
import pathlib
import sys
import time
from flowcept.flowcept_api.flowcept_controller import Flowcept

logging.basicConfig(
    level=logging.INFO,
    format="[WfBench][%%(asctime)s][%%(levelname)s] %%(message)s",
    datefmt="%%H:%%M:%%S",
    handlers=[logging.StreamHandler()]
)

workflow_name = sys.argv[1]
workflow_id = sys.argv[2]
out_files = ast.literal_eval(sys.argv[3])

logging.info("Flowcept Starting")
flowcept_agent = Flowcept(workflow_id=workflow_id, workflow_name=workflow_name, bundle_exec_id=workflow_id, start_persistence=False, save_workflow=True)

try:
    flowcept_agent.start()
except Exception:
    import traceback
    traceback.print_exc()

remaining_files = set(out_files)

while remaining_files:
    found_files = set()
    for f in remaining_files:
        if pathlib.Path(f).exists():
            found_files.add(f)
    remaining_files -= found_files
    if not remaining_files:
        break
    time.sleep(1)
    
try:
    flowcept_agent.stop()
    time.sleep(5)
except Exception:
    import traceback
    traceback.print_exc()

logging.info("Flowcept Completed")