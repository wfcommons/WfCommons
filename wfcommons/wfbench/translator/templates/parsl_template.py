#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import parsl

from pathlib import Path
from typing import List

from parsl.app.app import bash_app
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider
from parsl.data_provider.files import File
# NOTE: Uncomment the following lines to enable Parsl workflow monitoring
# from parsl.monitoring.monitoring import MonitoringHub
# from parsl.addresses import address_by_hostname

# NOTE: update the configuration below for your desired system (https://parsl.readthedocs.io/en/stable/userguide/configuring.html)
config = Config(
    executors=[
        HighThroughputExecutor(
            label="local_htex",
            worker_debug=True,
            cores_per_worker=1,
            max_workers_per_node=1,
            provider=LocalProvider(
                init_blocks=1,
                max_blocks=1,
            ),
        )
    ],
    strategy=None,

    # NOTE: Uncomment to enable Parsl workflow monitoring
    # monitoring=MonitoringHub(
    #     hub_address=address_by_hostname(),
    #     hub_port=55055,
    #     monitoring_debug=False,
    #     resource_monitoring_interval=10,
    # ),
)
parsl.clear()
parsl.load(config)

# Emit log lines to the screen
# parsl.set_stream_logger(level=logging.DEBUG)

# Write log to file, specify level of detail for logs
# FILENAME = "debug.log"
# parsl.set_file_logger(FILENAME, level=logging.DEBUG)

@bash_app
def generic_shell_app(cmd: str, inputs=None, outputs=None,
                      stdout="logs/stdout.txt", stderr="logs/stderr.txt"):
    """
    A wrapper for bash apps.
    This wrapper will create the necessary directories for the stdout and stderr files.
    """
    from pathlib import Path

    if inputs is None:
        inputs = []
    if outputs is None:
        outputs = []

    for i in inputs:
        input_path = Path(i.filepath)
        cmd = cmd.replace(input_path.name, i.filepath)
    for o in outputs:
        output_path = Path(o.filepath)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        cmd = cmd.replace(output_path.name, o.filepath)
    return cmd

current_workdir = Path.cwd()

file_map = {}

def get_parsl_files (filenames: List[str], is_output: bool = False) -> List[File]:
    """
    Get a list of Parsl File objects from a list of filenames.
    """
    parsl_files = []

    for filename in filenames:
        if filename not in file_map:
            file_folder = "data"
            if is_output:
                file_folder = "output"
            file_map[filename] = File(current_workdir.joinpath(f"{file_folder}/{filename}"))
        parsl_files.append(file_map[filename])

    return parsl_files

# Generated code goes here
