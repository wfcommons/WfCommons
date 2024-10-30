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

from parsl.app.app import bash_app
from parsl.config import Config
from parsl.data_provider.files import File

# temp for debugging
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider
from parsl.channels import LocalChannel


local_htex = Config(
    executors=[
        HighThroughputExecutor(
            label="htex_Local",
            worker_debug=True,
            cores_per_worker=1,
            provider=LocalProvider(
                channel=LocalChannel(),
                init_blocks=1,
                max_blocks=1,
            ),
        )
    ],
    strategy=None,
)

parsl.clear()
parsl.load(local_htex)

# end temp

@bash_app
def benchmark(task_id, percent_cpu, cpu_work, inputs=[], outputs=[]):
    output_file = f"\"{{\\\"{outputs[0]}\\\":169492}}\""
    return f'./wfbench {task_id} --percent-cpu {percent_cpu} --cpu-work {cpu_work} --out {output_file} {" ".join(inputs)}'


# Generated code goes here