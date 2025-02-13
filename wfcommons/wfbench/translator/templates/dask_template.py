#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import argparse
import json
import logging
import os
import pathlib
import random
import sys
import time

from dask.distributed import Client


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def build_dask_client():
    """
    Feel free to modify this to target your local dask configuration

    Lots of info there:
    https://docs.dask.org/en/stable/configuration.html
    https://dask.pydata.org/en/latest/scheduling.html
    """
    cpu_count = 2  
    threads_per_cpu = 2
    return Client(n_workers=cpu_count, threads_per_worker=threads_per_cpu)


class WorkflowTask:
    def __init__(self,
                 dag_id: str = None,
                 name: str = None,
                 command_arguments: list[str] = None,
                 inputs: list[str] = None,
                 outputs: list[str] = None,
                 simulate: bool = False,
                 randomizer: random.Random = random.Random(),
                 simulate_minimum_execution_time: float = 0.1,
                 simulate_maximum_execution_time: float = 1.1,
                 execution_time: float = None,  # This is an execution output
                 ):
        self.dag_id = dag_id
        self.name = name
        self.command_arguments = command_arguments
        self.inputs = inputs
        self.outputs = outputs
        self.simulate = simulate
        self.randomizer = randomizer
        self.simulate_minimum_execution_time = simulate_minimum_execution_time
        self.simulate_maximum_execution_time = simulate_maximum_execution_time
        self.execution_time = execution_time

    def simulate_execution(self):
        time.sleep(self.randomizer.uniform(self.simulate_minimum_execution_time,
                                           self.simulate_maximum_execution_time))


def execute_task(task: WorkflowTask, fut_inputs_list) -> WorkflowTask:
    """
    :param task: The task to be executed (it holds all relevant information)
    :param fut_inputs_list: Unused here but necessary for dask to build its own DAG
    :return:
    """
    logger.info("Executing task %s/%s: %s / in=%s / out=%s" % (task.name, task.dag_id, task.command_arguments, task.inputs, task.outputs))
    start = time.time()
    if task.simulate or task.command_arguments is None or len(task.command_arguments) == 0:
        logger.info("Simulating execution of task %s" % task.name)
        # Pretend we do something/Wait some time
        task.simulate_execution()
        for output in task.outputs:
            logger.debug("Simulating %s => %s" % (task.command_arguments, output))
            pathlib.Path(output).touch()
    else:
        command = " ".join(task.command_arguments)
        logger.info("Running command for task %s/%s: %s" % (task.name, task.dag_id, command))
        os.system(command)  # TODO Use subprocess?
    task.execution_time = time.time()-start
    logger.info("End of task %s/%s (%f)" % (task.name, task.dag_id, task.execution_time))
    return task


def run_workflow(client, simulate: bool, seed: int=42) -> list[WorkflowTask]:
# Generated code goes here
    return TASKS


def process_arguments():
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description='Runs a (translated) workflow using Dask')
    parser.add_argument("-sim", "--simulate",
                        help="Simulate all tasks (default: run the tasks for real)", action="store_true")
    parser.add_argument("-s", "--seed", help="Randomizer seed (used only when simulating tasks)", default=42)
    return parser.parse_args()


def to_json(obj):
    return json.dumps(obj, indent=2, default=lambda o: o.__dict__)


if __name__ == '__main__':
    args = process_arguments()
    with build_dask_client() as client:
        tasks = run_workflow(client, args.simulate, seed=int(args.seed))
    with open("run.json", "w") as fp:
        fp.write(to_json(tasks))
