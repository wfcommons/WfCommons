#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pytest

import pathlib
import sys
import json

from tests.test_helpers import _create_fresh_local_dir
from tests.test_helpers import _start_docker_container
from tests.test_helpers import _remove_local_dir_if_it_exists
from tests.test_helpers import _get_total_size_of_directory
from wfcommons import BlastRecipe
from wfcommons.common import Workflow
from wfcommons.wfbench import WorkflowBenchmark, BashTranslator

def _directory_content_as_expected(dirpath: pathlib.Path,
                                   workflow: Workflow,
                                   num_tasks: int,
                                   cpu_work: int,
                                   percent_cpu: float):
    workflow_file_name = f"{workflow.name.lower()}-{num_tasks}.json"
    workflow_file_is_there = (dirpath / workflow_file_name).exists()
    to_create_file_is_there = (dirpath / "to_create.txt").exists()
    return workflow_file_is_there and to_create_file_is_there


def _workflow_as_expected(dirpath: pathlib.Path,
                          workflow: Workflow,
                          num_tasks: int,
                          cpu_work: int,
                          percent_cpu: float):
    # Get the generated JSON
    json_path = dirpath / f"{workflow.name.lower()}-{num_tasks}.json"
    with json_path.open("r") as f:
        generated_json = json.load(f)

    # Check the number of tasks
    assert(len(workflow.tasks) == len(generated_json['workflow']['specification']['tasks']))

    # For each task check sanity
    for generated_task in generated_json['workflow']['specification']['tasks']:
        assert(generated_task['id'] in workflow.tasks)
        workflow_task = workflow.tasks[generated_task['id']]
        # Check input files
        assert(len(generated_task['inputFiles']) == len(workflow_task.input_files))
        for file in workflow_task.input_files:
            assert(file.file_id in generated_task['inputFiles'])
        # Check output files
        assert(len(generated_task['outputFiles']) == len(workflow_task.output_files))
        for file in workflow_task.input_files:
            assert(file.file_id in generated_task['inputFiles'])

    # TODO: Implement more sanity checks

    return True


def _to_create_file_as_expected(dirpath: pathlib.Path,
                          workflow: Workflow,
                          num_tasks: int,
                          cpu_work: int,
                          percent_cpu: float):

    # Build a set of output files
    outputfile_set = set({})
    for task_name in workflow.tasks:
        for f in workflow.tasks[task_name].output_files:
            outputfile_set.add(f.file_id)

    # Build a dict of input files that must be created files with sizes
    inputfile_dict = dict({})
    for task_name in workflow.tasks:
        for f in workflow.tasks[task_name].input_files:
            if f.file_id not in outputfile_set:
                inputfile_dict[f.file_id] = f.size

    # Open the "to_create.txt" file and go line by line and check names/sizes
    json_path = dirpath / "to_create.txt"
    with json_path.open("r") as f:
        for line in f.readlines():
            [filename, size] = line.strip().split(" ")
            assert(filename in inputfile_dict)
            assert(inputfile_dict[filename] == int(size))

    return True

def _actual_data_files_as_expected(dirpath: pathlib.Path,
                                   workflow: Workflow,
                                   data_spec):
    # Inspect the data footprint generated
    if isinstance(data_spec, int):
        total_bytes = _get_total_size_of_directory(str(dirpath / "data"))
        sys.stderr.write(f"Total observed data footprint is {total_bytes} bytes\n")
        assert (abs(data_spec * 1000 * 1000 - total_bytes) < 100)
    else:
        sys.stderr.write(f"Unsupported data spec FOR NOW\n")

    # Inspect the input/output files existence
    for task in workflow.tasks.values():
        for f in task.input_files:
            filename = f.file_id.split("/")[-1]
            assert (pathlib.Path.exists(dirpath / "data" / filename))
        for f in task.output_files:
            filename = f.file_id.split("/")[-1]
            assert (pathlib.Path.exists(dirpath / "data" / filename))

class TestWfBench:

    @pytest.mark.unit
    def test_create_from_recipe(self) -> None:
        """
        Very minimal testing here for creating from recipe
        """
        # create a workflow benchmark object to generate specifications based on a recipe
        benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=49)
        workflow: Workflow = benchmark.workflow
        sys.stderr.write(f"WORKFLOW = {workflow}\n")

        # Create the data_specification options
        fixed_total_footprint_in_mb = 5
        per_type_footprint = {}
        for task_type in ["blastall", "split_fasta", None]:
                per_type_footprint[task_type] = "1"

        for data_spec in [per_type_footprint]:
            benchmark.create_benchmark(_create_fresh_local_dir(f"/tmp/benchmark_{len(data_spec)}"), cpu_work=1, data=data_spec, percent_cpu=0.6)

            # Run the benchmark with the Bash translator
            # Create a local translation directory
            str_dirpath = "/tmp/bash_translated_benchmark/"
            _remove_local_dir_if_it_exists(str_dirpath)
            dirpath = pathlib.Path(str_dirpath)

            # Perform the translation
            sys.stderr.write("\nTranslating workflow...\n")
            translator = BashTranslator(benchmark.workflow)
            translator.translate(output_folder=dirpath)

            # Start the Docker container
            sys.stderr.write("Starting Docker container...\n")
            container = _start_docker_container("bash", str_dirpath, str_dirpath, str_dirpath + "bin/")

            # Run the workflow
            sys.stderr.write("Running workflow...\n")
            exit_code, output = container.exec_run(cmd="/bin/bash ./run_workflow.sh", stdout=True, stderr=True)

            # Kill the container
            container.remove(force=True)

            # Inspect the data after execution
            _actual_data_files_as_expected(dirpath, workflow, data_spec)


    @pytest.mark.unit
    def test_create_from_instance(self) -> None:
        workflow = BlastRecipe.from_num_tasks(500).build_workflow()
        benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
        dirpath = _create_fresh_local_dir("/tmp/benchmark/")
        path = benchmark.create_benchmark_from_synthetic_workflow(dirpath, workflow, cpu_work=100, percent_cpu=0.6)
        assert(_directory_content_as_expected(dirpath, workflow, 500, 100, 0.6))
        assert(_workflow_as_expected(dirpath, workflow, 500, 100, 0.6))
        assert(_to_create_file_as_expected(dirpath, workflow, 500, 100, 0.6))

