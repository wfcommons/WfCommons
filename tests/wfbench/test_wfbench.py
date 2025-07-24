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
import shutil
import json

from wfcommons import BlastRecipe
from wfcommons.common import Workflow
from wfcommons.wfbench import WorkflowBenchmark


class TestWfBench:

    @pytest.mark.unit
    def test_create_from_recipe(self) -> None:
        """
        Very minimal testing here for creating from recipe
        """
        # create a workflow benchmark object to generate specifications based on a recipe
        benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)

        # generate a specification based on performance characteristics
        dirpath = pathlib.Path("/tmp/benchmark/")
        if dirpath.exists():
            shutil.rmtree(dirpath)
        dirpath.mkdir(parents=True, exist_ok=True)
        path = benchmark.create_benchmark(dirpath, cpu_work=100, data=10, percent_cpu=0.6)


    @staticmethod
    def _directory_content_as_expected(dirpath: pathlib.Path,
                                       workflow: Workflow,
                                       num_tasks: int,
                                       cpu_work: int,
                                       percent_cpu: float):
        workflow_file_name = f"{workflow.name.lower()}-{num_tasks}.json"
        workflow_file_is_there = (dirpath / workflow_file_name).exists()
        to_create_file_is_there = (dirpath / "to_create.txt").exists()
        return workflow_file_is_there and to_create_file_is_there


    @staticmethod
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


    @staticmethod
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


    @pytest.mark.unit
    def test_create_from_instance(self) -> None:
        workflow = BlastRecipe.from_num_tasks(500).build_workflow()
        benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
        dirpath = pathlib.Path("/tmp/benchmark/")
        if dirpath.exists():
            shutil.rmtree(dirpath)
        dirpath.mkdir(parents=True, exist_ok=True)
        path = benchmark.create_benchmark_from_synthetic_workflow(dirpath, workflow, cpu_work=100, percent_cpu=0.6)
        assert(self._directory_content_as_expected(dirpath, workflow, 500, 100, 0.6))
        assert(self._workflow_as_expected(dirpath, workflow, 500, 100, 0.6))
        assert(self._to_create_file_as_expected(dirpath, workflow, 500, 100, 0.6))

