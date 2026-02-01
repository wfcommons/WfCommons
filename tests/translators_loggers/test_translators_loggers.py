#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025-2026 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import pytest
import shutil
import sys
import json
import time
import re
import os

from tests.test_helpers import _create_fresh_local_dir
from tests.test_helpers import _remove_local_dir_if_it_exists
from tests.test_helpers import _start_docker_container
from tests.test_helpers import _shutdown_docker_container_and_remove_image
from tests.test_helpers import _compare_workflows

from wfcommons import BlastRecipe
from wfcommons.common import Workflow, Task
from wfcommons.wfbench import WorkflowBenchmark
from wfcommons.wfbench import DaskTranslator
from wfcommons.wfbench import ParslTranslator
from wfcommons.wfbench import NextflowTranslator
from wfcommons.wfbench import AirflowTranslator
from wfcommons.wfbench import BashTranslator
from wfcommons.wfbench import TaskVineTranslator
from wfcommons.wfbench import MakeflowTranslator
from wfcommons.wfbench import CWLTranslator
from wfcommons.wfbench import PegasusTranslator
from wfcommons.wfbench import SwiftTTranslator

from wfcommons.wfinstances import PegasusLogsParser
from wfcommons.wfinstances.logs import TaskVineLogsParser
from wfcommons.wfinstances.logs import MakeflowLogsParser


def _create_workflow_benchmark() -> (WorkflowBenchmark, int):
    # Create a workflow benchmark object to generate specifications based on a recipe (in /tmp/, whatever)
    desired_num_tasks = 45
    benchmark_full_path = "/tmp/blast-benchmark-{desired_num_tasks}.json"
    shutil.rmtree(benchmark_full_path, ignore_errors=True)
    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=desired_num_tasks)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=10, data=10, percent_cpu=0.6)
    with open(f"/tmp/blast-benchmark-{desired_num_tasks}.json", "r") as f:
        generated_json = json.load(f)
        num_tasks = len(generated_json["workflow"]["specification"]["tasks"])
    return benchmark, num_tasks

#############################################################################
## Methods to perform additional setup required by the backend
#############################################################################

def noop(container):
    pass

def _additional_setup_taskvine(container):
    # Create the poncho package
    exit_code, output = container.exec_run(cmd=["bash", "-c",
                                                "source ~/conda/etc/profile.d/conda.sh && conda activate && poncho_package_create taskvine_poncho.json taskvine_poncho.tar.gz"],
                                           stdout=True, stderr=True)
    if exit_code != 0:
        raise Exception("Failed to setup TaskVine: cannot create poncho package")
    # Start a vine worker in the background
    exit_code, output = container.exec_run(
        cmd=["bash", "-c", "source ~/conda/etc/profile.d/conda.sh && conda activate && vine_worker localhost 9123"],
        detach=True, stdout=True, stderr=True)
   # Note that exit_code will always be None because of detach=True. So hopefully this works.
   # TODO?: check that the vine_worker is running (so as to abort early)

def _additional_setup_pegasus(container):
    # Start Condor
    exit_code, output = container.exec_run(cmd=["bash", "-c",
                                                "bash /home/wfcommons/start_condor.sh"],
                                           stdout=True, stderr=True)
    if exit_code != 0:
        raise Exception("Failed to setup Pegasus: cannot start HTCondor")
    # Run pegasus script
    exit_code, output = container.exec_run(cmd=["bash", "-c",
                                                "python3 ./pegasus_workflow.py"],
                                           stdout=True, stderr=True)
    if exit_code != 0:
        raise Exception("Failed to setup Pegasus: error while running the pegasus_workflow.py script")

def _additional_setup_swiftt(container):
    # Start a redis server in the background
    exit_code, output = container.exec_run(
        cmd=["bash", "-c", "redis-server"], detach=True, stdout=True, stderr=True)
    # Note that exit_code will always be None because of detach=True.

    # Check that the redis-server is up
    exit_code, output = container.exec_run(
        cmd=["bash", "-c", "redis-cli ping"], stdout=True, stderr=True)
    if output.decode().strip() != 'PONG':
        raise Exception("Failed to start redis-server...")

additional_setup_methods = {
    "dask": noop,
    "parsl": noop,
    "nextflow": noop,
    "airflow": noop,
    "bash": noop,
    "taskvine": _additional_setup_taskvine,
    "makeflow": noop,
    "cwl": noop,
    "pegasus": _additional_setup_pegasus,
    "swiftt": _additional_setup_swiftt,
}

#############################################################################
## Methods to run the various workflows, kill the container, and check sanity
#############################################################################

def run_workflow_dask(container, num_tasks, str_dirpath):
    exit_code, output = container.exec_run("python ./dask_workflow.py", stdout=True, stderr=True)
    # Check sanity
    assert (exit_code == 0)
    assert (output.decode().count("completed!")  == num_tasks)
    # TODO: Look at the (I think) generated run.json file on the container?

def run_workflow_parsl(container, num_tasks, str_dirpath):
    exit_code, output = container.exec_run("python ./parsl_workflow.py", stdout=True, stderr=True)
    ignored, output = container.exec_run(f"cat {str_dirpath}/runinfo/000/parsl.log", stdout=True, stderr=True)
    # Check sanity
    assert (exit_code == 0)
    assert ("completed" in output.decode())
    assert (output.decode().count("_complete_task") == num_tasks)

def run_workflow_nextflow(container, num_tasks, str_dirpath):
    # Run the workflow!
    exit_code, output = container.exec_run(f"nextflow run ./workflow.nf --pwd .", stdout=True, stderr=True)
    ignored, task_exit_codes = container.exec_run("find . -name .exitcode -exec cat {} \;", stdout=True, stderr=True)
    # Check sanity
    assert (exit_code == 0)
    assert (task_exit_codes.decode() == num_tasks * "0")

def run_workflow_airflow(container, num_tasks, str_dirpath):
    # Run the workflow! (use a specific working directory)
    # TODO: Remove the hardcoded Blast-Benchmark as it's ugly
    exit_code, output = container.exec_run(cmd=["sh", "-c", "cd /home/wfcommons/ && sudo /bin/bash /run_a_workflow.sh Blast-Benchmark"],
                                           stdout=True,
                                           stderr=True)
    # Check sanity
    assert (exit_code == 0)
    assert (output.decode().count("completed") == num_tasks * 2)

def run_workflow_bash(container, num_tasks, str_dirpath):
    # Run the workflow!
    exit_code, output = container.exec_run(cmd="/bin/bash ./run_workflow.sh", stdout=True, stderr=True)
    # Check sanity
    assert (exit_code == 0)
    assert (output.decode().count("completed") == num_tasks)

def run_workflow_taskvine(container, num_tasks, str_dirpath):
    # Run the workflow!
    exit_code, output = container.exec_run(cmd=["bash", "-c", "source ~/conda/etc/profile.d/conda.sh && conda activate && python3 ./taskvine_workflow.py"], stdout=True, stderr=True)
    # Check sanity
    assert (exit_code == 0)
    assert (output.decode().count("completed") == num_tasks)

def run_workflow_makeflow(container, num_tasks, str_dirpath):
    # Run the workflow (with full logging)
    exit_code, output = container.exec_run(cmd=["bash", "-c", "source ~/conda/etc/profile.d/conda.sh && conda activate && makeflow --log-verbose  --monitor=./monitor_data/ ./workflow.makeflow"], stdout=True, stderr=True)
    # Check sanity
    assert (exit_code == 0)
    num_completed_jobs = len(re.findall(r'job \d+ completed', output.decode()))
    assert (num_completed_jobs == num_tasks)

def run_workflow_cwl(container, num_tasks, str_dirpath):
    # Run the workflow!
    # Note that the input file is hardcoded and Blast-specific
    exit_code, output = container.exec_run(cmd="cwltool ./main.cwl --split_fasta_00000001_input ./data/workflow_infile_0001 ", stdout=True, stderr=True)
    # Check sanity
    assert (exit_code == 0)
    # this below is ugly (the 3 is for "workflow", "compile_output_files" and "compile_log_files",
    # and there is a 2* because there is a message for the job and for the step)
    assert (output.decode().count("completed success") == 3 + 2 *num_tasks)

def run_workflow_pegasus(container, num_tasks, str_dirpath):
    # Run the workflow!
    exit_code, output = container.exec_run(cmd="bash /home/wfcommons/run_workflow.sh", stdout=True, stderr=True)
    # Check sanity
    assert(exit_code == 0)
    assert("success" in output.decode())

def run_workflow_swiftt(container, num_tasks, str_dirpath):
    # Run the workflow!
    exit_code, output = container.exec_run(cmd="swift-t workflow.swift", stdout=True, stderr=True)
    # sys.stderr.write(output.decode())
    # Check sanity
    assert(exit_code == 0)
    assert (output.decode().count("completed!") == num_tasks)
    pass

run_workflow_methods = {
    "dask": run_workflow_dask,
    "parsl": run_workflow_parsl,
    "nextflow": run_workflow_nextflow,
    "airflow": run_workflow_airflow,
    "bash": run_workflow_bash,
    "taskvine": run_workflow_taskvine,
    "makeflow": run_workflow_makeflow,
    "cwl": run_workflow_cwl,
    "pegasus": run_workflow_pegasus,
    "swiftt": run_workflow_swiftt,
}

translator_classes = {
    "dask": DaskTranslator,
    "parsl": ParslTranslator,
    "nextflow": NextflowTranslator,
    "airflow": AirflowTranslator,
    "bash": BashTranslator,
    "taskvine": TaskVineTranslator,
    "makeflow": MakeflowTranslator,
    "cwl": CWLTranslator,
    "pegasus": PegasusTranslator,
    "swiftt": SwiftTTranslator,
}


class TestTranslators:

    @pytest.mark.parametrize(
        "backend",
        [
           "swiftt",
           "dask",
           "parsl",
           "nextflow",
           "airflow",
           "bash",
           "taskvine",
           "makeflow",
           "cwl",
           "pegasus",
        ])
    @pytest.mark.unit
    # @pytest.mark.skip(reason="tmp")
    def test_translator(self, backend) -> None:
        # Create workflow benchmark
        benchmark, num_tasks = _create_workflow_benchmark()

        # Create a local translation directory
        str_dirpath = "/tmp/" + backend + "_translated_workflow/"
        dirpath = pathlib.Path(str_dirpath)
        # dirpath = _create_fresh_local_dir(str_dirpath)
        _remove_local_dir_if_it_exists(str_dirpath)

        # Perform the translation
        sys.stderr.write(f"\n[{backend}] Translating workflow...\n")
        translator = translator_classes[backend](benchmark.workflow)
        translator.translate(output_folder=dirpath)

        # # Make the directory that holds the translation world-writable,
        # # so that docker commands won't fail
        # TODO: Explore whether this below makes tests runnable on Linux due to
        #       different Docker permission schemes, etc.
        # os.chmod(dirpath, 0o777)

        # Start the Docker container
        container = _start_docker_container(backend, str_dirpath, str_dirpath, str_dirpath + "bin/")

        # Do whatever necessary setup
        additional_setup_methods[backend](container)

        # Run the workflow
        sys.stderr.write(f"[{backend}] Running workflow...\n")
        start_time = time.time()
        run_workflow_methods[backend](container, num_tasks, str_dirpath)
        sys.stderr.write(f"[{backend}] Workflow ran in %.2f seconds\n" % (time.time() - start_time))

        # Run the log parser if any
        if backend == "pegasus":
            parser = PegasusLogsParser(dirpath / "work/wfcommons/pegasus/Blast-Benchmark/run0001/")
        elif backend == "taskvine":
            parser = TaskVineLogsParser(dirpath / "vine-run-info/most-recent/vine-logs", filenames_to_ignore=["cpu-benchmark","stress-ng", "wfbench"])
        elif backend == "makeflow":
            parser = MakeflowLogsParser(execution_dir = dirpath, resource_monitor_logs_dir = dirpath / "monitor_data/")
        else:
            parser = None

        if parser:
            sys.stderr.write(f"\n[{backend}] Parsing the logs...\n")
            reconstructed_workflow : Workflow = parser.build_workflow(f"reconstructed_workflow_{backend}")
            reconstructed_workflow.write_json(pathlib.Path("/tmp/reconstructed_workflow.json"))

            original_workflow : Workflow = benchmark.workflow

            _compare_workflows(original_workflow, reconstructed_workflow)

        # Shutdown the container (weirdly, container is already shutdown by now... not sure how)
        _shutdown_docker_container_and_remove_image(container)

