#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import pytest
import shutil
import docker
from docker.errors import ImageNotFound
import io
import tarfile
import os
import sys
import json
import time

from wfcommons import BlastRecipe
from wfcommons.wfbench import WorkflowBenchmark, TaskVineTranslator
from wfcommons.wfbench import DaskTranslator
from wfcommons.wfbench import ParslTranslator
from wfcommons.wfbench import NextflowTranslator
from wfcommons.wfbench import AirflowTranslator
from wfcommons.wfbench import BashTranslator


def _start_docker_container(backend, mounted_dir, working_dir, bin_dir, command=["sleep", "infinity"]):
    # Pulling the Docker image
    client = docker.from_env()
    image_name = f"wfcommons/wfcommons-testing-{backend}"

    try:
        image = client.images.get(image_name)
        sys.stderr.write(f"Image '{image_name}' is available locally\n")
    except ImageNotFound:
        sys.stderr.write(f"Pulling image '{image_name}'...\n")
        client.images.pull(image_name)

    # Launch the docker container to actually run the translated workflow
    sys.stderr.write("Starting Docker container...\n")
    container = client.containers.run(
        image_name,
        command=command,
        volumes={mounted_dir: {'bind': mounted_dir, 'mode': 'rw'}},
        working_dir=working_dir,
        tty=True,
        detach=True
    )

    # Installing WfCommons on container
    _install_WfCommons_on_container(container)

    # Copy over the wfbench and cpu-benchmark executables to where they should go on the container
    if bin_dir:
        exit_code, output = container.exec_run(["sh", "-c", "sudo cp -f `which wfbench` " + bin_dir],
                                               stdout=True, stderr=True)
        exit_code, output = container.exec_run(["sh", "-c", "sudo cp -f `which cpu-benchmark` " + bin_dir],
                                               stdout=True, stderr=True)

    return container

def _make_tarfile_of_wfcommons():
    source_dir = os.getcwd() # This assumes the testing is run from the root
    tar_stream = io.BytesIO()
    with tarfile.open(fileobj=tar_stream, mode='w') as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
    tar_stream.seek(0)
    return tar_stream

def _install_WfCommons_on_container(container):
    # sys.stderr.write("Installing WfCommons on the container...\n")
    # Copy the WfCommons code to it (removing stuff that should be removed)
    target_path = '/tmp/'  # inside container
    tar_data = _make_tarfile_of_wfcommons()
    container.put_archive(target_path, tar_data)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/build/", stdout=True, stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/*.egg-info/", stdout=True, stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/bin/cpu-benchmark.o", stdout=True,
                                           stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/bin/cpu-benchmark", stdout=True,
                                           stderr=True)

    # Install WfCommons on the container (to install wfbench and cpu-benchmark really)
    exit_code, output = container.exec_run("sudo python3 -m pip install . --break-system-packages",
                                           workdir="/tmp/WfCommons", stdout=True, stderr=True)

def _create_workflow_benchmark():
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
    # Start a vine worker in the background
    exit_code, output = container.exec_run(
        cmd=["bash", "-c", "source ~/conda/etc/profile.d/conda.sh && conda activate && vine_worker localhost 9123"],
        detach=True, stdout=True, stderr=True)


additional_setup_methods = {
    "dask": noop,
    "parsl": noop,
    "nextflow": noop,
    "airflow": noop,
    "bash": noop,
    "taskvine": _additional_setup_taskvine,
}

#############################################################################
## Methods to run the various workflows, kill the container, and check sanity
#############################################################################

def run_workflow_dask(container, num_tasks, str_dirpath):
    exit_code, output = container.exec_run("python ./dask_workflow.py", stdout=True, stderr=True)
    # Kill the container
    container.remove(force=True)
    # Check sanity
    assert (exit_code == 0)
    assert (output.decode().count("completed!")  == num_tasks)
    # TODO: Look at the (I think) generated run.json file on the container?

def run_workflow_parsl(container, num_tasks, str_dirpath):
    exit_code, output = container.exec_run("python ./parsl_workflow.py", stdout=True, stderr=True)
    ignored, output = container.exec_run(f"cat {str_dirpath}/runinfo/000/parsl.log", stdout=True, stderr=True)
    # Kill the container
    container.remove(force=True)
    # Check sanity
    assert (exit_code == 0)
    assert ("completed" in output.decode())
    assert (output.decode().count("_complete_task") == num_tasks)

def run_workflow_nextflow(container, num_tasks, str_dirpath):
    # Run the workflow!
    exit_code, output = container.exec_run(f"nextflow run ./workflow.nf --pwd .", stdout=True, stderr=True)
    ignored, task_exit_codes = container.exec_run("find . -name .exitcode -exec cat {} \;", stdout=True, stderr=True)
    # Kill the container
    container.remove(force=True)
    # Check sanity
    assert (exit_code == 0)
    assert (task_exit_codes.decode() == num_tasks * "0")

def run_workflow_airflow(container, num_tasks, str_dirpath):
    # Run the workflow! (use a specific working directory)
    # TODO: Remove the hardcoded Blast-Benchmark as it's ugly
    exit_code, output = container.exec_run(cmd=["sh", "-c", "cd /home/wfcommons/ && sudo /bin/bash /run_a_workflow.sh Blast-Benchmark"],
                                           stdout=True,
                                           stderr=True)
    # Kill the container
    container.remove(force=True)
    # Check sanity
    assert (exit_code == 0)
    assert (output.decode().count("completed") == num_tasks * 2)

def run_workflow_bash(container, num_tasks, str_dirpath):
    # Run the workflow!
    exit_code, output = container.exec_run(cmd="/bin/bash ./run_workflow.sh", stdout=True, stderr=True)
    # Kill the container
    container.remove(force=True)
    # Check sanity
    assert (exit_code == 0)
    assert (output.decode().count("completed") == num_tasks)

def run_workflow_taskvine(container, num_tasks, str_dirpath):

    # Run the workflow!
    exit_code, output = container.exec_run(cmd=["bash", "-c", "source ~/conda/etc/profile.d/conda.sh && conda activate && python3 ./taskvine_workflow.py"], stdout=True, stderr=True)
    # Kill the container
    container.remove(force=True)
    # # Check sanity
    assert (exit_code == 0)
    assert (output.decode().count("completed") == num_tasks)


run_workflow_methods = {
    "dask": run_workflow_dask,
    "parsl": run_workflow_parsl,
    "nextflow": run_workflow_nextflow,
    "airflow": run_workflow_airflow,
    "bash": run_workflow_bash,
    "taskvine": run_workflow_taskvine,
}

translator_classes = {
    "dask": DaskTranslator,
    "parsl": ParslTranslator,
    "nextflow": NextflowTranslator,
    "airflow": AirflowTranslator,
    "bash": BashTranslator,
    "taskvine": TaskVineTranslator,
}


class TestTranslators:

    @pytest.mark.parametrize(
        "backend",
        [
            "dask",
            "parsl",
            "nextflow",
            "airflow",
            "bash",
            "taskvine",
        ])
    @pytest.mark.unit
    # @pytest.mark.skip(reason="tmp")
    def test_translator(self, backend) -> None:
        # Create workflow benchmark
        benchmark, num_tasks = _create_workflow_benchmark()

        # Create a local translation directory
        str_dirpath = "/tmp/" + backend + "_translated_workflow/"
        dirpath = pathlib.Path(str_dirpath)
        if dirpath.exists():
            shutil.rmtree(dirpath)

        # Perform the translation
        sys.stderr.write("\nTranslating workflow...\n")
        translator = translator_classes[backend](benchmark.workflow)
        translator.translate(output_folder=dirpath)

        # Start the Docker container
        sys.stderr.write("Starting Docker container...\n")
        container = _start_docker_container(backend, str_dirpath, str_dirpath, str_dirpath + "bin/")

        # Do whatever necessary setup
        additional_setup_methods[backend](container)

        # Run the workflow
        sys.stderr.write("Running workflow...\n")
        start_time = time.time()
        run_workflow_methods[backend](container, num_tasks, str_dirpath)
        sys.stderr.write("Workflow ran in %.2f seconds\n" % (time.time() - start_time))