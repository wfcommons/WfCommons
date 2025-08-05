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
import io
import tarfile
import os
import sys
import json

from wfcommons import BlastRecipe
from wfcommons.wfbench import WorkflowBenchmark
from wfcommons.wfbench import DaskTranslator
from wfcommons.wfbench import ParslTranslator
from wfcommons.wfbench import NextflowTranslator
from wfcommons.wfbench import AirflowTranslator


def start_docker_container(backend, working_dir):
    # Pulling the Docker image
    client = docker.from_env()
    sys.stderr.write("Pulling Docker image...\n")
    image_name = f"wfcommons/wfcommons-testing-{backend}"

    try:
        image = client.images.get(image_name)
        sys.stderr.write(f"Image '{image_name}' is available locally\n")
    except ImageNotFound:
        sys.stderr.write(f"Pulling image '{image_name}'\n")
        client.images.pull(image_name)

    # Launch the docker container to actually run the translated workflow
    sys.stderr.write("Starting Docker container...\n")
    container = client.containers.run(
        image_name,
        "sleep infinity",
        volumes={working_dir: {'bind': working_dir, 'mode': 'rw'}},
        working_dir=working_dir,
        tty=True,
        detach=True
    )
    return container

def make_tarfile_of_wfcommons():
    source_dir = os.getcwd() # This assumes the testing is run from the root
    tar_stream = io.BytesIO()
    with tarfile.open(fileobj=tar_stream, mode='w') as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
    tar_stream.seek(0)
    return tar_stream

def install_WfCommons_on_container(container):
    sys.stderr.write("Installing WfCommons on the container...\n")
    # Copy the WfCommons code to it (removing stuff that should be removed)
    target_path = '/tmp/'  # inside container
    tar_data = make_tarfile_of_wfcommons()
    container.put_archive(target_path, tar_data)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/build/", stdout=True, stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/*.egg-info/", stdout=True, stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/bin/cpu-benchmark.o", stdout=True,
                                           stderr=True)
    exit_code, output = container.exec_run("sudo /bin/rm -rf /tmp/WfCommons/bin/cpu-benchmark", stdout=True,
                                           stderr=True)

    # Install WfCommons on the container (to install wfbench and cpu-benchmark really)
    exit_code, output = container.exec_run("sudo python3 -m pip install -e . --break-system-packages",
                                           workdir="/tmp/WfCommons", stdout=True, stderr=True)

def create_workflow_benchmark():
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

class TestTranslators:

    @pytest.mark.unit
    # @pytest.mark.skip(reason="tmp")
    def test_dask_translator(self) -> None:

        # Create workflow benchmark
        benchmark, num_tasks = create_workflow_benchmark()

        # Create a local translation directory
        str_dirpath = "/tmp/dask_translated_workflow/"
        dirpath = pathlib.Path(str_dirpath)
        if dirpath.exists():
            shutil.rmtree(dirpath)

        # Perform the translation
        sys.stderr.write("Translating workflow...\n")
        translator = DaskTranslator(benchmark.workflow)
        translator.translate(output_folder=dirpath)

        # Pulling the Docker image
        container = start_docker_container("dask", str_dirpath)

        # Installing WfCommons on container
        install_WfCommons_on_container(container)

        # Copy over the wfbench and cpu-benchmark executables to where they should go
        exit_code, output = container.exec_run("sudo cp -f /tmp/WfCommons/bin/wfbench " + str_dirpath + "bin/", stdout=True, stderr=True)
        exit_code, output = container.exec_run("sudo cp -f /tmp/WfCommons/bin/cpu-benchmark " + str_dirpath + "bin/", stdout=True, stderr=True)

        # Run the workflow!
        sys.stderr.write("Running the Dask workflow on the container...\n")
        exit_code, output = container.exec_run("python ./dask_workflow.py", stdout=True, stderr=True)
        num_completed_tasks = output.decode().count("completed!")  # TODO: This is pretty lame

        # Kill the container
        container.remove(force=True)

        # Do sanity checks
        sys.stderr.write("Checking sanity...\n")
        assert(exit_code == 0)
        assert(num_completed_tasks == num_tasks)
        # TODO: Look at the (I think) generated run.json file on the container


    @pytest.mark.unit
    # @pytest.mark.skip(reason="tmp")
    def test_parsl_translator(self) -> None:

        # Create workflow benchmark
        benchmark, num_tasks = create_workflow_benchmark()

        # Create a local translation directory
        str_dirpath = "/tmp/parsl_translated_workflow/"
        dirpath = pathlib.Path(str_dirpath)
        if dirpath.exists():
            shutil.rmtree(dirpath)

        # Perform the translation
        sys.stderr.write("Translating workflow...\n")
        translator = ParslTranslator(benchmark.workflow)
        translator.translate(output_folder=dirpath)

        # Pulling the Docker image
        container = start_docker_container("parsl", str_dirpath)

        # Installing WfCommons on container
        install_WfCommons_on_container(container)

        # Copy over the wfbench and cpu-benchmark executables to where they should go
        exit_code, output = container.exec_run("sudo cp -f /tmp/WfCommons/bin/wfbench " + str_dirpath + "bin/", stdout=True, stderr=True)
        exit_code, output = container.exec_run("sudo cp -f /tmp/WfCommons/bin/cpu-benchmark " + str_dirpath + "bin/", stdout=True, stderr=True)

        # Run the workflow!
        sys.stderr.write("Running the Parsl workflow on the container...\n")
        exit_code, output = container.exec_run("python ./parsl_workflow.py", stdout=True, stderr=True)

        exit_code, output = container.exec_run(f"cat {str_dirpath}/runinfo/000/parsl.log", stdout=True, stderr=True)
        num_completed_tasks = output.decode().count("_complete_task")

        # Kill the container
        container.remove(force=True)

        # Do sanity checks
        sys.stderr.write("Checking sanity...\n")
        assert(exit_code == 0)
        assert("completed" in output.decode())
        assert(num_completed_tasks == num_tasks)

    @pytest.mark.unit
    # @pytest.mark.skip(reason="tmp")
    def test_nextflow_translator(self) -> None:

        # Create workflow benchmark
        benchmark, num_tasks = create_workflow_benchmark()

        # Create a local translation directory
        str_dirpath = "/tmp/nextflow_translated_workflow/"
        dirpath = pathlib.Path(str_dirpath)
        if dirpath.exists():
            shutil.rmtree(dirpath)

        # Perform the translation
        sys.stderr.write("Translating workflow...\n")
        translator = NextflowTranslator(benchmark.workflow)
        translator.translate(output_folder=dirpath)

        # Pulling the Docker image
        container = start_docker_container("nextflow", str_dirpath)

        # Installing WfCommons on container
        install_WfCommons_on_container(container)

        # Copy over the wfbench and cpu-benchmark executables to where they should go
        exit_code, output = container.exec_run("sudo cp -f /tmp/WfCommons/bin/wfbench " + str_dirpath + "bin/",
                                               stdout=True, stderr=True)
        exit_code, output = container.exec_run("sudo cp -f /tmp/WfCommons/bin/cpu-benchmark " + str_dirpath + "bin/",
                                               stdout=True, stderr=True)

        # Run the workflow!
        sys.stderr.write("Running the Nextflow workflow on the container...\n")
        exit_code, output = container.exec_run(f"nextflow run ./workflow.nf --pwd .", stdout=True, stderr=True)
        ignored, task_exit_codes = container.exec_run("find . -name .exitcode -exec cat {} \;", stdout=True, stderr=True)

        # Kill the container
        container.remove(force=True)

        # Do sanity checks
        sys.stderr.write("Checking sanity...\n")
        assert (exit_code == 0)
        assert (task_exit_codes.decode() == num_tasks * "0")


    @pytest.mark.unit
    @pytest.mark.skip(reason="tmp")
    def test_airflow_translator(self) -> None:

        # Create workflow benchmark
        benchmark, num_tasks = create_workflow_benchmark()

        # Create a local translation directory
        str_dirpath = "/tmp/airflow_translated_workflow/"
        dirpath = pathlib.Path(str_dirpath)
        if dirpath.exists():
            shutil.rmtree(dirpath)

        # Perform the translation
        sys.stderr.write("Translating workflow...\n")
        translator = AirflowTranslator(benchmark.workflow)
        translator.translate(output_folder=dirpath)

        # Pulling the Docker image
        container = start_docker_container("airflow", str_dirpath)

        # Installing WfCommons on container
        install_WfCommons_on_container(container)

        # Copy over the wfbench and cpu-benchmark executables to where they should go
        # exit_code, output = container.exec_run("sudo cp -f /tmp/WfCommons/bin/wfbench " + str_dirpath + "bin/",
        #                                        stdout=True, stderr=True)
        # exit_code, output = container.exec_run("sudo cp -f /tmp/WfCommons/bin/cpu-benchmark " + str_dirpath + "bin/",
        #                                        stdout=True, stderr=True)
        #
        # # Run the workflow!
        # sys.stderr.write("Running the Airflow workflow on the container...\n")
        # exit_code, output = container.exec_run(f"nextflow run ./workflow.nf --pwd .", stdout=True, stderr=True)
        # ignored, task_exit_codes = container.exec_run("find . -name .exitcode -exec cat {} \;", stdout=True, stderr=True)
        #
        # # Kill the container
        # container.remove(force=True)
        #
        # # Do sanity checks
        # sys.stderr.write("Checking sanity...\n")
        # assert (exit_code == 0)
        # assert (task_exit_codes.decode() == num_tasks * "0")
