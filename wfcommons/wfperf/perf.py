#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import glob
import json
import logging
import os
import pathlib
import subprocess
import uuid

from logging import Logger
from typing import Dict, Optional, Union, List, Type

from ..wfchef.wfchef_abstract_recipe import WfChefWorkflowRecipe
from ..wfgen import WorkflowGenerator

this_dir = pathlib.Path(__file__).resolve().parent


class WorkflowBenchmark:
    """

    :param recipe:
    :type recipe: Type[WfChefWorkflowRecipe]
    :param num_tasks:
    :type num_tasks: int
    :param logger: The logger where to log information/warning or errors.
    :type logger: Optional[Logger]
    """

    def __init__(self, recipe: Type[WfChefWorkflowRecipe], num_tasks: int, logger: Optional[Logger] = None) -> None:
        """Create an object that represents a workflow benchmark."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        self.recipe = recipe
        self.num_tasks = num_tasks

    def create(self,
               save_dir: pathlib.Path,
               percent_cpu: float = 0.6,
               data_footprint: int = 1,
               test_mode: str = "seqwr",
               mem_total_size: str = "1000000G",
               block_size: str = "1K",
               scope: str = "global",
               max_prime: int = 10000,
               file_block_size: int = 16384,
               rw_ratio: float = 1.5,
               max_time: int = 150,
               create: bool = True,
               path: Optional[pathlib.Path] = None) -> pathlib.Path:

        save_dir = save_dir.resolve()
        save_dir.mkdir(exist_ok=True, parents=True)

        if create:
            self.logger.debug("Generating workflow")
            generator = WorkflowGenerator(self.recipe.from_num_tasks(self.num_tasks))
            workflow = generator.build_workflow()
            name = f"{workflow.name.split('-')[0]}-Benchmark"
            workflow_savepath = save_dir.joinpath(f"{name}-{self.num_tasks}").with_suffix(".json")
            workflow.write_json(str(workflow_savepath))
            wf = json.loads(workflow_savepath.read_text())
        else:
            wf = json.loads(path.read_text())

        # Creating the lock files
        self.logger.debug("Creating lock files.")
        with save_dir.joinpath("cores.txt").open("w+") as fp, save_dir.joinpath("cores.txt.lock").open("w+") as fpl:
            pass

        lock = save_dir.joinpath("cores.txt.lock")
        cores = save_dir.joinpath("cores.txt")

        data = "--data" if data_footprint else ""

        # Setting the parameters for the arguments section of the JSON
        params = [data,
                  f"--file-test-mode={test_mode}",
                  f"--file-total-size={data_footprint}G",
                  f"--file-block-size={file_block_size}",
                  f"--file-rw-ratio={rw_ratio}",
                  "--file-num=1",
                  "--forced-shutdown=0",
                  f"--path-lock={lock}",
                  f"--path-cores={cores}",
                  f"--memory-block-size={block_size}",
                  f"--memory-scope={scope}",
                  f"--memory-total-size={mem_total_size}",
                  f"--cpu-max-prime={max_prime}",
                  f"--percent-cpu={percent_cpu}",
                  f"--time={max_time}"]

        wf["name"] = name
        for job in wf["workflow"]["jobs"]:
            job["files"] = []
            job.setdefault("command", {})
            job["command"]["program"] = f"{this_dir.joinpath('wfperf_benchmark.py')}"
            job["command"]["arguments"] = [job["name"]]
            job["command"]["arguments"].extend(params)
            if "runtime" in job:
                del job["runtime"]

        num_sys_files, num_total_files = input_files(wf)

        self.logger.debug(f"Number of input files to be created by the system: {num_sys_files}")
        self.logger.debug(f"Total number of files used by the workflow: {num_total_files}")
        file_size = round(data_footprint / num_total_files)
        self.logger.debug(f"Every input/output file is of size: {file_size}")
        add_io_to_json(wf, file_size)

        self.logger.debug("Generating system files.")
        generate_sys_data(num_sys_files, file_size)

        self.logger.debug("Removing system files.")
        cleanup_sys_files()

        json_path = save_dir.joinpath(f"{name}-{self.num_tasks}").with_suffix(".json")
        self.logger.info(f"Saving benchmark workflow: {json_path}")
        json_path.write_text(json.dumps(wf, indent=4))

        return json_path

    def run(self, json_path: pathlib.Path, save_dir: pathlib.Path):
        """Run the benchmark workflow locally (for test purposes only).
        """
        self.logger.debug("Running")
        try:
            wf = json.loads(json_path.read_text())
            with save_dir.joinpath(f"run.txt").open("w+") as fp:
                procs: List[subprocess.Popen] = []
                for item in wf["workflow"]["jobs"]:
                    executable = item["command"]["program"]
                    arguments = item["command"]["arguments"]
                    program = [
                        "time", "python", executable,
                        item["name"].split("_")[0],
                        *arguments,
                    ]
                    folder = pathlib.Path(this_dir.joinpath(f"wfperf_execution/{uuid.uuid4()}"))
                    folder.mkdir(exist_ok=True, parents=True)
                    os.chdir(str(folder))
                    procs.append(subprocess.Popen(program, stdout=fp, stderr=fp))
                    os.chdir("../..")

                for proc in procs:
                    proc.wait()

        except Exception as e:
            import traceback
            traceback.print_exc()
            raise FileNotFoundError("Not able to find the executable.")


def generate_sys_data(num_files: int, file_total_size: int, test_mode: str = "seqwr"):
    """Generate workflow's input data
    """
    params = [
        f"--file-num={num_files}",
        f"--file-total-size={file_total_size}G",
        f"--file-test-mode={test_mode}"
    ]
    proc = subprocess.Popen(["sysbench", "fileio", *params, "prepare"], stdout=subprocess.PIPE)
    proc.wait()
    out, _ = proc.communicate()
    if not out:
        raise FileNotFoundError("Couldn't create files. Check parameters.")

    # have to change the name back to test_file on bash
    for t in glob.glob("test_file.*"):
        os.rename(t, f"sys_{t}")


def cleanup_sys_files():
    """Remove files already used
    """
    for t in glob.glob("*test_file.*"):
        new = t.split("_")[1:]
        new = "_".join(new)
        os.rename(t, new)
    proc = subprocess.Popen(["sysbench", "fileio", "cleanup"], stdout=subprocess.PIPE)
    proc.wait()
    out, _ = proc.communicate()
    if not out:
        raise FileNotFoundError("Couldn't delete files.")


def add_io_to_json(wf: Dict[str, Dict], file_size: int) -> None:
    """Add input and output files to JSON
    """
    i = 0
    all_jobs = {
        job["name"]: job
        for job in wf["workflow"]["jobs"]
    }

    for job in wf["workflow"]["jobs"]:
        job.setdefault("files", [])
        job["files"].append(
            {
                "link": "output",
                "name": f"{job['name']}_test_file.0",
                "size": file_size
            }
        )

        parents = [parent for parent in job["parents"]]
        if not parents:
            job["files"].append(
                {
                    "link": "input",
                    "name": f"sys_test_file.{i}",
                    "size": file_size
                }
            )
            i += 1
        else:
            for parent in parents:
                job["files"].extend(
                    [
                        {
                            "link": "input",
                            "name": item["name"],
                            "size": item["size"]
                        }
                        for item in all_jobs[parent]["files"] if item["link"] == "output"
                    ]
                )


def input_files(wf: Dict[str, Dict]):
    """Calculate total number of files needed
    """
    tasks_need_input = 0
    tasks_dont_need_input = 0

    for job in wf["workflow"]["jobs"]:
        parents = [parent for parent in job["parents"]]
        if not parents:
            tasks_need_input += 1
        else:
            tasks_dont_need_input += 1

    total_num_files = tasks_need_input * 2 + tasks_dont_need_input

    return tasks_need_input, total_num_files
