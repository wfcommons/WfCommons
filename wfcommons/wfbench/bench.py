#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import ast
import glob
import json
import logging
import os
import pathlib
import re
import subprocess
import time
import shortuuid
import sys

from logging import Logger
from typing import Dict, Optional, List, Set, Tuple, Type, Union

from ..common import File, Task, Workflow

from ..wfchef.wfchef_abstract_recipe import WfChefWorkflowRecipe
from ..wfgen import WorkflowGenerator

this_dir = pathlib.Path(__file__).resolve().parent

logging.basicConfig(stream=sys.stdout, level=logging.INFO)


class WorkflowBenchmark:
    """Generate a workflow benchmark instance based on a workflow recipe (WfChefWorkflowRecipe)

    :param recipe: A workflow recipe.
    :type recipe: Type[WfChefWorkflowRecipe]
    :param num_tasks: Total number of tasks in the benchmark workflow.
    :type num_tasks: int
    :param with_flowcept:
    :type with_flowcept: bool
    :param logger: The logger where to log information/warning or errors.
    :type logger: Optional[Logger]
    """

    def __init__(self,
                 recipe: Type[WfChefWorkflowRecipe],
                 num_tasks: int,
                 with_flowcept: bool = False,
                 logger: Optional[Logger] = None) -> None:
        """Create an object that represents a workflow benchmark generator."""
        self.logger: Logger = logging.getLogger(
            __name__) if logger is None else logger
        self.recipe: Type[WfChefWorkflowRecipe] = recipe
        self.num_tasks = num_tasks
        self.with_flowcept = with_flowcept
        self.workflow: [Workflow|None] = None

    def create_benchmark_from_input_file(self,
                                         save_dir: pathlib.Path,
                                         input_file: pathlib.Path,
                                         lock_files_folder: Optional[pathlib.Path] = None,
                                         rundir:Optional[pathlib.Path] = None) -> pathlib.Path:
        """Create a workflow benchmark.

        :param save_dir: Folder to generate the workflow benchmark JSON instance and input data files.
        :type save_dir: pathlib.Path
        :param input_file: 
        :type input_file: pathlib.Path
        :param lock_files_folder:
        :type lock_files_folder: Optional[pathlib.Path]
        :param rundir: If you would like for the files to be created/saved in a different directory.
        :type rundir: Optional[pathlib.Path]

        :return: The path to the workflow benchmark JSON instance.
        :rtype: pathlib.Path
        """
        params = json.loads(input_file.read_text())
        return self.create_benchmark(save_dir, lock_files_folder=lock_files_folder, rundir=rundir, **params)

    def create_benchmark_from_synthetic_workflow(
            self,
            save_dir: pathlib.Path,
            workflow: Workflow,
            percent_cpu: Union[float, Dict[str, float]] = 0.6,
            cpu_work: Union[int, Dict[str, int]] = None,
            gpu_work: Union[int, Dict[str, int]] = None,
            num_chunks: Optional[int] = 10,
            time: Optional[int] = None,
            mem: Optional[float] = None,
            lock_files_folder: Optional[pathlib.Path] = None,
            rundir: Optional[pathlib.Path] = None) -> pathlib.Path:
        """Create a workflow benchmark from a synthetic workflow

        :param save_dir: Folder to generate the workflow benchmark JSON instance and input data files.
        :type save_dir: pathlib.Path
        :param workflow: The (synthetic) workflow to use as a benchmark.
        :type workflow: Workflow
        :param percent_cpu: The maximum percentage of CPU threads.
        :type percent_cpu: Union[float, Dict[str, float]]
        :param cpu_work: Maximum CPU work per workflow task.
        :type cpu_work: Union[int, Dict[str, int]]
        :param gpu_work: Maximum GPU work per workflow task.
        :type gpu_work: Union[int, Dict[str, int]]
        :param num_chunks: Number of chunks for pipelining I/O and computation for each task execution.
        :type num_chunks: Optional[int]
        :param time: Time limit for running each task (in seconds).
        :type time: Optional[int]
        :param mem: Maximum amount of memory consumption per task (in MB).
        :type mem: Optional[float]
        :param lock_files_folder:
        :type lock_files_folder: Optional[pathlib.Path]
        :param rundir: If you would like for the files to be created/saved in a different directory.
        :type rundir: Optional[pathlib.Path]

        :return: The path to the workflow benchmark JSON instance.
        :rtype: pathlib.Path
        """

        self.workflow = workflow
        # the benchmark tasks (i.e., their commands) are part of the execution section
        self.workflow.include_execution = True
        self.workflow.makespan = 0.0

        save_dir = save_dir.resolve()
        save_dir.mkdir(exist_ok=True, parents=True)

        json_path = save_dir.joinpath(
            f"{self.workflow.name.lower()}-{self.num_tasks}").with_suffix(".json")

        # tasks are grouped by category (or by name, as tasks often do not have a category)
        task_max_runtimes = {}
        for task in self.workflow.tasks.values():
            task_type = self._task_type(task)
            if task_type not in task_max_runtimes or task.runtime > task_max_runtimes[task_type]:
                task_max_runtimes[task_type] = task.runtime
        max_runtime = max(task_max_runtimes.values(), default=0)

        # if no cpu_work is provided, use the maximum runtime of each task as a reference
        if cpu_work is None:
            cpu_work = {task_type: runtime * 1000 for task_type, runtime in task_max_runtimes.items()}
            if not max_runtime:
                self.logger.warning("The workflow tasks have no runtimes (e.g., the workflow instance has no "
                                    "execution data): the benchmark tasks will not perform any CPU work "
                                    "unless cpu_work is provided.")

        cores, lock = self._creating_lock_files(lock_files_folder)

        for task in self.workflow.tasks.values():
            # all runtimes are zero if the workflow has no execution data
            runtime_factor = task.runtime / max_runtime if max_runtime else 1
            task_max_runtime = task_max_runtimes[self._task_type(task)]
            task_runtime_factor = task.runtime / task_max_runtime if task_max_runtime else 1
            # scale argument parameters to achieve a runtime distribution
            task_percent_cpu = self._task_value(percent_cpu, task) * task_runtime_factor if isinstance(percent_cpu, dict) else percent_cpu * runtime_factor
            task_cores = int(10 * task_percent_cpu)  # set number of cores to cpu threads in wfbench
            task_percent_cpu = max(0.1, task_percent_cpu)  # set minimum to 0.1 which is equivalent to 1 thread in wfbench
            task_percent_cpu = round(task_percent_cpu, 2)
            if cpu_work is not None:
                task_cpu_work = self._task_value(cpu_work, task) * task_runtime_factor if isinstance(cpu_work, dict) else cpu_work * runtime_factor
                task_cpu_work = int(task_cpu_work)
            else:
                task_cpu_work = None
            if gpu_work is not None:
                task_gpu_work = self._task_value(gpu_work, task) * task_runtime_factor if isinstance(gpu_work, dict) else gpu_work * runtime_factor
                task_gpu_work = int(task_gpu_work)
            else:
                task_gpu_work = None
            task_memory = int(mem * runtime_factor) if mem else None
            self._set_argument_parameters(
                task,
                percent_cpu=task_percent_cpu,
                cpu_work=task_cpu_work,
                gpu_work=task_gpu_work,
                time=time,
                num_chunks=num_chunks,
                mem=task_memory,
                lock_files_folder=lock_files_folder,
                cores=cores,
                lock=lock,
                rundir=rundir
            )
            task.cores = task_cores + 1
            if task_memory:
                task.memory = task_memory * 1024 * 1024  # megabytes to bytes

        # create data footprint
        for task in self.workflow.tasks.values():
            output_files = {file.file_id: file.size for file in task.output_files}
            task.args.append(f"--output-files {output_files}")

            input_files = [file.file_id for file in task.input_files]
            task.args.append(f"--input-files {input_files}")

        workflow_input_files: List[File] = self._rename_files_to_wfbench_format()

        for i, file in enumerate(workflow_input_files):
            file_path = save_dir.joinpath(file.file_id)
            if not file_path.is_file():
                print(
                    f"Creating {str(file_path)} ({file.size} bytes) ... file {i+1} out of {len(workflow_input_files)}",
                    end='\r'
                )
                with open(save_dir.joinpath("to_create.txt"), "a+") as fp:
                    fp.write(f"{file.file_id} {file.size}\n")
                self.logger.debug(f"Created file: {str(file_path)}")

        self.logger.info(f"Saving benchmark workflow: {json_path}")
        self.workflow.write_json(json_path)

        return json_path

    def _rename_files_to_wfbench_format(self) -> List[File]:
        """
        Rename the files in the workflow to the wfbench format.

        :return: A list of the input files that need to be generated (with their new names).
        :rtype: List[File]
        """
        new_file_names: Dict = {}
        renamed_files: Set[int] = set()
        task_output_counter = 0
        workflow_inputs: List[File] = []
    
        for task in self.workflow.tasks.values():
            # task IDs may contain characters that are not safe in file names (e.g., '/' or ':')
            file_prefix = re.sub(r"[^0-9A-Za-z_.-]", "_", task.task_id)
            output_files = sorted(task.output_files, key=lambda x: -len(x.file_id))
            for file in output_files:
                if file.file_id in new_file_names:
                    raise ValueError(f"File name {file.file_id} already exists")
                
                task_output_counter += 1
                # extension = ''.join(pathlib.Path(file.file_id).suffixes)
                new_name = f"{file_prefix}_outfile_{task_output_counter:04d}" #{extension}
                new_file_names[file.file_id] = new_name
                file.file_id = new_name
                renamed_files.add(id(file))

        for task in self.workflow.tasks.values():
            input_files = sorted(task.input_files, key=lambda x: -len(x.file_id))
            for file in input_files:
                if id(file) in renamed_files:
                    # file object shared with the task that produces it, already renamed
                    continue
                if file.file_id in new_file_names:
                    # file is an output file of another task (or an input file for the workflow
                    # used by another task) and receives the corresponding name
                    file.file_id = new_file_names[file.file_id]
                else:
                    # file is an input file for the workflow and needs to be generated
                    workflow_inputs.append(file)
                    # extension = ''.join(pathlib.Path(file.file_id).suffixes)
                    new_name = f"workflow_infile_{len(workflow_inputs):04d}"#{extension}
                    new_file_names[file.file_id] = new_name
                    file.file_id = new_name
                renamed_files.add(id(file))

        # rename the files in the input/output files arguments of the tasks (only)
        for task in self.workflow.tasks.values():
            task.args = [self._rename_files_in_argument(arg, new_file_names) for arg in task.args]
    
        return workflow_inputs

    @staticmethod
    def _rename_files_in_argument(argument: str, new_file_names: Dict[str, str]) -> str:
        """
        Rename the files listed in an '--input-files' or '--output-files' task argument.

        :param argument: A task argument.
        :type argument: str
        :param new_file_names: The new names of the files, indexed by their current names.
        :type new_file_names: Dict[str, str]

        :return: The argument with the renamed files (other arguments are returned unchanged).
        :rtype: str
        """
        if argument.startswith("--output-files "):
            output_files = ast.literal_eval(argument[len("--output-files "):])
            output_files = {new_file_names.get(file_name, file_name): size for file_name, size in output_files.items()}
            return f"--output-files {output_files}"
        if argument.startswith("--input-files "):
            input_files = ast.literal_eval(argument[len("--input-files "):])
            input_files = [new_file_names.get(file_name, file_name) for file_name in input_files]
            return f"--input-files {input_files}"
        return argument

    @staticmethod
    def _task_type(task: Task) -> str:
        """
        Get the type of a task, i.e., its category, or its name if the task has no category.

        :param task: The task.
        :type task: Task

        :return: The task type.
        :rtype: str
        """
        return task.category if task.category is not None else task.name

    @staticmethod
    def _task_value(values: Union[float, Dict[str, float]], task: Task) -> float:
        """
        Get the value of a benchmark parameter for a task, if values are provided per task
        category (or per task name, as tasks often do not have a category).

        :param values: A value, or values indexed by task category or task name.
        :type values: Union[float, Dict[str, float]]
        :param task: The task.
        :type task: Task

        :return: The value of the parameter for the task.
        :rtype: float
        """
        if not isinstance(values, dict):
            return values
        if task.category is not None and task.category in values:
            return values[task.category]
        if task.name in values:
            return values[task.name]
        raise KeyError(f"No value provided for task '{task.task_id}' (category: '{task.category}', "
                       f"name: '{task.name}').")

    def create_benchmark(self,
                         save_dir: pathlib.Path,
                         percent_cpu: Union[float, Dict[str, float]] = 0.6,
                         cpu_work: Union[int, Dict[str, int]] = None,
                         gpu_work: Union[int, Dict[str, int]] = None,
                         num_chunks: Optional[int] = 10,
                         time: Optional[int] = None,
                         data: Optional[int] = 0,
                         mem: Optional[float] = None,
                         lock_files_folder: Optional[pathlib.Path] = None,
                         regenerate: Optional[bool] = True,
                         rundir: Optional[pathlib.Path] = None,
                         ) -> pathlib.Path:
        """Create a workflow benchmark.

        :param save_dir: Folder to generate the workflow benchmark JSON instance and input data files.
        :type save_dir: pathlib.Path
        :param percent_cpu: The percentage of CPU threads.
        :type percent_cpu: Union[float, Dict[str, float]]
        :param cpu_work: CPU work per workflow task.
        :type cpu_work: Union[int, Dict[str, int]]
        :param gpu_work: GPU work per workflow task.
        :type gpu_work: Union[int, Dict[str, int]]
        :param num_chunks: Number of chunks for pipelining I/O and computation for each task execution.
        :type num_chunks: Optional[int]
        :param time: Time limit for running each task (in seconds).
        :type time: Optional[int]
        :param data: Total workflow data footprint (in MB).
        :type data: Optional[Union[int, Dict[str, str]]]
        :param mem: Maximum amount of memory consumption per task (in MB).
        :type mem: Optional[float]
        :param lock_files_folder:
        :type lock_files_folder: Optional[pathlib.Path]
        :param regenerate: Whether to regenerate the workflow tasks
        :type regenerate: Optional[bool]
        :param rundir: If you would like for the files to be created/saved in a different directory.
        :type rundir: Optional[pathlib.Path]

        :return: The path to the workflow benchmark JSON instance.
        :rtype: pathlib.Path
        """
        save_dir = save_dir.resolve()
        save_dir.mkdir(exist_ok=True, parents=True)

        if not self.workflow or regenerate:
            self.logger.debug("Generating workflow")
            generator = WorkflowGenerator(
                self.recipe.from_num_tasks(self.num_tasks))
            self.workflow = generator.build_workflow()
            self.workflow.name = f"{self.workflow.name.split('-')[0]}-Benchmark"
        # the benchmark tasks (i.e., their commands) are part of the execution section
        self.workflow.include_execution = True
        self.workflow.makespan = 0.0
        json_path = save_dir.joinpath(
            f"{self.workflow.name.lower()}-{self.num_tasks}").with_suffix(".json")

        if self.with_flowcept:
            self.workflow.workflow_id = str(shortuuid.uuid())

        cores, lock = self._creating_lock_files(lock_files_folder)
        for task in self.workflow.tasks.values():
            self._set_argument_parameters(
                task,
                percent_cpu=percent_cpu,
                cpu_work=cpu_work,
                gpu_work=gpu_work,
                time=time,
                num_chunks=num_chunks,
                mem=mem,
                lock_files_folder=lock_files_folder,
                cores=cores,
                lock=lock,
                rundir=rundir,
            )
            task.input_files = []
            task.output_files = []
        
        self._create_data_footprint(data)
        
        # TODO: add a flag to allow the file names to be changed 
        workflow_input_files: List[File] = self._rename_files_to_wfbench_format()

        for i, file in enumerate(workflow_input_files):
            file_path = save_dir.joinpath(file.file_id)
            if not file_path.is_file():
                print(
                    f"Creating {str(file_path)} ({file.size} bytes) ... file {i+1} out of {len(workflow_input_files)}",
                    end='\r'
                )
                with open(save_dir.joinpath("created_input_files.txt"), "a+") as fp:
                    fp.write(f"{file.file_id} {file.size}\n")
                self.logger.debug(f"Created file: {str(file_path)}")

                with open(save_dir.joinpath(file.file_id), 'wb') as fp:
                    fp.write(os.urandom(file.size))
        
        self.logger.info(f"Saving benchmark workflow: {json_path}")
        self.workflow.write_json(json_path)

        return json_path

    def _creating_lock_files(self, lock_files_folder: Optional[pathlib.Path]) -> Tuple[pathlib.Path | None, pathlib.Path | None]:
        """
        Creating the lock files
        """
        if not lock_files_folder:
            return None, None
        try:
            lock_files_folder.mkdir(exist_ok=True, parents=True)
            self.logger.debug(
                f"Creating lock files at: {lock_files_folder.resolve()}")
            lock = lock_files_folder.joinpath("cores.txt.lock")
            cores = lock_files_folder.joinpath("cores.txt")
            with lock.open("w+"), cores.open("w+"):
                pass
            return lock, cores
        except (FileNotFoundError, OSError) as e:
            self.logger.warning(f"Could not find folder to create lock files: {lock_files_folder.resolve()}\n"
                                f"You will need to create them manually: 'cores.txt.lock' and 'cores.txt'")
            return None, None

    def _set_argument_parameters(self,
                                 task: Task,
                                 percent_cpu: Union[float, Dict[str, float]],
                                 cpu_work: Union[int, Dict[str, int]],
                                 gpu_work: Union[int, Dict[str, int]],
                                 time: Optional[int],
                                 num_chunks: Optional[int],
                                 mem: Optional[float],
                                 lock_files_folder: Optional[pathlib.Path],
                                 cores: Optional[pathlib.Path],
                                 lock: Optional[pathlib.Path],
                                 rundir: Optional[pathlib.Path]) -> None:
        """
        Setting the parameters for the arguments section of the JSON
        """
        params = []

        cpu_params = self._generate_task_cpu_params(task, percent_cpu, cpu_work, lock_files_folder, cores, lock)
        params.extend(cpu_params)
        gpu_params = self._generate_task_gpu_params(task, gpu_work)
        params.extend(gpu_params)
        params.extend([f"--num-chunks {num_chunks}"])

        if mem:
            params.extend([f"--mem {mem}"])

        if time:
            params.extend([f"--time {time}"])

        if rundir:
            params.extend([f"--rundir {rundir}"])

        if self.with_flowcept:
            params.extend(["--with-flowcept"])

            if self.workflow.workflow_id:
                params.extend([f"--workflow_id {self.workflow.workflow_id}"])

        # the benchmark has not been executed yet: clear the measured execution data of the task
        task.runtime = 0
        task.avg_cpu = None
        task.bytes_read = None
        task.bytes_written = None
        task.memory = None
        task.energy = None
        task.avg_power = None
        task.machines = None

        task.program = "wfbench"
        task.args = [f"--name {task.task_id}"]
        task.args.extend(params)

    def _generate_task_cpu_params(self,
                                  task: Task,
                                  percent_cpu: Union[float, Dict[str, float]],
                                  cpu_work: Union[int, Dict[str, int]],
                                  lock_files_folder: Optional[pathlib.Path],
                                  cores: Optional[pathlib.Path],
                                  lock: Optional[pathlib.Path]) -> List[str]:
        """
        Setting cpu arguments if cpu benchmark requested
        """
        if not cpu_work:
            return []

        _percent_cpu = self._task_value(percent_cpu, task)
        _cpu_work = self._task_value(cpu_work, task)

        params = [f"--percent-cpu {_percent_cpu}", f"--cpu-work {int(_cpu_work)}"]

        if lock_files_folder:
            params.extend([f"--path-lock {lock}",
                           f"--path-cores {cores}"])
        return params

    def _generate_task_gpu_params(self, task: Task, gpu_work: Union[int, Dict[str, int]]) -> List[str]:
        """
        Setting gpu arguments if gpu benchmark requested
        """
        if not gpu_work:
            return []
        _gpu_work = self._task_value(gpu_work, task)

        return [f"--gpu-work {_gpu_work}"]


    def _create_data_footprint(self, data: int) -> None:
        """
        task's data footprint provided as an int
        """
        num_sys_files, num_total_files = self._calculate_input_files()
        self.logger.debug(
            f"Number of input files to be created by the system: {num_sys_files}")
        self.logger.debug(
            f"Total number of files used by the workflow: {num_total_files}")
        file_size = round(data * 1000000 / num_total_files)  # MB to B
        self.logger.debug(
            f"Every input/output file is of size: {file_size}")

        for task in self.workflow.tasks.values():
            output = {f"{task.task_id}_output.txt": file_size}
            task.args.extend([f"--output-files {output}"])

        self._add_output_files(file_size)
        self._add_input_files(file_size)


    def _output_files(self, data: Dict[str, str]) -> Dict[str, Dict[str, int]]:
        """
        Calculate, for each task, total number of output files needed.
        This method is used when the user is specifying the input file sizes.

        :param data:
        :type data: Dict[str, str]

        :return: 
        :rtype: Dict[str, Dict[str, int]]
        """
        output_files = {}
        for task in self.workflow.tasks.values():
            output_files.setdefault(task.task_id, {})
            if not self.workflow.tasks_children[task.task_id]:
                output_files[task.task_id][task.task_id] = int(data[task.category])
            else:
                for child_name in self.workflow.tasks_children[task.task_id]:
                    child = self.workflow.tasks[child_name]
                    output_files[task.task_id][child.task_id] = int(
                        data[child.category])

        return output_files

    def _calculate_input_files(self):
        """
        Calculate total number of files needed.
        This method is used if the user provides total data footprint.
        """
        tasks_need_input = 0
        tasks_dont_need_input = 0

        for task in self.workflow.tasks.values():
            parents = self.workflow.tasks_parents[task.task_id]
            if not parents:
                tasks_need_input += 1
            else:
                tasks_dont_need_input += 1

        total_num_files = tasks_need_input * 2 + tasks_dont_need_input

        return tasks_need_input, total_num_files

    def _add_output_files(self, output_file_size: int) -> None:
        """
        Add output files when input data was offered by the user.

        :param output_file_size: file size in MB
        :type output_file_size: int
        """
        for task in self.workflow.tasks.values():
            task.output_files.append(
                File(f"{task.task_id}_output.txt", output_file_size))

    def _add_input_files(self, input_file_size: int) -> None:
        """
        Add input files when input data was offered by the user.

        :param input_file_size: a file size in MB
        :type input_file_size: int
        """
        for task in self.workflow.tasks.values():
            inputs = []
            if not self.workflow.tasks_parents[task.task_id]:
                task.input_files.append(
                    File(f"{task.task_id}_input.txt", input_file_size))
                inputs.append(f'{task.task_id}_input.txt')
            else:
                for parent in self.workflow.tasks_parents[task.task_id]:
                    task.input_files.append(
                        File(f"{parent}_output.txt", input_file_size))
                    inputs.append(f"{parent}_output.txt")

            task.args.append(f"--input-files {inputs}")

    # def _generate_data_for_root_nodes(self, save_dir: pathlib.Path, data: Union[int, Dict[str, str]]) -> None:
    #     """
    #     Generate workflow's input data for root nodes based on user's input.
    #
    #     :param save_dir:
    #     :type save_dir: pathlib.Path
    #     :param data:
    #     :type data: Dict[str, str]
    #     """
    #     for task in self.workflow.tasks.values():
    #         if not self.workflow.tasks_parents[task.task_id]:
    #             file_size = data[task.category] if isinstance(
    #                 data, Dict) else data
    #             file = save_dir.joinpath(f"{task.task_id}_input.txt")
    #             if not file.is_file():
    #                 with open(file, 'wb') as fp:
    #                     fp.write(os.urandom(int(file_size)))
    #                 self.logger.debug(f"Created file: {str(file)}")

    def generate_input_file(self, path: pathlib.Path) -> None:
        """
        Generates input file where customization of cpu percentage, cpu work, gpu work, data size

        :param path:
        :type path: pathlib.Path
        """
        generator = WorkflowGenerator(
            self.recipe.from_num_tasks(self.num_tasks))
        workflow = generator.build_workflow()

        defaults = {
            "percent_cpu": 0.6,
            "cpu_work": 1000,
            "gpu_work": 100,
            "data": 10
        }
        inputs = {
            "percent_cpu": {},
            "cpu_work": {},
            "gpu_work": {},
            "data": {}

        }
        for node in workflow.nodes:
            task: Task = workflow.nodes[node]['task']
            task_type = task.name

            for key in inputs.keys():
                inputs[key].setdefault(task_type, defaults[key])

        path.parent.mkdir(exist_ok=True, parents=True)
        path.write_text(json.dumps(inputs, indent=2))
        input("Please fill up the input file and press ENTER to continue...")


def generate_sys_data(num_files: int, tasks: Dict[str, int], save_dir: pathlib.Path) -> List[str]:
    """Generate workflow's input data

    :param num_files: number of each file to be generated.
    :type num_files: int
    :param tasks: Dictionary with the name of the tasks and their data sizes.
    :type tasks: Dict[str, int]
    :param save_dir: Folder to generate the workflow benchmark's input data files.
    :type save_dir: pathlib.Path
    """
    names = []
    for _ in range(num_files):
        for name, size in tasks.items():
            # name = f'{name}_input.txt'
            names.append(name)
            file = f"{save_dir.joinpath(name)}"
            with open(file, 'wb') as fp:
                fp.write(os.urandom(size))
            print(f"Created file: {file}")

    return names

def cleanup_sys_files() -> None:
    """Remove files already used"""
    input_files = glob.glob("*input*.txt")
    output_files = glob.glob("*output.txt")
    all_files = input_files + output_files
    for t in all_files:
        os.remove(t)


# Function to clean and adjust the list entries
def clean_entry(entry):

    if entry.startswith('--out '):
        # Replace --out "..." with --out=...
        return entry.replace('--out "', '--out=').replace('}"', '}')
    else:
        # Remove extra double quotes
        entry = entry.replace(' ', '=')
        return entry.strip('"')
