#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging
import os
import re
import sys
import site
import pathlib
import shutil
import textwrap
from abc import ABC, abstractmethod
from typing import Dict, Iterable, List, Optional, Tuple, Union

from ...common import File, Task, Workflow
from ...wfinstances.instance import Instance


this_dir = pathlib.Path(__file__).resolve().parent

class Translator(ABC):
    """
    An abstract class of WfFormat parser for creating workflow benchmark applications.

    :param workflow: Workflow benchmark object or path to the workflow benchmark JSON instance.
    :type workflow: Union[Workflow, pathlib.Path]
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: logging.Logger
    """

    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 logger: Optional[logging.Logger] = None) -> None:
        """Create an object of the translator."""
        self.logger = logging.getLogger(__name__) if logger is None else logger
        
        if isinstance(workflow, Workflow):
            self.workflow = workflow
        else:
            instance = Instance(workflow, logger=logger)
            self.workflow: Workflow = instance.workflow

        # find all tasks
        self.tasks = {}
        for task_id, task in self.workflow.nodes(data="task"):
            self.tasks[task_id] = task

        # find root, parents, and children tasks (parents and children are listed in edge order)
        self.task_parents = {task_id: [] for task_id in self.tasks}
        self.task_children = {task_id: [] for task_id in self.tasks}
        for parent, child in self.workflow.edges:
            self.task_children[parent].append(child)
            self.task_parents[child].append(parent)
        self.root_task_names = [task_id for task_id in self.tasks if not self.task_parents[task_id]]

        self._check_benchmark_commands()
        self._dependency_markers = None

    @abstractmethod
    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into an actual workflow application.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """

    def _check_benchmark_commands(self) -> None:
        """
        Check that every task has a command to run. The command is optional in WfFormat, but
        translators generate applications that run the (WfBench) command of each task.

        :raises ValueError: If some tasks do not have a command.
        """
        missing = [task_id for task_id, task in self.tasks.items() if not task.program]
        if missing:
            listed = ", ".join(f"'{task_id}'" for task_id in missing[:5])
            if len(missing) > 5:
                listed += f", ... ({len(missing)} tasks in total)"
            raise ValueError(f"Cannot translate workflow '{self.workflow.name}': task(s) {listed} have no command "
                             f"to run. Translators require a workflow benchmark, in which every task runs "
                             f"'wfbench'. Create one with WorkflowBenchmark (e.g., "
                             f"WorkflowBenchmark.create_benchmark_from_synthetic_workflow() for an existing "
                             f"workflow or instance), and translate the generated benchmark instead.")

    @staticmethod
    def _unique_identifiers(names: Dict[str, str],
                            invalid_chars: str = r"[^0-9A-Za-z_]",
                            first_char: Optional[str] = r"[A-Za-z_]",
                            reserved: Iterable[str] = ()) -> Dict[str, str]:
        """
        Map names (e.g., task IDs or task names, which may contain characters such as ':', '/',
        '#', '.', or '-') to unique identifiers that are valid in the target language. A name that
        is already a valid identifier is kept unchanged, so that the translation of a workflow with
        well-behaved names does not change. Otherwise, each invalid character is replaced by '_',
        the identifier is prefixed with '_' if it does not start with a valid character, and a
        numeric suffix is appended if it collides with another identifier.

        :param names: The names to map, indexed by a unique key (e.g., the task ID).
        :type names: Dict[str, str]
        :param invalid_chars: Regular expression matching the characters not allowed in an identifier.
        :type invalid_chars: str
        :param first_char: Regular expression matching the allowed first characters (None for any).
        :type first_char: Optional[str]
        :param reserved: Identifiers that must not be used (e.g., keywords of the target language).
        :type reserved: Iterable[str]

        :return: The identifier of each name, indexed by the same keys.
        :rtype: Dict[str, str]
        """
        invalid_re = re.compile(invalid_chars)
        first_re = re.compile(first_char) if first_char else None
        reserved = set(reserved)
        identifiers = {}
        used = set()

        def _is_valid(identifier: str) -> bool:
            return (identifier != "" and not invalid_re.search(identifier) and identifier not in reserved
                    and (first_re is None or first_re.match(identifier) is not None))

        # valid names keep their own identifier (unless already taken)
        for key, name in names.items():
            name = str(name)
            if _is_valid(name) and name not in used:
                identifiers[key] = name
                used.add(name)

        # other names are sanitized and disambiguated deterministically
        for key, name in names.items():
            if key in identifiers:
                continue
            identifier = invalid_re.sub("_", str(name))
            if identifier == "" or (first_re is not None and first_re.match(identifier) is None):
                identifier = f"_{identifier}"
            candidate = identifier
            suffix = 1
            while candidate in used or candidate in reserved:
                candidate = f"{identifier}_{suffix}"
                suffix += 1
            identifiers[key] = candidate
            used.add(candidate)

        return identifiers

    @staticmethod
    def _shell_safe_argument(argument: str) -> str:
        """
        Make a task argument safe to be used as a word in a shell command. The '--name' value of a
        WfBench task is its task ID, which may start with '#' (a shell comment) or '-' (an option
        for the argument parser). Such a value is attached to the flag ('--name=<value>'), which
        does not require any quoting and therefore works in every generated language.

        :param argument: A task argument.
        :type argument: str

        :return: The argument to use in a shell command.
        :rtype: str
        """
        if argument.startswith("--name ") and argument[len("--name "):].startswith(("#", "-")):
            return f"--name={argument[len('--name '):]}"
        return argument

    def _dependencies_without_files(self) -> List[Tuple[str, str]]:
        """
        Find the task dependencies that are not expressed by files, i.e., the (parent, child)
        dependencies for which the child does not consume any file produced by the parent.

        :return: List of (parent, child) task IDs.
        :rtype: List[Tuple[str, str]]
        """
        dependencies = []
        for parent_id, children in self.task_children.items():
            parent_outputs = {file.file_id for file in self.tasks[parent_id].output_files}
            for child_id in children:
                if not parent_outputs.intersection(file.file_id for file in self.tasks[child_id].input_files):
                    dependencies.append((parent_id, child_id))
        return dependencies

    def _get_dependency_markers(self) -> Tuple[Dict[str, List[File]], Dict[str, List[File]]]:
        """
        For workflow systems that infer task dependencies from files only, create empty marker
        files that express the dependencies that are not expressed by files (e.g., in workflows
        without files). A marker file is an additional (zero-size) output file of the parent task,
        and an additional input file of its children. A task without output files also gets a
        marker file, so that it is a target of the workflow. Workflows in which all dependencies
        are expressed by files do not get any marker file.

        :return: The additional input files and the additional output files, indexed by task ID.
        :rtype: Tuple[Dict[str, List[File]], Dict[str, List[File]]]
        """
        if self._dependency_markers is None:
            dependencies = self._dependencies_without_files()
            marked_tasks = {parent_id for parent_id, _ in dependencies}
            marked_tasks.update(task_id for task_id, task in self.tasks.items() if not task.output_files)
            markers = self._create_marker_files(marked_tasks) if marked_tasks else {}

            extra_inputs = {task_id: [] for task_id in self.tasks}
            for parent_id, child_id in dependencies:
                extra_inputs[child_id].append(markers[parent_id])
            extra_outputs = {task_id: [markers[task_id]] if task_id in markers else [] for task_id in self.tasks}
            if markers:
                self.logger.info(f"Using {len(markers)} empty marker file(s) to express task dependencies "
                                 f"that are not expressed by files.")
            self._dependency_markers = (extra_inputs, extra_outputs)

        return self._dependency_markers

    def _create_marker_files(self, task_ids: Iterable[str]) -> Dict[str, File]:
        """
        Create a (zero-size) dependency marker file for each task, with a unique file ID.

        :param task_ids: The IDs of the tasks.
        :type task_ids: Iterable[str]

        :return: The marker file of each task, indexed by task ID.
        :rtype: Dict[str, File]
        """
        file_ids = {file.file_id for task in self.tasks.values() for file in task.input_files + task.output_files}
        safe_names = self._unique_identifiers({task_id: task_id for task_id in self.tasks},
                                              invalid_chars=r"[^0-9A-Za-z_.-]", first_char=None)
        markers = {}
        for task_id in task_ids:
            marker_id = f"{safe_names[task_id]}.wfbench_done"
            while marker_id in file_ids:
                marker_id = f"_{marker_id}"
            markers[task_id] = File(marker_id, 0)
        return markers

    def _task_input_files(self, task: Task) -> List[File]:
        """
        Get the input files of a task, including dependency marker files (see _get_dependency_markers).

        :param task: The task.
        :type task: Task

        :return: The input files of the task.
        :rtype: List[File]
        """
        return task.input_files + self._get_dependency_markers()[0][task.task_id]

    def _task_output_files(self, task: Task) -> List[File]:
        """
        Get the output files of a task, including dependency marker files (see _get_dependency_markers).

        :param task: The task.
        :type task: Task

        :return: The output files of the task.
        :rtype: List[File]
        """
        return task.output_files + self._get_dependency_markers()[1][task.task_id]

    def _warn_dependencies_without_files(self, system_name: str) -> None:
        """
        Warn that the task dependencies not expressed by files cannot be enforced by a workflow
        system that infers task dependencies from files.

        :param system_name: The name of the workflow system.
        :type system_name: str
        """
        dependencies = self._dependencies_without_files()
        if dependencies:
            examples = ", ".join(f"'{parent}' -> '{child}'" for parent, child in dependencies[:3])
            self.logger.warning(f"{system_name} infers task dependencies from files: {len(dependencies)} "
                                f"task dependencies that are not expressed by files (e.g., {examples}) "
                                f"will not be enforced in the generated workflow.")

    def _copy_binary_files(self, output_folder: pathlib.Path) -> None:
        """
        Copy binary files to workflow benchmark's bin folder.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        bin_folder = output_folder.joinpath("bin")
        bin_folder.mkdir(exist_ok=True)

        shutil.copy(shutil.which("wfbench"), bin_folder)
        # Fix the shebang in wfbench
        with open(bin_folder / "wfbench", "r") as f:
            lines = f.readlines()

        if lines and lines[0].startswith("#!"):
            lines[0] = "#!/usr/bin/env python3\n"
            with open(bin_folder / "wfbench", "w") as f:
                f.writelines(lines)


    def _generate_input_files(self, output_folder: pathlib.Path) -> None:
        """
        Generate workflow input files into workflow benchmark's data folder.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        workflow_input_files = []
        input_files = set()
        output_files = set()

        for task in self.workflow.tasks.values():
            for file in task.input_files:
                input_files.add(file)
            for file in task.output_files:
                output_files.add(file)
        
        workflow_input_files = input_files - output_files

        data_folder = output_folder.joinpath("data")
        data_folder.mkdir(exist_ok=True)
        for file in workflow_input_files:
            # file IDs may contain '/'
            file_path = data_folder.joinpath(file.file_id)
            file_path.parent.mkdir(parents=True, exist_ok=True)
            with open(file_path, "wb") as fp:
                fp.write(os.urandom(int(file.size)))

        # wfbench does not create the folders of the output files it writes
        for file in output_files:
            data_folder.joinpath(file.file_id).parent.mkdir(parents=True, exist_ok=True)

    def _write_output_file(self, contents: str, output_file_path: pathlib.Path) -> None:
        """
        Write the translated content to a file.

        :param contents: Contents to be written to the output file.
        :type contents: str
        :param output_file_path: The path of the output file.
        :type output_file_path: pathlib.Path
        """
        # file will be written to the same folder as for the original JSON instance.
        with open(output_file_path, "w") as out:
            out.write(contents)
        self.logger.info(f"Translated content written to '{output_file_path}'")

    def _find_root_tasks(self) -> list[Task]:
        """
        Find the workflow's root (i.e., parentless) tasks.

        :return: List of task
        :rtype: list[Task]
        """
        return [self.tasks[task_name] for task_name in self.root_task_names]

    def _find_children(self, task_name: str) -> list[Task]:
        """
        Find the children for a specific task.

        :param task_name: The task name.
        :type task_name: str

        :return: List of task's children.
        :rtype: list[Task]
        """
        self.logger.debug(f"Finding children for task '{task_name}'")
        return [self.tasks[task_id] for task_id in self.task_children[task_name]]

    def _find_parents(self, task_name: str) -> list[Task]:
        """
        Find the parents for a specific task.

        :param task_name: The task name.
        :type task_name: str

        :return: List of task's parents.
        :rtype: list[Task]
        """
        self.logger.debug(f"Finding parents for task '{task_name}'")
        return [self.tasks[task_id] for task_id in self.task_parents[task_name]]

    def _merge_codelines(self, template_file_path: str, wf_codelines: str) -> str:
        """
        Incorporate generated workflow codelines into the template.

        :param template_file_path: The path to the template file.
        :type template_file_path: str

        :param wf_codelines: The generated workflow codelines.
        :type wf_codelines: str
        
        :return: Incorporated workflow codelines.
        :rtype: str
        """
        with open(this_dir.joinpath(template_file_path)) as fp:
            run_workflow_code = fp.read()
            return run_workflow_code.replace("# Generated code goes here", wf_codelines)

    def _flowcept_init_python(self, workflow_id: str, workflow_name: str) -> str:
        """

        :param workflow_id: 
        :type workflow_id: str

        :param workflow_name:
        :type workflow_name: str

        :return:
        :rtype: str
        """
        code = textwrap.dedent(f"""
        from flowcept.flowcept_api.flowcept_controller import Flowcept
        flowcept_agent = Flowcept(workflow_id="{workflow_id}", workflow_name="{workflow_name}", bundle_exec_id="{workflow_id}", start_persistence=False)
        flowcept_agent.start()
        """)
        return code

    def _flowcept_stop_python(self) -> str:
        """
        
        :return: 
        :rtype: str
        """
        return "flowcept_agent.stop()"
