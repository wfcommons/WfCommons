#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2026 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import jsonschema
import logging
import math
import pathlib

from logging import Logger
from typing import Any, Dict, List, Optional, Set

from ..version import __schema_version__

SCHEMA_FILE_NAME = "wfcommons-schema.json"
SCHEMA_URL = f"https://raw.githubusercontent.com/wfcommons/wfformat/main/{SCHEMA_FILE_NAME}"
SCHEMA_FETCH_TIMEOUT_IN_SECONDS = 30

# tolerance used when comparing declared metrics against recomputed ones
METRICS_RELATIVE_TOLERANCE = 1e-9


class SchemaValidator:
    """
    Validate JSON files against WfCommons schema (WfFormat). If schema file path
    is not provided, it will use the copy of the WfFormat schema distributed with
    WfCommons, which matches the schema version written by WfCommons.

    :param schema_file_path: JSON schema file path.
    :type schema_file_path: Optional[pathlib.Path]
    :param logger: The logger where to log information/warning or errors.
    :type logger: Optional[Logger]
    """

    def __init__(self, schema_file_path: Optional[pathlib.Path] = None, logger: Optional[Logger] = None) -> None:
        """Create an object of the schema validator class."""
        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        self.schema = self._load_schema(schema_file_path)
        self.errors: List[str] = []
        self.warnings: List[str] = []

    def validate_instance(self, data: Dict[str, Any]) -> None:
        """
        Perform syntax validation against the schema, and semantic validation.

        :param data: Workflow instance in JSON format.
        :type data: Dict[str, Any]
        """
        self.errors = []
        self.warnings = []
        self._syntax_validation(data)
        self._semantic_validation(data)

    def _load_schema(self, schema_file_path: Optional[pathlib.Path] = None) -> Dict[str, Any]:
        """
        Load the schema file. If schema file path is not provided, it will use the copy of
        the WfFormat schema distributed with WfCommons, and if not available it will fetch
        the latest schema from the GitHub repository.

        :param schema_file_path: JSON schema file path.
        :type schema_file_path: Optional[pathlib.Path]

        :return: The JSON schema.
        :rtype: Dict[str, Any]
        """
        if schema_file_path:
            self.logger.info(f'Using schema file: {schema_file_path}')
            with open(schema_file_path, encoding="utf-8") as f:
                return json.load(f)

        # copy of the schema distributed with WfCommons
        schema_path = pathlib.Path(__file__).parent.joinpath(SCHEMA_FILE_NAME)
        if schema_path.exists():
            self.logger.debug(f'Using schema file: {schema_path}')
            with open(schema_path, encoding="utf-8") as f:
                return json.load(f)

        # fetching latest schema file from GitHub repository
        import requests

        response = requests.get(SCHEMA_URL, timeout=SCHEMA_FETCH_TIMEOUT_IN_SECONDS)
        response.raise_for_status()
        schema = response.json()
        self.logger.info(f"Using latest schema file from GitHub repository: {SCHEMA_URL}")
        if __schema_version__ not in schema.get("properties", {}).get("schemaVersion", {}).get("enum", []):
            self.logger.warning(f"The latest schema does not accept version {__schema_version__}, "
                                f"which is the version written by WfCommons.")
        return schema

    def _error(self, message: str) -> None:
        self.errors.append(message)
        self.logger.error(message)

    def _warning(self, message: str) -> None:
        self.warnings.append(message)
        self.logger.warning(message)

    def _syntax_validation(self, data: Dict[str, Any]) -> None:
        """
        Validate the JSON workflow execution instance against the schema.

        :param data: Workflow instance in JSON format.
        :type data: Dict[str, Any]
        """
        # honor the draft the schema declares, rather than assuming one
        validator_class = jsonschema.validators.validator_for(self.schema)
        validator = validator_class(self.schema, format_checker=jsonschema.FormatChecker())

        has_error = False
        for error in sorted(validator.iter_errors(data), key=str):
            location = ' > '.join([str(e) for e in error.absolute_path]) or "<instance>"
            msg = f"{location}: {error.message}"
            # "format" is an annotation rather than an assertion in JSON Schema
            if error.validator == "format":
                self._warning(msg)
            else:
                self._error(msg)
                has_error = True

        if has_error:
            raise RuntimeError('The workflow instance has syntax errors.')

    def _semantic_validation(self, data: Dict[str, Any]) -> None:
        """
        Validate the semantics of the JSON workflow execution instance.

        :param data: Workflow instance in JSON format.
        :type data: Dict[str, Any]
        """
        errors_before = len(self.errors)

        workflow = data["workflow"]
        specification = workflow["specification"]
        execution = workflow.get("execution")

        tasks = specification["tasks"]
        task_ids = self._collect_ids(tasks, "id", "task")
        file_ids = self._collect_ids(specification.get("files", []), "id", "file")

        self._validate_dependencies(tasks, task_ids)
        self._validate_file_references(tasks, specification, file_ids)
        level_widths = self._compute_level_widths(tasks)
        self._validate_specification_metrics(specification, level_widths)

        if execution is None:
            self.logger.debug("Skipping execution processing.")
        else:
            self._validate_execution(execution, task_ids)
            self._validate_execution_metrics(execution)

        self.logger.debug(f'The workflow has {len(task_ids)} tasks.')
        self.logger.debug(f'The workflow has {len(file_ids)} files.')

        if len(self.errors) > errors_before:
            raise RuntimeError('The workflow instance has semantic errors.')

    def _collect_ids(self, items: List[Dict[str, Any]], key: str, kind: str) -> Set[str]:
        identifiers = set()
        for item in items:
            identifier = item.get(key)
            if identifier is None:
                continue
            if identifier in identifiers:
                self._error(f"The {kind} \"{identifier}\" is declared more than once.")
            identifiers.add(identifier)
        return identifiers

    def _validate_dependencies(self, tasks: List[Dict[str, Any]], task_ids: Set[str]) -> None:
        parents_of = {}
        children_of = {}
        for task in tasks:
            task_id = task["id"]
            for kind, key in (("parent", "parents"), ("child", "children")):
                references = task.get(key, [])
                for reference in references:
                    if reference == task_id:
                        self._error(f"Task \"{task_id}\" declares itself as its own {kind}.")
                    elif reference not in task_ids:
                        self._error(f"{kind.capitalize()} task \"{reference}\" is not declared in the list "
                                    f"of workflow tasks.")
                if len(references) != len(set(references)):
                    self._warning(f"Task \"{task_id}\" declares the same {kind} task more than once.")
            parents_of[task_id] = set(task.get("parents", []))
            children_of[task_id] = set(task.get("children", []))

        # both ends of every dependency have to agree
        for relation, reverse, kind, reverse_kind in ((parents_of, children_of, "parent", "child"),
                                                      (children_of, parents_of, "child", "parent")):
            for task_id, references in relation.items():
                for reference in references:
                    if reference != task_id and reference in reverse and task_id not in reverse[reference]:
                        self._warning(f"Task \"{task_id}\" declares \"{reference}\" as a {kind}, but task "
                                      f"\"{reference}\" does not declare \"{task_id}\" as a {reverse_kind}.")

    def _validate_file_references(self,
                                  tasks: List[Dict[str, Any]],
                                  specification: Dict[str, Any],
                                  file_ids: Set[str]) -> None:
        files_declared = "files" in specification
        referenced_file_ids = set()
        for task in tasks:
            for key in ("inputFiles", "outputFiles"):
                for file_id in task.get(key, []):
                    referenced_file_ids.add(file_id)
                    if files_declared and file_id not in file_ids:
                        self._error(f"File \"{file_id}\", used by task \"{task['id']}\" as one of its {key}, "
                                    f"is not declared in the list of workflow files.")

        if not files_declared and referenced_file_ids:
            self._warning(f"{len(referenced_file_ids)} file(s) are used by tasks while the workflow "
                          f"declares no list of files.")

        unused_file_ids = file_ids - referenced_file_ids
        if unused_file_ids:
            self._warning(f"{len(unused_file_ids)} file(s) are declared but not used by any task.")

    def _compute_level_widths(self, tasks: List[Dict[str, Any]]) -> Optional[List[int]]:
        """
        Compute the width of each level of the workflow DAG, where a task with no parents
        is at level 0 and every other task is at one plus the maximum level of any of its
        parents. Reports cycles as errors.
        """
        task_ids = [task["id"] for task in tasks]
        index_of = {task_id: index for index, task_id in enumerate(task_ids)}
        parent_counts = [0] * len(task_ids)
        children = [[] for _ in task_ids]

        for task in tasks:
            for parent in set(task.get("parents", [])):
                if parent in index_of and parent != task["id"]:
                    parent_counts[index_of[task["id"]]] += 1
                    children[index_of[parent]].append(index_of[task["id"]])

        ready = [index for index, count in enumerate(parent_counts) if count == 0]
        level_widths = []
        number_processed = 0
        while ready:
            level_widths.append(len(ready))
            number_processed += len(ready)
            next_ready = []
            for index in ready:
                for child_index in children[index]:
                    parent_counts[child_index] -= 1
                    if parent_counts[child_index] == 0:
                        next_ready.append(child_index)
            ready = next_ready

        if number_processed != len(task_ids):
            unresolved = sorted(task_ids[index] for index, count in enumerate(parent_counts) if count > 0)
            preview = ", ".join(f"\"{task_id}\"" for task_id in unresolved[:5])
            self._error(f"The workflow graph contains a cycle: {len(unresolved)} task(s) are part of, "
                        f"or depend on, a cycle ({preview}).")
            return None

        return level_widths

    def _validate_execution(self, execution: Dict[str, Any], task_ids: Set[str]) -> None:
        machines_declared = "machines" in execution
        machine_names = self._collect_ids(execution.get("machines", []), "nodeName", "machine")
        executed_task_ids = set()

        for task in execution["tasks"]:
            task_id = task["id"]
            if task_id in executed_task_ids:
                self._error(f"The task \"{task_id}\" is declared more than once in the execution section.")
            executed_task_ids.add(task_id)

            if task_id not in task_ids:
                self._error(f"Task \"{task_id}\" of the execution section is not declared in the list "
                            f"of workflow tasks.")

            for machine_name in task.get("machines", []):
                if machines_declared and machine_name not in machine_names:
                    self._error(f"Machine \"{machine_name}\", used by task \"{task_id}\", is not declared "
                                f"in the list of machines.")
                elif not machines_declared:
                    self._warning(f"Task \"{task_id}\" ran on machine \"{machine_name}\" while the workflow "
                                  f"declares no list of machines.")

        missing_task_ids = task_ids - executed_task_ids
        if missing_task_ids:
            self._warning(f"{len(missing_task_ids)} task(s) have no entry in the execution section.")

        self.logger.debug(f'The workflow has {len(machine_names)} machines.')

    def _validate_specification_metrics(self,
                                        specification: Dict[str, Any],
                                        level_widths: Optional[List[int]]) -> None:
        metrics = specification.get("metrics")
        if not metrics:
            return

        files = specification.get("files", [])
        expected = {
            "numberOfTasks": len(specification["tasks"]),
            "numberOfFiles": len(files),
            "sumOfFileSizesInBytes": sum(f.get("sizeInBytes", 0) for f in files)
        }
        if level_widths:
            expected["numberOfLevels"] = len(level_widths)
            expected["minimumWidth"] = min(level_widths)
            expected["maximumWidth"] = max(level_widths)

        self._compare_metrics(metrics, expected, "specification")

    def _validate_execution_metrics(self, execution: Dict[str, Any]) -> None:
        metrics = execution.get("metrics")
        if not metrics:
            return

        tasks = execution["tasks"]
        expected = {
            "sumTaskRuntimesInSeconds": sum(t["runtimeInSeconds"] for t in tasks if "runtimeInSeconds" in t),
            "totalNumBytesRead": sum(t["readBytes"] for t in tasks if "readBytes" in t),
            "totalNumBytesWritten": sum(t["writtenBytes"] for t in tasks if "writtenBytes" in t)
        }

        self._compare_metrics(metrics, expected, "execution")

    def _compare_metrics(self, metrics: Dict[str, Any], expected: Dict[str, Any], section: str) -> None:
        for key, expected_value in expected.items():
            if key in metrics and not math.isclose(metrics[key], expected_value,
                                                   rel_tol=METRICS_RELATIVE_TOLERANCE):
                self._error(f"The {section} metric \"{key}\" is declared as {metrics[key]}, but the "
                            f"instance data gives {expected_value}.")
