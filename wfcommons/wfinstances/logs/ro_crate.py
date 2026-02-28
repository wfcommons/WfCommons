#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import itertools
import math
import os
import pathlib

from datetime import datetime, timezone
from logging import Logger
from typing import List, Optional

from .abstract_logs_parser import LogsParser
from ...common.file import File
from ...common.machine import Machine
from ...common.task import Task, TaskType
from ...common.workflow import Workflow


class ROCrateLogsParser(LogsParser):
    """
    Parse RO Crate directory to generate workflow instance. This parser has some limitations when it comes to non-file
    dependencies between tasks. It determines these via ParameterConnection type objects in the ro-crate-metadata.json,
    which contain the "instrument" (the cwl file they execute) of the parent and child task. However, since tasks
    can share an "instrument", the parser creates dependencies between every task pair matching the parent and child
    "instrument"s, assuming they're all related.

    :param crate_dir: RO crate directory (contains ro-crate-metadata.json).
    :type crate_dir: pathlib.Path
    :param description: Workflow instance description.
    :type description: Optional[str]
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Optional[Logger]
    """

    def __init__(self,
                 crate_dir: pathlib.Path,
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None,
                 steps_to_ignore: Optional[list[str]]=None) -> None:
        """Create an object of the RO crate parser."""

        # TODO: Decide if these should be RO crate or Streamflow or whatev
        super().__init__('Streamflow-ROCrate', 'https://w3id.org/workflowhub/workflow-ro-crate/1.0', description, logger)

        # Sanity check
        if steps_to_ignore is None:
            steps_to_ignore = []
        if not crate_dir.is_dir():
            raise OSError(f'The provided path does not exist or is not a folder: {crate_dir}')

        metadata: pathlib.Path = crate_dir / 'ro-crate-metadata.json'
        if not metadata.is_file():
            raise OSError(f'Unable to find ro-crate-metadata.json file in: {crate_dir}')
        self.metadata = metadata

        self.crate_dir: pathlib.Path = crate_dir

        self.file_objects = {}

        self.task_id_name_map: dict[str, str] = {}

        self.steps_to_ignore = steps_to_ignore


    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow instance based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: Optional[str]

        :return: A workflow instance object.
        :rtype: Workflow
        """
        self.workflow_name = workflow_name

        # create base workflow instance object
        self.workflow = Workflow(name=self.workflow_name,
                                 description=self.description,
                                 runtime_system_name=self.wms_name,
                                 runtime_system_url=self.wms_url)

        with open(self.metadata, 'r') as f:
            self.data = json.load(f)
        self.graph_data = self.data.get('@graph', [])

        # Dictionary of ro-crate objects by "@id"
        self.lookup = {item["@id"]: item for item in self.graph_data}

        # Find id of the main workflow
        overview = self.lookup.get("./")
        main_workflow_id = overview.get("mainEntity").get("@id")

        create_actions = list(filter((lambda x: x.get('@type') == "CreateAction"), self.graph_data))
        self._create_tasks(create_actions, main_workflow_id)



        return self.workflow

    def _create_tasks(self, create_actions, main_workflow_id):
        # Object to track dependencies between tasks based on files
        files = {}
        # Object to track task's "instrument" for further dependencies
        instruments = {}

        for create_action in create_actions:

            # Handle overall workflow create_action then skip
            if create_action["name"] == f"Run of workflow/{main_workflow_id}":
                self._process_main_workflow(create_action)
                continue

            create_action['name'] = create_action['name'].removeprefix("Run of workflow/")

            # Below would remove the "file.cwl#" tag, which runs the risk
            # of non-uniqueness of action names perhaps
            # create_action['name'] = create_action['name'].split('#', 1)[-1]

            # Check if we should ignore this step
            if create_action["name"] in self.steps_to_ignore:
                continue

            # Get all input & output for the create_action
            input = [obj['@id'] for obj in create_action['object']]
            output = [obj['@id'] for obj in create_action['result']]

            # Filter for actual files
            input_files = self._filter_file_ids(input)
            output_files = self._filter_file_ids(output)


            task = Task(name=create_action['name'],
                        task_id=create_action['name'],
                        # task_id=create_action['name'] + "_" + create_action['@id'],
                        task_type=TaskType.COMPUTE,
                        runtime=self._time_diff(create_action['startTime'], create_action['endTime']),
                        executed_at=create_action['startTime'],
                        input_files=self._get_file_objects(input_files),
                        output_files=self._get_file_objects(output_files),
                        logger=self.logger)
            self.workflow.add_task(task)
            self.task_id_name_map[create_action['@id']] = create_action['name']
            # self.task_id_name_map[create_action['@id']] = create_action['name'] + "_" + create_action['@id']

            # For each file, track which task(s) it is in/output for
            for infile in input_files:
                if infile not in files:
                    files[infile] = {}
                if 'in' not in files[infile]:
                    files[infile]['in'] = []

                files[infile]['in'].append(create_action['@id'])

            for outfile in output_files:
                if outfile not in files:
                    files[outfile] = {}
                if 'out' not in files[outfile]:
                    files[outfile]['out'] = []

                files[outfile]['out'].append(create_action['@id'])

            # For each task, track which 'instrument' it uses
            instrument = create_action['instrument']['@id']
            if instrument not in instruments:
                instruments[instrument] = []
            instruments[instrument].append(create_action['@id'])

        self._add_dependencies(files, instruments)

    def _add_dependencies(self, files, instruments):
        for file in files.values():
            for parent in file.get('out', []):
                for child in file.get('in', []):
                    # self.workflow.add_dependency(parent, child)
                    self.workflow.add_dependency(self.task_id_name_map[parent], self.task_id_name_map[child])

        # Assumes
        parameter_connections = list(filter((lambda x: x.get('@type') == "ParameterConnection"), self.graph_data))
        for parameter_connection in parameter_connections:
            # parameter_connection["sourceParameter"] is either a single dict or a list of dicts,
            # which is bad design but whatever
            source_parameters = parameter_connection["sourceParameter"]
            if not isinstance(source_parameters, list):
                source_parameters = [source_parameters]

            for item in source_parameters:
                source = item["@id"]
                source = source.rsplit("#", 1)[0]   # Trim to get instrument

                target = parameter_connection["targetParameter"]["@id"]
                target = target.rsplit("#", 1)[0]   # Trim to get instrument

                for parent in instruments.get(source, []):
                    for child in instruments.get(target, []):
                        self.workflow.add_dependency(self.task_id_name_map[parent], self.task_id_name_map[child])


    def _time_diff(self, start_time, end_time):
        diff = datetime.fromisoformat(end_time) - datetime.fromisoformat(start_time)
        return diff.total_seconds()

    def _get_file_objects(self, files):
        # Given a list of "@id"s, returns the File objs.
        output = []
        for file in files:
            if file not in self.file_objects:
                self.file_objects[file] = File(file_id=file,
                                               size=os.path.getsize(f"{self.crate_dir}/{file}"),
                                               logger=self.logger)
            output.append(self.file_objects[file])
        return output

    def _filter_file_ids(self, ids):
        # Given a list of "@id"s, returns those with the File type as well as unpacks PropertyValue into Files.
        file_ids = list(filter(lambda x: self.lookup.get(x)['@type'] == 'File', ids))

        property_value_ids = list(filter(lambda x: self.lookup.get(x)['@type'] == 'PropertyValue', ids))
        for property_value_id in property_value_ids:
            property_values = self.lookup.get(property_value_id)['value']

            # Filter out values without "@id"s (i.e. int values, etc.)
            pv_contained_ids = list(filter(lambda x: isinstance(x, dict) and "@id" in x, property_values))
            pv_contained_ids = [obj["@id"] for obj in pv_contained_ids]

            # Recurse to verify everything's a file
            pv_filtered_ids = self._filter_file_ids(pv_contained_ids)

            # Filter duplicates while adding
            file_ids = list(set(file_ids + pv_filtered_ids))

        return file_ids
    def _process_main_workflow(self, main_workflow):
        self.workflow.makespan = self._time_diff(main_workflow['startTime'], main_workflow['endTime'])
        self.workflow.executed_at = main_workflow['startTime']