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
import glob
import pathlib
import csv

from datetime import datetime, timezone, timedelta
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
    :param steps_to_ignore: Names of CWL steps that should be ignored in the translation
    :type steps_to_ignore: Optional[list[str]]
    :param file_extensions_to_ignore: File extensions that should be ignored in the translation
    :type file_extensions_to_ignore: Optional[list[str]]
    :param instruments_to_ignore: Names of instruments that should be ignored in the translation
    :type instruments_to_ignore: Optional[list[str]]
    """

    def __init__(self,
                 crate_dir: pathlib.Path,
                 description: Optional[str] = None,
                 logger: Optional[Logger] = None,
                 steps_to_ignore: Optional[list[str]] = None,
                 file_extensions_to_ignore: Optional[list[str]] = None,
                 instruments_to_ignore: Optional[list[str]] = None,
                 ) -> None:
        """Create an object of the RO crate parser."""

        super().__init__('Streamflow-ROCrate', 'https://w3id.org/workflowhub/workflow-ro-crate/1.0', description, logger)

        # Sanity check
        if not crate_dir.is_dir():
            raise OSError(f'The provided path does not exist or is not a folder: {crate_dir}')
        self.crate_dir: pathlib.Path = crate_dir

        # Find the metadata file
        metadata: pathlib.Path = crate_dir / 'ro-crate-metadata.json'
        if not metadata.is_file():
            raise OSError(f'Unable to find ro-crate-metadata.json file in: {self.crate_dir}')
        self.metadata = metadata

        # Find the Nextflow execution trace if any
        nextflow_execution_trace_files = list(self.crate_dir.glob("results/*/pipeline_info/execution_trace_*.txt"))
        if len(nextflow_execution_trace_files) == 1:
            self.nextflow_execution_trace_file : pathlib.Path = nextflow_execution_trace_files[0]
        else:
            self.nextflow_execution_trace_file = None


        self.file_objects = {}

        self.task_id_name_map: dict[str, str] = {}
        self.data_file_id_name_map: dict[str, str] = {}

        self.steps_to_ignore : list[str] = steps_to_ignore or []
        self.file_extensions_to_ignore : list[str] = file_extensions_to_ignore or []
        self.instruments_to_ignore : list[str] = instruments_to_ignore or []


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

        # Dictionary of application data files
        self._construct_data_file_id_name_map()

        # Task runtime dictionary in case there was a Nextflow trace file
        self.nextflow_trace_times = {}
        if self.nextflow_execution_trace_file:
            self.nextflow_trace_times = self._load_trace(self.nextflow_execution_trace_file)
        
        # File size directory in case this was a Nextflow-generated RO-Crate
        self.nextflow_publish_map = {}
        for item in self.graph_data:
            if item.get('@type') == 'CreateAction' and item['@id'].startswith('#publish/'):
                obj = item['object']
                res = item['result']
                obj_id = obj['@id'] if isinstance(obj, dict) else obj[0]['@id']
                res_id = res['@id'] if isinstance(res, dict) else res[0]['@id']
                self.nextflow_publish_map[obj_id] = res_id

        # Find id of the main workflow
        overview = self.lookup.get("./")
        main_workflow_id = overview.get("mainEntity").get("@id")

        create_actions = list(filter((lambda x: x.get('@type') == "CreateAction"), self.graph_data))
        self._create_tasks(create_actions, main_workflow_id)

        return self.workflow


    def _construct_data_file_id_name_map(self):
        for item in self.graph_data:
            if item["@type"] != "File":
                continue
            id = item["@id"]
            #if "alternateName" not in item:
            #    continue
            #alternate_name = item["alternateName"]
            alternate_name = item.get("alternateName", id)
            self.data_file_id_name_map[id] = alternate_name


    def _create_tasks(self, create_actions, main_workflow_id):
        # Object to track dependencies between tasks based on files
        files = {}
        # Object to track task's "instrument" for further dependencies
        instruments = {}

        for create_action in create_actions:
            # Handle overall workflow create_action then skip
            if create_action["name"] == f"Run of workflow/{main_workflow_id}":
                # Streamflow-generated RO-Crate
                self._process_main_workflow(create_action)
                continue
            elif "Nextflow workflow run" in create_action["name"]:
                # Nextflow-generated RO-Crate
                self._process_main_workflow(create_action)
                continue

            # Ignore all (Nextflow-generated) publish tasks if any
            if create_action['@id'].startswith('#publish/'):
                continue

            create_action['name'] = create_action['name'].removeprefix("Run of workflow/")
            # print("***************************************")
            # print("DEALING WITH TASK:", create_action['name'])

            # Below would remove the "file.cwl#" tag, which runs the risk
            # of non-uniqueness of action names perhaps
            # create_action['name'] = create_action['name'].split('#', 1)[-1]

            # Check if we should ignore this step
            if create_action["name"] in self.steps_to_ignore:
                continue

            # Get all input & output for the create_action
            input = [obj['@id'] if isinstance(obj, dict) else obj for obj in create_action['object']]
            output = [obj['@id'] if isinstance(obj, dict) else obj for obj in create_action['result']]

            # Filter for actual files
            input_files = self._filter_file_ids(input)
            output_files = self._filter_file_ids(output)

            # Figure out start/end times for the task
            start_time = create_action.get('startTime')
            end_time = create_action.get('endTime')
            if not start_time or not end_time:
                start_time, end_time = self.nextflow_trace_times.get(create_action['name'], (None, None))

            task = Task(name=create_action['name'],
                        # task_id=create_action['name'],
                        task_id=create_action['name'] + "_" + create_action['@id'],
                        task_type=TaskType.COMPUTE,
                        runtime=self._time_diff(start_time, end_time),
                        executed_at=create_action.get('startTime',''),
                        input_files=self._get_file_objects(input_files),
                        output_files=self._get_file_objects(output_files),
                        logger=self.logger)

            self.workflow.add_task(task)
            # self.task_id_name_map[create_action['@id']] = create_action['name']
            self.task_id_name_map[create_action['@id']] = create_action['name'] + "_" + create_action['@id']

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
            if create_action.get('instrument'):
                instrument = create_action['instrument']['@id']
                if instrument not in instruments:
                    instruments[instrument] = []
                instruments[instrument].append(create_action['@id'])

        self._add_dependencies(files, instruments)

    def _add_dependencies(self, files, instruments):

        # File dependencies
        for file in files.values():
            for parent in file.get('out', []):
                for child in file.get('in', []):
                    self.workflow.add_dependency(self.task_id_name_map[parent], self.task_id_name_map[child])

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

                if source in self.instruments_to_ignore or target in self.instruments_to_ignore:
                    continue
                # print("source", source, "----> target", target)

                for parent in instruments.get(source, []):
                    for child in instruments.get(target, []):
                        self.workflow.add_dependency(self.task_id_name_map[parent], self.task_id_name_map[child])

    def _time_diff(self, start_time, end_time):
        if not start_time or not end_time:
            return 0.0
        diff = datetime.fromisoformat(end_time) - datetime.fromisoformat(start_time)
        return diff.total_seconds()

    def _get_file_objects(self, files):
        # Given a list of "@id"s, returns the File objs.
        output = []
        for file in files:
            if file not in self.file_objects:
                #self.file_objects[file] = File(file_id=self.data_file_id_name_map[file],
                #                               size=os.path.getsize(f"{self.crate_dir}/{file}"),
                #                               logger=self.logger)
                if file not in self.data_file_id_name_map:
                    # File is referenced but not in the map — use its @id as the name
                    self.logger.warning(f"File not in data_file_id_name_map, using @id as name: {file}") if self.logger else None
                    file_name = file
                else:
                    file_name = self.data_file_id_name_map[file]

                # Figure out the file size

                # Straight from the RO-Crate?
                try:
                    obj = self.lookup.get(file, {})
                    size = int(obj.get('contentSize', 0))
                except (OSError, ValueError):
                    size = 0
                    
                # From the Nextflow-specific "#publish" tasks?
                if not size:
                    resolved = self.nextflow_publish_map.get(file, file)
                    try:
                        size = os.path.getsize(self.crate_dir / resolved)
                    except (OSError, ValueError):
                        work_path = self._resolve_task_file_path(file)
                        size = os.path.getsize(work_path) if work_path else 0

                # Perhaps a (Nextflow-specific) absolute path?
                if not size:
                    if file.startswith('file://'):
                        abs_path = pathlib.Path(file[len('file://'):])
                        try:
                            size = os.path.getsize(abs_path)
                        except (OSError, ValueError):
                            size = 0

                # Create the file object
                self.file_objects[file] = File(file_id=file_name,
                                            size=size,
                                            logger=self.logger)
            output.append(self.file_objects[file])
        return output

    def _filter_file_ids(self, ids):

        file_ids = list(filter(lambda x: (self.lookup.get(x) or {}).get('@type') in ('File','CreativeWork'), ids))
        # Ignore the files that start with http:// or https://
        file_ids = [x for x in file_ids if not x.startswith("http://")]
        file_ids = [x for x in file_ids if not x.startswith("https://")]
        property_value_ids = list(filter(lambda x: (self.lookup.get(x) or {}).get('@type') == 'PropertyValue', ids))
        for property_value_id in property_value_ids:
            property_values = self.lookup.get(property_value_id)['value']
            # If the lookup fails, ignore
            if not property_values:
                    continue
            if not isinstance(property_values, list):
                property_values = [property_values]

            # Filter out values without "@id"s (i.e. int values, etc.)
            pv_contained_ids = list(filter(lambda x: isinstance(x, dict) and "@id" in x, property_values))
            pv_contained_ids = [obj["@id"] for obj in pv_contained_ids]

            # Recurse to verify everything's a file
            pv_filtered_ids = self._filter_file_ids(pv_contained_ids)

            # Filter duplicates while adding
            file_ids = list(set(file_ids + pv_filtered_ids))

        # Removing files based on file extensions
        to_return = []
        for file_id in file_ids:
            to_ignore = False
            for suffix in self.file_extensions_to_ignore:
                if self.data_file_id_name_map[file_id].endswith(suffix):
                    to_ignore = True
                    break
            if not to_ignore:
                to_return.append(file_id)

        return to_return

    import csv
    from datetime import datetime, timedelta

    def _load_trace(self, trace_file: pathlib.Path) -> dict:
        """
        Parse a Nextflow execution trace file.
        Returns dict of task_name -> (start_iso, end_iso).
        """
        times = {}
        with open(trace_file, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                name = row['name'].strip()
                submit_str = row['submit'].strip()
                duration_str = row['duration'].strip()
                if not submit_str or not duration_str:
                    continue
                try:
                    start = datetime.strptime(submit_str, '%Y-%m-%d %H:%M:%S.%f')
                    duration_sec = self._parse_duration(duration_str)
                    end = start + timedelta(seconds=duration_sec)
                    # Store as ISO strings to match _time_diff expectations
                    times[name] = (start.isoformat(), end.isoformat())
                except Exception:
                    continue
        return times

    def _parse_duration(self, duration_str: str) -> float:
        """Parse Nextflow duration strings like '2m 40s', '519ms', '1h 3m 2s'."""
        import re
        total = 0.0
        for value, unit in re.findall(r'(\d+(?:\.\d+)?)\s*(h|m|s|ms)', duration_str):
            value = float(value)
            if unit == 'h':
                total += value * 3600
            elif unit == 'm':
                total += value * 60
            elif unit == 's':
                total += value
            elif unit == 'ms':
                total += value / 1000
        return total
    
    def _resolve_task_file_path(self, file_id: str) -> Optional[pathlib.Path]:
        """Resolve a #task/<hash>/<filename> id to its work/ directory path."""
        if not file_id.startswith('#task/'):
            return None
        # Strip '#task/'
        rest = file_id[len('#task/'):]
        # rest is now '<hash>/<filename>'
        hash_and_name = rest.split('/', 1)
        if len(hash_and_name) != 2:
            return None
        task_hash, filename = hash_and_name
        # Nextflow work dir structure: work/<first2chars>/<rest>/
        work_path = self.crate_dir / 'work' / task_hash[:2] / task_hash[2:] / filename
        if work_path.exists():
            return work_path
        return None

    def _process_main_workflow(self, main_workflow):
        self.workflow.makespan = self._time_diff(main_workflow['startTime'], main_workflow['endTime'])
        self.workflow.executed_at = main_workflow['startTime']
