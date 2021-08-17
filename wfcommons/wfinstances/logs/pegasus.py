#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import fnmatch
import json
import os
import yaml
import xml.etree.ElementTree

from datetime import datetime
from logging import Logger
from typing import Optional

from .abstract_logs_parser import LogsParser
from ...common.file import File, FileLink
from ...common.machine import Machine, MachineSystem
from ...common.task import Task, TaskType
from ...common.workflow import Workflow


class PegasusLogsParser(LogsParser):
    """
    Parse Pegasus submit directory to generate workflow instance.

    :param submit_dir: Pegasus submit directory.
    :type submit_dir: str
    :param legacy: Whether the submit directory is from a Pegasus 4.x version.
    :type legacy: bool
    :param description: Workflow instance description.
    :type description: str
    :param ignore_auxiliary: Ignore auxiliary jobs.
    :type ignore_auxiliary: bool
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self,
                 submit_dir: str,
                 description: Optional[str] = None,
                 ignore_auxiliary: Optional[bool] = True,
                 legacy: Optional[bool] = False,
                 logger: Optional[Logger] = None) -> None:
        """Create an object of the pegasus log parser."""
        super().__init__('Pegasus', 'https://pegasus.isi.edu', description, logger)

        # Sanity check
        if not os.path.isdir(submit_dir):
            raise OSError('The provided path does not exist or is not a directory: {}'.format(submit_dir))

        if ignore_auxiliary:
            self.logger.warning('Ignoring Pegasus auxiliary jobs.')

        self.submit_dir = submit_dir
        self.legacy = legacy
        self.ignore_auxiliary = ignore_auxiliary
        self.files_map = {}

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Create workflow instance based on the workflow execution logs.

        :param workflow_name: The workflow name.
        :type workflow_name: str

        :return: A workflow instance object.
        :rtype: Workflow
        """
        self.workflow_name = workflow_name

        if self.legacy:
            # parse braindump file
            self._parse_braindump()
            # parse DAX file
            self._parse_dax()

        else:
            # parse workflow YAML
            self._parse_workflow()

        # parse DAG file
        self._parse_dag()

        return self.workflow

    def _parse_braindump(self):
        """Parse the Pegasus braindump.txt file"""
        braindump_file = '{}/braindump.txt'.format(self.submit_dir)

        if not os.path.exists(braindump_file):
            raise OSError('Unable to find braindump file: {}'.format(braindump_file))

        with open(braindump_file) as f:
            for line in f:
                if line.startswith('planner_version'):
                    wms_version = line.split()[1]
                elif line.startswith('pegasus_wf_name'):
                    self.instance_name = line.split()[1]
                elif line.startswith('timestamp'):
                    executed_at = line.split()[1]

        # sanity checks
        if not wms_version:
            self.logger.warning('Unable to determine pegasus version.')
        if not self.instance_name:
            self.logger.warning('Unable to determine instance name from "pegasus_wf_name".')
        if not executed_at:
            self.logger.warning('Unable to determine execution time from "timestamp".')

        # create base workflow instance object
        self.workflow = Workflow(name=self.workflow_name,
                                 description=self.description,
                                 wms_name=self.wms_name,
                                 wms_version=wms_version,
                                 wms_url=self.wms_url,
                                 executed_at=executed_at)

    def _parse_workflow(self):
        """Parse the Workflow file."""
        workflows_file_list = self._fetch_all_files('yml', '*workflow')

        if len(workflows_file_list) < 1:
            raise OSError('Unable to find workflow file.')

        workflow_file = workflows_file_list[0]

        with open(workflow_file) as f:
            data = yaml.load(f, Loader=yaml.SafeLoader)
            self.logger.info('Processing Pegasus workflow file: {}'.format(os.path.basename(workflow_file)))

            # create base workflow instance object
            self.workflow = Workflow(name=self.workflow_name,
                                     description=self.description,
                                     wms_name=self.wms_name,
                                     wms_version=data['pegasus'],
                                     wms_url=self.wms_url,
                                     executed_at=data['x-pegasus']['createdOn'])

            for j in data['jobs']:
                task_name = '{}_{}'.format(j['name'], j['id'])

                list_files = [File(
                    name=f['lfn'],
                    size=0,
                    link=FileLink(f['type']),
                    logger=self.logger
                ) for f in j['uses']]

                self.workflow.add_node(
                    task_name,
                    task=Task(
                        name=task_name,
                        task_id=j['id'],
                        category=j['name'],
                        task_type=TaskType.COMPUTE,
                        runtime=0,
                        args=j['arguments'],
                        cores=0,
                        files=list_files,
                        logger=self.logger
                    )
                )

    def _parse_dax(self):
        """Parse the DAX file."""
        dax_list = self._fetch_all_files("dax", "")
        if len(dax_list) < 1:
            dax_list = self._fetch_all_files("xml", "")
            if len(dax_list) < 1:
                raise OSError('The directory contains no ".dax" or ".xml" file')

        dax_file = dax_list[0]
        self.logger.info('Processing Pegasus DAX file: {}'.format(os.path.basename(dax_file)))

        line_num = 0
        temp_file = None
        with open(dax_file) as f:
            for line in f:
                line_num += 1
                if line.startswith('<?xml'):
                    if line_num == 1:
                        break
                    else:
                        temp_file = open('.pegasus-parser-tmp', 'w')

                if line.startswith('</adag>'):
                    temp_file.write(line)
                    dax_file = temp_file.name
                    temp_file.close()
                    break

                if temp_file:
                    temp_file.write(line)

        try:
            e = xml.etree.ElementTree.parse(dax_file).getroot()
            for j in e.findall('{http://pegasus.isi.edu/schema/DAX}job'):
                task_name = str(j.get('name')) + '_' + str(j.get('id'))

                list_files = [File(
                    name=f.get('name') if not f.get('name') is None else f.get('file'),
                    size=0,
                    link=FileLink(f.get('link')),
                    logger=self.logger
                ) for f in j.findall('{http://pegasus.isi.edu/schema/DAX}uses')]

                self.workflow.add_node(
                    task_name,
                    task=Task(
                        name=task_name,
                        task_id=str(j.get('id')),
                        category=str(j.get('name')),
                        task_type=TaskType.COMPUTE,
                        runtime=0,
                        args=[],
                        cores=0,
                        files=list_files,
                        logger=self.logger
                    )
                )

        except xml.etree.ElementTree.ParseError as ex:
            self.logger.warning(str(ex))

        # removing temporary file
        if temp_file:
            os.remove(dax_file)

    def _fetch_all_files(self, extension: str, file_name: str = ""):
        """
        Fetch all files from the directory and its hierarchy

        :param extension: file extension to be searched for
        :type extension: str
        :param file_name: file_name to be searched
        :type file_name: str

        :return: List of file names that match
        :rtype: List[str]
        """
        files = []
        for root, dirnames, filenames in os.walk(self.submit_dir):
            if len(file_name) == 0:
                for filename in fnmatch.filter(filenames, '*.{}'.format(extension)):
                    files.append(os.path.join(root, filename))
            else:
                for filename in fnmatch.filter(filenames, '{}.{}'.format(file_name, extension)):
                    files.append(os.path.join(root, filename))
        return files

    def _parse_dag(self):
        """Parse the DAG file."""
        dags_list = self._fetch_all_files("dag", "")
        if len(dags_list) < 1:
            raise OSError('The directory contains no ".dag" file')

        dag_file = dags_list[0]
        self.logger.info('Processing Pegasus DAG file: {}'.format(os.path.basename(dag_file)))

        num_tasks = 0
        tasks_set = set()

        with open(dag_file) as f:
            for line in f:
                if line.startswith('JOB'):
                    num_tasks += 1
                    task_name = line.split()[1]

                    # find the task
                    task = None
                    for node in self.workflow.nodes.data():
                        if node[0].lower() == task_name.lower():
                            task = node[1]['task']
                            break

                    if not task and not self.ignore_auxiliary:
                        task = Task(
                            name=task_name,
                            task_type=TaskType.COMPUTE,
                            runtime=0,
                            args=[],
                            cores=0,
                            files=[],
                            logger=self.logger
                        )
                        self.workflow.add_node(task_name, task)
                    else:
                        self._parse_meta_file(task_name)

                    # Parsing job stdout file
                    if task:
                        tasks_set.add(task.name)
                        self._parse_job_output(task)

                elif line.startswith('PARENT'):
                    # Typically, parent/child references are at the end of the DAG file
                    s = line.split()
                    parent = s[1]
                    child = s[3]
                    for node in self.workflow.nodes.data():
                        if node[0].lower() == child.lower() and parent in tasks_set:
                            task = node[1]['task']
                            self.workflow.add_edge(parent, task.name, weight=0)
                            break

        # parse workflow files
        for node in self.workflow.nodes.data():
            task = node[1]['task']
            for f in task.files:
                if f.name in self.files_map:
                    f.size = int(self.files_map[f.name])

        # parse workflow makespan
        dagman_file = dag_file + '.dagman.out'
        self.logger.debug('Processing Pegasus DAGMan output file.')
        with open(dagman_file) as f:
            lines = f.readlines()
            try:
                s = datetime.strptime(' '.join(lines[0].split()[0:2]), '%m/%d %H:%M:%S')
                e = datetime.strptime(' '.join(lines[-1].split()[0:2]), '%m/%d %H:%M:%S')
            except ValueError:
                s = datetime.strptime(' '.join(lines[0].split()[0:2]), '%m/%d/%y %H:%M:%S')
                e = datetime.strptime(' '.join(lines[-1].split()[0:2]), '%m/%d/%y %H:%M:%S')
            self.workflow.makespan = (e - s).total_seconds()

        self.logger.debug('Found {} jobs.'.format(num_tasks))

    def _parse_meta_file(self, task_name):
        """
        Parse the Pegasus meta file (generated from pegasus-integrity)

        :param task_name: Task file name.
        :type task_name: str
        """
        meta_list = self._fetch_all_files("meta", task_name)
        if not meta_list:
            self.logger.warning('Job {} has no meta record (skipping meta analysis).'.format(task_name))
            return

        with open(meta_list[0]) as metadata:
            m = json.load(metadata)
            for f in m:
                if f['_id'] not in self.files_map:
                    self.files_map[f['_id']] = f['_attributes']['size']

    def _parse_job_output(self, task):
        """
        Parse the kickstart job output file (e.g., .out.000).

        :param task: Task object.
        :type task: Task
        """
        output_list = self._fetch_all_files('out.*', task.name)
        if len(output_list) == 0:
            self.logger.warning('Job has no kickstart record. Skipping it.')
            if task.name.lower().startswith(('stage_', 'create_dir', 'cleanup', 'clean_up', 'register_')):
                task.type = TaskType.AUXILIARY
            return

        if len(output_list) > 1:
            self.logger.debug('Job "{}" has multiple runs. Parsing last attempt.'.format(task.name))
        output_file = output_list[-1]

        # setting task type if transfer
        if task.name.lower().startswith(('stage_in_', 'stage_out')):
            task.type = TaskType.TRANSFER

        # parsing job output file
        self.logger.debug('Parsing Job output file: {}'.format(output_file))
        if self.legacy:
            self._parse_job_output_legacy(task, output_file)
        else:
            self._parse_job_output_latest(task, output_file)

        # parsing meta file
        self._parse_meta_file(task.name)

        # parsing .sub file to get job priorities
        sub_list = self._fetch_all_files("sub", task.name)
        if not sub_list:
            self.logger.warning('Job {} has no .sub record. Skipping it.'.format(task.name))
        else:
            with open(sub_list[0]) as f:
                for line in f:
                    if line.startswith('priority'):
                        task.priority = int(line.split()[2])

    def _parse_job_output_latest(self, task, output_file):
        """
        Parse the kickstart job output file in YAML format (e.g., .out.000).

        :param task: Task object.
        :type task: Task
        :param output_file: Output file name.
        :type output_file: str
        """
        tmp_file = '.pegasus-parser-tmp'
        with open(tmp_file, 'w') as t:
            with open(output_file) as f:
                for line in f:
                    if not line.startswith('---------------'):
                        t.write(line)

        with open(tmp_file) as f:
            data = yaml.load(f, Loader=yaml.FullLoader)[0]

            task.program = data['transformation']

            if data['transformation'].startswith('pegasus:') or task.name.lower().startswith('chmod_'):
                task.type = TaskType.AUXILIARY

            mainjob = data['mainjob']
            task.runtime = float(mainjob['duration'])
            task.memory = int(mainjob['usage']['maxrss'])
            total_time = float(mainjob['usage']['utime']) + float(mainjob['usage']['stime'])
            if total_time > 0:
                task.avg_cpu = float('%.4f' % (100 * (total_time / task.runtime)))

            bytes_read = 0
            bytes_written = 0

            # get job memory and I/O information
            if mainjob['procs']:
                for p in mainjob['procs']:
                    bytes_read += max(int(p['rbytes']), int(p['rchar']))
                    bytes_written += max(int(p['wbytes']), int(p['wchar']))
                if bytes_read > 0:
                    task.bytes_read = bytes_read
                if bytes_written > 0:
                    task.bytes_written = bytes_written

            # machine
            task.machine = Machine(
                name=data['machine']['uname_nodename'],
                cpu={
                    'count': data['machine']['cpu_count'],
                    'speed': data['machine']['cpu_speed'],
                    'vendor': data['machine']['cpu_vendor']
                },
                system=MachineSystem(data['machine']['uname_system']),
                architecture=data['machine']['uname_machine'],
                memory=data['machine']['ram_total'],
                release=data['machine']['uname_release']
            )

        os.remove(tmp_file)

    def _parse_job_output_legacy(self, task, output_file):
        """
        Parse the kickstart job output file in XML format (e.g., .out.000).

        :param task: Task object.
        :type task: Task
        :param output_file: Output file name.
        :type output_file: str
        """
        runtime = 0
        total_time = 0
        bytes_read = 0
        bytes_written = 0
        memory = 0
        args = []

        # clean output file from PBS logs
        line_num = 0
        temp_file = None

        with open(output_file) as f:
            for line in f:
                line_num += 1
                if line.startswith('<?xml'):
                    if line_num == 1:
                        break
                    else:
                        temp_file = open('.pegasus-parser-tmp', 'w')

                if line.startswith('</invocation>'):
                    temp_file.write(line)
                    output_file = temp_file.name
                    temp_file.close()
                    break

                if temp_file:
                    temp_file.write(line)

        try:
            e = xml.etree.ElementTree.parse(output_file).getroot()
            # main job information
            task.program = e.get('transformation')
            if e.get('transformation').startswith('pegasus:') or task.name.lower().startswith('chmod_'):
                task.type = TaskType.AUXILIARY

            for mj in e.findall('{http://pegasus.isi.edu/schema/invocation}mainjob'):
                runtime += float(mj.get('duration'))

                # get average cpu utilization
                for u in mj.findall('{http://pegasus.isi.edu/schema/invocation}usage'):
                    total_time += float(u.get('utime')) + float(u.get('stime'))

                # get job arguments
                for av in mj.findall('{http://pegasus.isi.edu/schema/invocation}argument-vector'):
                    for a in av.findall('{http://pegasus.isi.edu/schema/invocation}arg'):
                        args.append(a.text)

                # get job memory information
                for p in mj.findall('{http://pegasus.isi.edu/schema/invocation}proc'):
                    memory += float(p.get('rsspeak'))
                    bytes_read += max(int(p.get('rbytes')), int(p.get('rchar')))
                    bytes_written += max(int(p.get('wbytes')), int(p.get('wchar')))

                # machine
                machine = {}
                for m in e.findall('{http://pegasus.isi.edu/schema/invocation}machine'):
                    for u in m.findall('{http://pegasus.isi.edu/schema/invocation}uname'):
                        machine['system'] = MachineSystem(u.get('system'))
                        machine['architecture'] = u.get('machine')
                        machine['release'] = u.get('release')
                        machine['nodeName'] = u.get('nodename')
                    for u in m.findall('{http://pegasus.isi.edu/schema/invocation}linux'):
                        for r in u.findall('{http://pegasus.isi.edu/schema/invocation}ram'):
                            machine['memory'] = int(r.get('total'))
                        for c in u.findall('{http://pegasus.isi.edu/schema/invocation}cpu'):
                            machine['cpu'] = {
                                'count': int(c.get('count')),
                                'speed': int(c.get('speed')),
                                'vendor': c.get('vendor')
                            }
                    task.machine = Machine(
                        name=machine['nodeName'],
                        cpu={
                            'count': machine['cpu']['count'],
                            'speed': machine['cpu']['speed'],
                            'vendor': machine['cpu']['vendor']
                        },
                        system=machine['system'],
                        architecture=machine['architecture'],
                        memory=machine['memory'],
                        release=machine['release']
                    )

            task.runtime = runtime
            if total_time > 0:
                task.avg_cpu = float('%.4f' % (100 * (total_time / runtime)))
            if memory > 0:
                task.memory = memory
            if bytes_read > 0:
                task.bytes_read = bytes_read
            if bytes_written > 0:
                task.bytes_written = bytes_written
            if len(args) > 0:
                task.args = args

        except xml.etree.ElementTree.ParseError as ex:
            # parse create_dir output file
            if task.name.lower().startswith(('stage_', 'create_dir')):
                task.type = TaskType.AUXILIARY
                with open(output_file) as f:
                    s = None
                    e = None
                    for line in f:
                        values = line.split()
                        if not s:
                            s = datetime.strptime('{} {}'.format(values[0], values[1]), '%Y-%m-%d %H:%M:%S,%f')
                        e = datetime.strptime('{} {}'.format(values[0], values[1]), '%Y-%m-%d %H:%M:%S,%f')
                    task.runtime = (e - s).total_seconds()

            else:
                self.logger.warning('{}: {}'.format(task.name, str(ex)))

        # cleaning temporary file
        if temp_file:
            os.remove(output_file)
