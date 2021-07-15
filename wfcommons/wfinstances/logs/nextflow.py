import fnmatch
import json
import os
import itertools
import glob
import csv
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

class NextflowLogsParser(LogsParser):
    """ Parse Nextflow submit directory to generate workflow trace."""
    def __init__(self,
                 submit_dir: str, #
                 description: Optional[str] = None,
                 logger:Optional[Logger] = None) -> None:
        super.__init__('Nextflow','https://www.nextflow.io', description, logger)

        if not os.path.isdir(submit_dir):
            raise OSError('The provided path does not exist or is not a directory: {}'.format(submit_dir))
        self.submit_dir = submit_dir
        self.files_map = {}
        self.text_files = None
        self.line_count = None

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
               Create workflow trace based on the workflow execution logs.

               :param workflow_name: The workflow name.
               :type workflow_name: str

               :return: A workflow trace object.
               :rtype: Workflow
        """
        self.workflow_name = workflow_name
        self.workflow = Workflow(name=self.workflow_name,
                                 description=self.description,
                                 executed_at=self.executed_at,
                                 makespan=self.makespan
                                 )

        self._parse_workflow_file()

        # parse makeflow log file
        self._parse_makeflow_log_file()

        # parse resource monitor files
        self._parse_resource_monitor_logs()

        return self.workflow

    def parseFiles(self):
        self.text_files = glob.glob(self.submit_dir + "/execution_trace.txt", recursive=True)
        # text_files is now the file "execution_trace.txt
        file = open(self.text_files, "r")
        self.line_count = 0
        for line in file:
            if line != "\n":
                self.line_count += 1
        file.close()
        data = [12][self.line_count]
        table = open(self.text_files,"r")
        #parses the table format as a 2D array
        for i in len(15):
            for j in len(self.line_count):
                temp = table.readline()
                string_list = temp.split(" ")
                for x in string_list:
                    data[i][j] = string_list[x]
        for i in len(data):
            temp = {
                "name": data[i][3],
                "type": "compute",
                "runtime": None,
                "parents": [],
                "files": [],
                "avgCPU": data[i][10],
                "arguments": [],
                "bytesRead": data[i][13],
                "bytesWritten": data[i][14],
                "memory": data[i][12]
            }
            self.files_map[i] = temp

    def _parse_workflow_file(self):
        """Parse the nextlflow workflow file and build the workflow structure."""
        task_id_counter = 1

        with open(self.mf_file) as f:
            outputs = []
            inputs = []
            for line in f:
                if ':' in line:
                    outputs = line.split(':')[0].split()
                    inputs = line.split(':')[1].split()

                    for file in itertools.chain(outputs, inputs):
                        if not file in self.files_map:
                            self.files_map[file] = {'task_name': None, 'children': [], 'file': []}

                elif len(line.strip()) > 0:
                    # task execution command
                    prefix = line.replace('./', '').replace('perl', '').strip().split()[1 if 'LOCAL' in line else 0]
                    task_name = "{}_ID{:06d}".format(prefix, task_id_counter)
                    task_id_counter += 1

                    # create list of task files
                    list_files = []
                    list_files.extend(self._create_files(outputs, FileLink.OUTPUT, task_name))
                    list_files.extend(self._create_files(inputs, FileLink.INPUT, task_name))

                    # create task
                    args = ' '.join(line.replace('LOCAL', '').replace('perl', '').strip().split())
                    task = Task(name=task_name,
                                task_type=TaskType.COMPUTE,
                                runtime=0,
                                args=args.split(),
                                cores=1,
                                files=list_files,
                                logger=self.logger)
                    self.workflow.add_node(task_name, task=task)
                    self.args_map[args] = task

        # adding edges
        for file in self.files_map:
            for child in self.files_map[file]['children']:
                if self.files_map[file]['task_name']:
                    self.workflow.add_edge(self.files_map[file]['task_name'], child)

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

