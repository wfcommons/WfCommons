#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import uuid

from abc import ABC, abstractmethod
from os import path
from typing import Any, Dict, List, Optional

from ...common.file import File, FileLink
from ...common.job import Job, JobType
from ...common.workflow import Workflow
from ...utils import generate_rvs


class WorkflowRecipe(ABC):
    def __init__(self, name: str, data_footprint: Optional[int], num_jobs: Optional[int]) -> None:
        """
        :param name:
        :param data_footprint:
        :param num_jobs:
        """
        self.name: str = name
        self.data_footprint = data_footprint
        self.num_jobs = num_jobs
        self.workflows: List[Workflow] = []
        self.jobs_files: Dict[str, List[File]] = {}
        self.job_id_counter = 1

    @abstractmethod
    def _workflow_recipe(self) -> Dict:
        """
        """

    @classmethod
    @abstractmethod
    def from_num_jobs(self, num_jobs: int) -> 'WorkflowRecipe':
        """
        :param num_jobs: The upper bound for the total number of jobs in the worklfow.
        :type num_jobs: int
        """

    @abstractmethod
    def build_workflow(self, workflow_name: str = None) -> Workflow:
        """
        :param workflow_name:
        """

    def _generate_job(self, job_name: str, job_id: str, input_files: List[File] = None,
                      files_recipe: Dict[FileLink, Dict[str, int]] = None) -> Job:
        """
        :param job_name:
        :param job_id:
        :param input_files:
        :param files_recipe:
        """
        job_recipe = self._workflow_recipe()[job_name]

        # runtime
        runtime: float = float(format(generate_rvs(job_recipe['runtime']['distribution'],
                                                   job_recipe['runtime']['min'],
                                                   job_recipe['runtime']['max']), '.3f'))

        # linking previous generated output files as input files
        self.jobs_files[job_id] = []
        if input_files:
            for f in input_files:
                if f.link == FileLink.OUTPUT:
                    self.jobs_files[job_id].append(File(name=f.name, size=f.size, link=FileLink.INPUT))

        # generate additional in/output files
        self._generate_files(job_id, job_recipe['input'], FileLink.INPUT, files_recipe)
        self._generate_files(job_id, job_recipe['output'], FileLink.OUTPUT, files_recipe)

        return Job(
            name=job_id,
            job_type=JobType.COMPUTE,
            runtime=runtime,
            machine=None,
            args=[],
            cores=1,
            avg_cpu=None,
            bytes_read=None,
            bytes_written=None,
            memory=None,
            energy=None,
            avg_power=None,
            priority=None,
            files=self.jobs_files[job_id]
        )

    def _generate_job_name(self, prefix: str) -> str:
        """
        :param prefix:
        """
        job_name = "{}_{:08d}".format(prefix, self.job_id_counter)
        self.job_id_counter += 1
        return job_name

    def _generate_files(self, job_id: str, recipe: Dict[str, Any], link: FileLink,
                        files_recipe: Dict[FileLink, Dict[str, int]] = None):
        """
        :param job_id:
        :param recipe:
        :param link:
        :param files_recipe:
        """
        extension_list: List[str] = []
        for f in self.jobs_files[job_id]:
            if f.link == link:
                extension_list.append(path.splitext(f.name)[1] if '.' in f.name else f.name)

        for extension in recipe:
            if extension not in extension_list:
                num_files = 1
                if files_recipe and link in files_recipe and extension in files_recipe[link]:
                    num_files = files_recipe[link][extension]
                for _ in range(0, num_files):
                    self.jobs_files[job_id].append(self._generate_file(extension, recipe, link))

    def _generate_file(self, extension: str, recipe: Dict[str, Any], link: FileLink):
        """
        :param extension:
        :param recipe:
        :param link:
        """
        return File(name=str(uuid.uuid4()) + extension,
                    link=link,
                    size=int(generate_rvs(recipe[extension]['distribution'],
                                          recipe[extension]['min'],
                                          recipe[extension]['max'])))

    def _get_files_by_job_and_link(self, job_id: str, link: FileLink) -> List[File]:
        """
        :param job_id:
        :param link:
        """
        files_list: List[File] = []
        for f in self.jobs_files[job_id]:
            if f.link == link:
                files_list.append(f)
        return files_list
