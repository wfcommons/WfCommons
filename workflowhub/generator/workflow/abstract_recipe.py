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
from typing import Any, Dict, List, Optional

from ...common.file import File, FileLink
from ...common.job import Job
from ...common.workflow import Workflow
from ...utils import generate_rvs


class WorkflowRecipe(ABC):
    def __init__(self, name: str = None) -> None:
        self.name: str = name
        self.workflows: List[Workflow] = []
        self.jobs_files: Dict[str, List[File]] = {}

    @abstractmethod
    def _workflow_recipe(self) -> Dict:
        pass

    @abstractmethod
    def build_workflow(self) -> Workflow:
        pass

    def _generate_job(self, job_name: str, job_id: str, input_files: Optional[List[File]]) -> Job:
        """
        :param job_name:
        :param job_id:
        """
        job_recipe = self._workflow_recipe()[job_name]

        # runtime
        runtime: float = generate_rvs(job_recipe['runtime']['distribution'],
                                      job_recipe['runtime']['min'],
                                      job_recipe['runtime']['max'])

        # files
        self.jobs_files[job_id] = []
        if input_files:
            for f in input_files:
                if f.link != FileLink.INPUT:
                    self.jobs_files[job_id].append(File(name=f.name, size=f.size, link=FileLink.INPUT))
        else:
            self._generate_files(job_id, job_recipe['input'], FileLink.INPUT)
        self._generate_files(job_id, job_recipe['output'], FileLink.OUTPUT)

        return Job(
            name=job_id,
            job_type='compute',
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

    def _generate_files(self, job_id: str, recipe: Dict[str, Any], link: FileLink):
        """
        :param job_id:
        :param recipe:
        :param link:
        """
        for extension in recipe:
            self.jobs_files[job_id].append(self._generate_file(extension, recipe, link))

    def _generate_file(self, extension: str, recipe: Dict[str, Any], link: FileLink):
        """
        :param extension:
        :param recipe:
        :param link:
        """
        file_name = str(uuid.uuid4()) + extension
        file_size = int(generate_rvs(recipe[extension]['distribution'],
                                     recipe[extension]['min'],
                                     recipe[extension]['max']))
        print("FILE: {} ({}) - {}".format(file_name, link, str(file_size)))
        return File(name=file_name, size=file_size, link=link)

    def _get_files_by_job_and_link(self, job_name: str, link: FileLink) -> List[File]:
        """
        :param job_name:
        :param link:
        """
        files_list: List[File] = []
        for f in self.jobs_files[job_name]:
            if f.link == link:
                files_list.append(f)
        return files_list
