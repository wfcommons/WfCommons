#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from typing import Dict, List, Optional

from .abstract_recipe import WorkflowRecipe
from ...common.file import FileLink
from ...common.job import Job
from ...common.workflow import Workflow


class SeismologyRecipe(WorkflowRecipe):
    def __init__(self,
                 num_pairs: Optional[int],
                 data_size: Optional[int],
                 num_jobs: Optional[int]
                 ) -> None:
        """
        :param num_pairs: The number of pair of signals to estimate earthquake STFs.
        :param data_size:
        :param num_jobs:
        """
        super().__init__("Seismology", data_size, num_jobs)

        self.num_pairs: int = num_pairs

    @classmethod
    def from_num_pairs(cls, num_pairs: int) -> 'SeismologyRecipe':
        """
        :param num_pairs: The number of pair of signals to estimate earthquake STFs.
        """
        if num_pairs < 2:
            raise ValueError("The number of pairs should be at least 2.")

        return cls(num_pairs=num_pairs, data_size=None, num_jobs=None)

    def build_workflow(self, workflow_name: str = None) -> Workflow:
        """
        Build a synthetic trace of a Seismology workflow.
        :param workflow_name: workflow name
        """
        workflow = Workflow(name=self.name + "-synthetic-trace" if not workflow_name else workflow_name, makespan=None)
        self.job_id_counter: int = 1
        sg1iterdecon_jobs: List[Job] = []

        for _ in range(0, self.num_pairs):
            # sG1IterDecon job
            job_name = self._generate_job_name("sG1IterDecon")
            sg1iterdecon_job = self._generate_job('sG1IterDecon', job_name, files_recipe={FileLink.INPUT: {".lht": 2}})
            sg1iterdecon_jobs.append(sg1iterdecon_job)
            workflow.add_node(job_name, job=sg1iterdecon_job)

        # wrapper_siftSTFByMisfit job
        input_files = []
        for j in sg1iterdecon_jobs:
            for f in j.files:
                if f.link == FileLink.OUTPUT:
                    input_files.append(f)

        job_name = self._generate_job_name('wrapper_siftSTFByMisfit')
        wrapper_job = self._generate_job('wrapper_siftSTFByMisfit', job_name, input_files)
        workflow.add_node(job_name, job=wrapper_job)
        for j in sg1iterdecon_jobs:
            workflow.add_edge(j.name, wrapper_job.name)

        self.workflows.append(workflow)
        return workflow

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the Seismology workflow.
        """
        return {
            "sG1IterDecon": {
                "runtime": {
                    "min": 0.087,
                    "max": 5.615,
                    "distribution": {
                        "name": "alpha",
                        "params": [
                            2.8535577839854487e-09,
                            -0.6968250029499959,
                            1.0879675561652093
                        ]
                    }
                },
                "input": {
                    ".lht": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.4312877358390993,
                                -5.198983213279189e-26,
                                2.042713032349936
                            ]
                        },
                        "min": 1024,
                        "max": 16012
                    }
                },
                "output": {
                    ".stf": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                1.3444573433417438e-05,
                                -2922.408647942764,
                                6738.674391937242
                            ]
                        },
                        "min": 1144,
                        "max": 17016
                    }
                }
            },
            "wrapper_siftSTFByMisfit": {
                "runtime": {
                    "min": 0.089,
                    "max": 1.351,
                    "distribution": {
                        "name": "alpha",
                        "params": [
                            201.5475957749603,
                            -19.73655588794835,
                            6194.796244344304
                        ]
                    }
                },
                "input": {
                    ".stf": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                1.3444573433417438e-05,
                                -2922.408647942764,
                                6738.674391937242
                            ]
                        },
                        "min": 1144,
                        "max": 17016
                    },
                    ".py": {
                        "distribution": "None",
                        "min": 0,
                        "max": 0
                    },
                    "siftSTFByMisfit": {
                        "distribution": "None",
                        "min": 1386,
                        "max": 1386
                    }
                },
                "output": {
                    ".gz": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                10.99999999999958,
                                4.1986078970425886e-13
                            ]
                        },
                        "min": 63471,
                        "max": 687098
                    }
                }
            }
        }
