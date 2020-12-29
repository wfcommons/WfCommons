#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from typing import Dict, Optional

from .abstract_recipe import WorkflowRecipe
from ...common.file import FileLink
from ...common.task import Task
from ...common.workflow import Workflow


class SRASearchRecipe(WorkflowRecipe):
    """An SRA Search workflow recipe class for creating synthetic workflow traces.

    :param num_accession: The number of NCBI accession numbers.
    :type num_accession: int
    :param data_footprint: The upper bound for the workflow total data footprint (in bytes).
    :type data_footprint: int
    :param num_tasks: The upper bound for the total number of tasks in the workflow.
    :type num_tasks: int
    :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
    :type runtime_factor: float
    :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
    :type input_file_size_factor: float
    :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
    :type output_file_size_factor: float
    """

    def __init__(self,
                 num_accession: Optional[int] = 2,
                 data_footprint: Optional[int] = 0,
                 num_tasks: Optional[int] = 3,
                 runtime_factor: Optional[float] = 1.0,
                 input_file_size_factor: Optional[float] = 1.0,
                 output_file_size_factor: Optional[float] = 1.0
                 ) -> None:
        """Create an object of the SRA Search workflow recipe."""
        super().__init__("Seismology",
                         data_footprint,
                         num_tasks,
                         runtime_factor,
                         input_file_size_factor,
                         output_file_size_factor)

        self.num_accession: int = num_accession

    @classmethod
    def from_num_tasks(cls,
                       num_tasks: int,
                       runtime_factor: Optional[float] = 1.0,
                       input_file_size_factor: Optional[float] = 1.0,
                       output_file_size_factor: Optional[float] = 1.0
                       ) -> 'SRASearchRecipe':
        """
        Instantiate an SRA Search workflow recipe that will generate synthetic workflows up to
        the total number of tasks provided.

        :param num_tasks: The upper bound for the total number of tasks in the workflow (at least 6).
        :type num_tasks: int
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

        :return: An SRA Search workflow recipe object that will generate synthetic workflows up
                 to the total number of tasks provided.
        :rtype: SRASearchRecipe
        """
        if num_tasks < 6:
            raise ValueError("The upper bound for the number of tasks should be at least 6.")

        return cls(num_accession=int((num_tasks - 2) / 2),
                   data_footprint=None,
                   num_tasks=num_tasks,
                   runtime_factor=runtime_factor,
                   input_file_size_factor=input_file_size_factor,
                   output_file_size_factor=output_file_size_factor)

    @classmethod
    def from_num_accession(cls,
                           num_accession: int,
                           runtime_factor: Optional[float] = 1.0,
                           input_file_size_factor: Optional[float] = 1.0,
                           output_file_size_factor: Optional[float] = 1.0
                           ) -> 'SRASearchRecipe':
        """
        Instantiate an SRA Search workflow recipe that will generate synthetic workflows using
        the defined number of pairs.

        :param num_accession: The number of NCBI accession numbers.
        :type num_accession: int
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

        :return: An SRA Search workflow recipe object that will generate synthetic workflows
                 using the defined number of pairs.
        :rtype: SRASearchRecipe
        """
        if num_accession < 2:
            raise ValueError("The number of accessions should be at least 2.")

        return cls(num_accession=num_accession,
                   data_footprint=None,
                   num_tasks=None,
                   runtime_factor=runtime_factor,
                   input_file_size_factor=input_file_size_factor,
                   output_file_size_factor=output_file_size_factor)

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Generate a synthetic workflow trace of an SRA Search workflow.

        :param workflow_name: The workflow name
        :type workflow_name: int

        :return: A synthetic workflow trace object.
        :rtype: Workflow
        """
        workflow = Workflow(name=self.name + '-synthetic-trace' if not workflow_name else workflow_name, makespan=None)
        self.task_id_counter: int = 1

        # single bowtie2-build task
        task_name = self._generate_task_name('bowtie2-build')
        bowtie_build_task = self._generate_task('bowtie2-build', task_name,
                                                files_recipe={FileLink.OUTPUT: {'.bt2': 6}})
        workflow.add_node(task_name, task=bowtie_build_task)
        bowtie_files = []
        for f in bowtie_build_task.files:
            if f.link == FileLink.OUTPUT:
                bowtie_files.append(f)

        bowtie_tasks = []
        for _ in range(0, self.num_accession):
            # fasterq-dump task
            task_name = self._generate_task_name('fasterq-dump')
            fasterq_task = self._generate_task('fasterq-dump', task_name,
                                               files_recipe={FileLink.OUTPUT: {'.fastq': 2}})
            workflow.add_node(task_name, task=fasterq_task)

            # bowtie2 task
            input_files = bowtie_files.copy()
            input_files.extend(fasterq_task.files)
            task_name = self._generate_task_name('bowtie2')
            bowtie_task = self._generate_task('bowtie2', task_name, input_files)
            workflow.add_node(task_name, task=bowtie_task)

            workflow.add_edge(fasterq_task.name, bowtie_task.name)
            workflow.add_edge(bowtie_build_task.name, bowtie_task.name)

            bowtie_tasks.append(bowtie_task)

        # merge tasks
        input_files = []
        parents = []
        merge_tasks = []
        for bowtie_task in bowtie_tasks:
            if len(parents) < 25:
                parents.append(bowtie_task)
                for f in bowtie_task.files:
                    if f.link == FileLink.OUTPUT:
                        input_files.append(f)
            else:
                merge_tasks.append(self._add_merge_task(workflow, input_files, parents))
                parents = []
                input_files = []

        if not merge_tasks or len(parents) > 0:
            merge_tasks.append(self._add_merge_task(workflow, input_files, parents))

        parents = []
        input_files = []
        if len(merge_tasks) > 1:
            for m in merge_tasks:
                parents.append(m)
                for f in m.files:
                    if f.link == FileLink.OUTPUT:
                        input_files.append(f)
            self._add_merge_task(workflow, input_files, parents)

        return workflow

    def _add_merge_task(self, workflow, input_files, parents) -> Task:
        """Create a merge task.

        :param workflow: Workflow object instance.
        :rtype workflow: Workflow
        :param input_files: List of input files for the task.
        :rtype input_files: List[File]
        :param parents: List of parent tasks.
        :rtype parents: List[Task]

        :return: A merge task object.
        """
        task_name = self._generate_task_name('merge')
        merge_task = self._generate_task('merge', task_name, input_files)
        workflow.add_node(task_name, task=merge_task)
        for p in parents:
            workflow.add_edge(p.name, merge_task.name)
        return merge_task

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the SRA Search workflow. Recipes can be
        generated by using the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`.

        :return: A recipe in the form of a dictionary in which keys are task prefixes.
        :rtype: Dict[str, Any]
        """
        return {
            "bowtie2-build": {
                "runtime": {
                    "min": 0.316,
                    "max": 33.288,
                    "distribution": {
                        "name": "wald",
                        "params": [
                            -0.17799345401850686,
                            0.7891642364923879
                        ]
                    }
                },
                "input": {
                    ".fna": {
                        "distribution": "None",
                        "min": 98721,
                        "max": 98721
                    }
                },
                "output": {
                    ".bt2": {
                        "distribution": {
                            "name": "norm",
                            "params": [
                                0.1,
                                0.27080128015453203
                            ]
                        },
                        "min": 17,
                        "max": 4226858
                    }
                }
            },
            "fasterq-dump": {
                "runtime": {
                    "min": 1.173,
                    "max": 4343.699,
                    "distribution": {
                        "name": "fisk",
                        "params": [
                            3.9199026441968337,
                            -0.30421363341992647,
                            0.6309055633022219
                        ]
                    }
                },
                "input": {

                },
                "output": {
                    ".fastq": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                0.00015637238065969544,
                                -0.2508819176275605,
                                1.2531514622818194
                            ]
                        },
                        "min": 26168,
                        "max": 1981081890
                    }
                }
            },
            "bowtie2": {
                "runtime": {
                    "min": 2.377,
                    "max": 121.024,
                    "distribution": {
                        "name": "arcsine",
                        "params": [
                            -0.07358918218072019,
                            1.0735891821807204
                        ]
                    }
                },
                "input": {
                    ".bt2": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
                            ]
                        },
                        "min": 17,
                        "max": 4226858
                    },
                    ".fastq": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                0.00015637238065969544,
                                -0.2508819176275605,
                                1.2531514622818194
                            ]
                        },
                        "min": 26168,
                        "max": 1981081890
                    }
                },
                "output": {
                    ".bam": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750425928578598,
                                1.0,
                                1.3654380389187343e-20
                            ]
                        },
                        "min": 314,
                        "max": 317
                    },
                    ".bai": {
                        "distribution": "None",
                        "min": 24,
                        "max": 24
                    }
                }
            },
            "merge": {
                "runtime": {
                    "min": 0.039,
                    "max": 0.305,
                    "distribution": {
                        "name": "beta",
                        "params": [
                            0.12845511317891758,
                            0.19967338970784937,
                            -0.09952333240693918,
                            1.0995233324069393
                        ]
                    }
                },
                "input": {
                    ".bam": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750425928578598,
                                1.0,
                                1.3654380389187343e-20
                            ]
                        },
                        "min": 314,
                        "max": 317
                    },
                    ".bai": {
                        "distribution": "None",
                        "min": 24,
                        "max": 24
                    },
                    ".gz": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -1.000000000000004,
                                2.0000000000000044
                            ]
                        },
                        "min": 1540,
                        "max": 6175
                    }
                },
                "output": {
                    ".gz": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.25734524051070895,
                                0.2142857142857143,
                                0.3299110686952953
                            ]
                        },
                        "min": 1540,
                        "max": 10427
                    }
                }
            }
        }
