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
from ...common.workflow import Workflow


class BLASTRecipe(WorkflowRecipe):
    """A BLAST workflow recipe class for creating synthetic workflow traces.

    :param num_subsample: The number of subsample the reference file will be split.
    :type num_subsample: int
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
                 num_subsample: Optional[int] = 2,
                 data_footprint: Optional[int] = 0,
                 num_tasks: Optional[int] = 5,
                 runtime_factor: Optional[float] = 1.0,
                 input_file_size_factor: Optional[float] = 1.0,
                 output_file_size_factor: Optional[float] = 1.0
                 ) -> None:
        """Create an object of the BLAST workflow recipe."""
        super().__init__("BLAST",
                         data_footprint,
                         num_tasks,
                         runtime_factor,
                         input_file_size_factor,
                         output_file_size_factor)

        self.num_subsample: int = num_subsample

    @classmethod
    def from_num_tasks(cls,
                       num_tasks: int,
                       runtime_factor: Optional[float] = 1.0,
                       input_file_size_factor: Optional[float] = 1.0,
                       output_file_size_factor: Optional[float] = 1.0
                       ) -> 'BLASTRecipe':
        """
        Instantiate a BLAST workflow recipe that will generate synthetic workflows up to
        the total number of tasks provided.

        :param num_tasks: The upper bound for the total number of tasks in the workflow (at least 5).
        :type num_tasks: int
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

        :return: A BLAST workflow recipe object that will generate synthetic workflows up
                 to the total number of tasks provided.
        :rtype: BLASTRecipe
        """
        if num_tasks < 5:
            raise ValueError("The upper bound for the number of tasks should be at least 5.")

        return cls(num_subsample=int(num_tasks - 3),
                   data_footprint=None,
                   num_tasks=num_tasks,
                   runtime_factor=runtime_factor,
                   input_file_size_factor=input_file_size_factor,
                   output_file_size_factor=output_file_size_factor)

    @classmethod
    def from_num_subsample(cls,
                           num_subsample: int,
                           runtime_factor: Optional[float] = 1.0,
                           input_file_size_factor: Optional[float] = 1.0,
                           output_file_size_factor: Optional[float] = 1.0
                           ) -> 'BLASTRecipe':
        """
        Instantiate a BLAST workflow recipe that will generate synthetic workflows using
        the defined number of subsample.

        :param num_subsample: The number of subsample the reference file will be split.
        :type num_subsample: int
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

        :return: A BLAST workflow recipe object that will generate synthetic workflows
                 using the defined number of subsample.
        :rtype: BLASTRecipe
        """
        if num_subsample < 2:
            raise ValueError("The number of subsample should be at least 2.")

        return cls(num_subsample=num_subsample,
                   data_footprint=None,
                   num_tasks=None,
                   runtime_factor=runtime_factor,
                   input_file_size_factor=input_file_size_factor,
                   output_file_size_factor=output_file_size_factor)

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Generate a synthetic workflow trace of a BLAST workflow.

        :param workflow_name: The workflow name
        :type workflow_name: int

        :return: A synthetic workflow trace object.
        :rtype: Workflow
        """
        workflow = Workflow(name=self.name + '-synthetic-trace' if not workflow_name else workflow_name, makespan=None)
        self.task_id_counter: int = 1

        # split_fasta task
        task_name = self._generate_task_name('split_fasta')
        split_fasta_task = self._generate_task('split_fasta', task_name,
                                               files_recipe={FileLink.OUTPUT: {'.fasta': self.num_subsample}})
        workflow.add_node(task_name, task=split_fasta_task)

        out_files = []
        err_files = []
        blastall_tasks = []

        for f in split_fasta_task.files:
            if f.link == FileLink.OUTPUT:
                task_name = self._generate_task_name('blastall')
                blastall_task = self._generate_task('blastall', task_name, [f])
                workflow.add_node(task_name, task=blastall_task)
                workflow.add_edge(split_fasta_task.name, task_name)
                blastall_tasks.append(task_name)
                for file in blastall_task.files:
                    if file.link == FileLink.OUTPUT and '.out' in file.name:
                        out_files.append(file)
                    if file.link == FileLink.OUTPUT and '.err' in file.name:
                        err_files.append(file)

        # cat_blast task
        task_name = self._generate_task_name('cat_blast')
        cat_blast_task = self._generate_task('cat_blast', task_name,
                                             input_files=out_files,
                                             files_recipe={FileLink.OUTPUT: {'.out': self.num_subsample}})
        workflow.add_node(task_name, task=cat_blast_task)

        # cat task
        task_name = self._generate_task_name('cat')
        cat_task = self._generate_task('cat', task_name,
                                       input_files=err_files,
                                       files_recipe={FileLink.OUTPUT: {'.err': self.num_subsample}})
        workflow.add_node(task_name, task=cat_task)

        for t in blastall_tasks:
            workflow.add_edge(t, cat_blast_task.name)
            workflow.add_edge(t, cat_task.name)

        return workflow

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the BLAST workflow. Recipes can be
        generated by using the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`.

        :return: A recipe in the form of a dictionary in which keys are task prefixes.
        :rtype: Dict[str, Any]
        """
        return {
            "split_fasta": {
                "runtime": {
                    "min": 0.051992,
                    "max": 3.160018,
                    "distribution": {
                        "name": "arcsine",
                        "params": [
                            -0.2258070520586602,
                            1.2258070520586604
                        ]
                    }
                },
                "input": {
                    ".fasta": {
                        "distribution": {
                            "name": "arcsine",
                            "params": [
                                -0.2258070520586602,
                                1.2258070520586604
                            ]
                        },
                        "min": 203,
                        "max": 201389
                    },
                    "split_fasta": {
                        "distribution": "None",
                        "min": 1,
                        "max": 1
                    }
                },
                "output": {
                    ".fasta": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
                            ]
                        },
                        "min": 6,
                        "max": 2015
                    }
                }
            },
            "blastall": {
                "runtime": {
                    "min": 8.116334,
                    "max": 1799.556624,
                    "distribution": {
                        "name": "trapz",
                        "params": [
                            1.0,
                            1.0,
                            -0.10500000000000001,
                            1.2
                        ]
                    }
                },
                "input": {
                    "blastall": {
                        "distribution": "None",
                        "min": 7688,
                        "max": 7688
                    },
                    ".fasta": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
                            ]
                        },
                        "min": 6,
                        "max": 2015
                    },
                    "nt": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
                            ]
                        },
                        "min": 0,
                        "max": 5117704493
                    }
                },
                "output": {
                    ".out": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                2.465535551931572e-05,
                                -0.7452662890705088,
                                1.7477663092998088
                            ]
                        },
                        "min": 5,
                        "max": 17952
                    },
                    ".err": {
                        "distribution": "None",
                        "min": 0,
                        "max": 0
                    }
                }
            },
            "cat_blast": {
                "runtime": {
                    "min": 0.034811,
                    "max": 16.689957,
                    "distribution": {
                        "name": "arcsine",
                        "params": [
                            -0.2258070520586602,
                            1.2258070520586604
                        ]
                    }
                },
                "input": {
                    "cat_blast": {
                        "distribution": "None",
                        "min": 1,
                        "max": 1
                    },
                    ".out": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                2.465535551931572e-05,
                                -0.7452662890705088,
                                1.7477663092998088
                            ]
                        },
                        "min": 5,
                        "max": 17952
                    }
                },
                "output": {
                    "None": {
                        "distribution": {
                            "name": "arcsine",
                            "params": [
                                -0.2258070520586602,
                                1.2258070520586604
                            ]
                        },
                        "min": 454,
                        "max": 565948
                    }
                }
            },
            "cat": {
                "runtime": {
                    "min": 0.009596,
                    "max": 0.021895,
                    "distribution": {
                        "name": "arcsine",
                        "params": [
                            -0.2258070520586602,
                            1.2258070520586604
                        ]
                    }
                },
                "input": {
                    ".err": {
                        "distribution": "None",
                        "min": 0,
                        "max": 0
                    }
                },
                "output": {
                    ".err": {
                        "distribution": "None",
                        "min": 0,
                        "max": 0
                    }
                }
            }
        }
