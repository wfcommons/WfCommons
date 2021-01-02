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


class BWARecipe(WorkflowRecipe):
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
        """Create an object of the BWA workflow recipe."""
        super().__init__("BWA",
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
                       ) -> 'BWARecipe':
        """
        Instantiate a BWA workflow recipe that will generate synthetic workflows up to
        the total number of tasks provided.

        :param num_tasks: The upper bound for the total number of tasks in the workflow (at least 6).
        :type num_tasks: int
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

        :return: A BWA workflow recipe object that will generate synthetic workflows up
                 to the total number of tasks provided.
        :rtype: BWARecipe
        """
        if num_tasks < 6:
            raise ValueError("The upper bound for the number of tasks should be at least 6.")

        return cls(num_subsample=int(num_tasks - 4),
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
                           ) -> 'BWARecipe':
        """
        Instantiate a BWA workflow recipe that will generate synthetic workflows using
        the defined number of subsample.

        :param num_subsample: The number of subsample the reference file will be split.
        :type num_subsample: int
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

        :return: A BWA workflow recipe object that will generate synthetic workflows
                 using the defined number of subsample.
        :rtype: BWARecipe
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
        Generate a synthetic workflow trace of a BWA workflow.

        :param workflow_name: The workflow name
        :type workflow_name: int

        :return: A synthetic workflow trace object.
        :rtype: Workflow
        """
        workflow = Workflow(name=self.name + '-synthetic-trace' if not workflow_name else workflow_name, makespan=None)
        self.task_id_counter: int = 1

        # fastq_reduce task
        task_name = self._generate_task_name('fastq_reduce')
        fastq_reduce_task = self._generate_task('fastq_reduce', task_name,
                                                files_recipe={FileLink.OUTPUT: {'.fastq': self.num_subsample}})
        workflow.add_node(task_name, task=fastq_reduce_task)

        # bwa_index task
        task_name = self._generate_task_name('bwa_index')
        bwa_index_task = self._generate_task('bwa_index', task_name,
                                             files_recipe={FileLink.OUTPUT: {'.fasta': self.num_subsample}})
        workflow.add_node(task_name, task=bwa_index_task)

        index_out_files = []
        for f in bwa_index_task.files:
            if f.link == FileLink.OUTPUT:
                index_out_files.append(f)

        out_files = []
        err_files = []
        bwa_tasks = []

        for f in fastq_reduce_task.files:
            if f.link == FileLink.OUTPUT:
                task_name = self._generate_task_name('bwa')
                files_list = index_out_files.copy()
                bwa_task = self._generate_task('bwa', task_name, files_list.append(f))
                workflow.add_node(task_name, task=bwa_task)
                workflow.add_edge(fastq_reduce_task.name, task_name)
                workflow.add_edge(bwa_index_task.name, task_name)
                bwa_tasks.append(task_name)
                for file in bwa_task.files:
                    if file.link == FileLink.OUTPUT and '.out' in file.name:
                        out_files.append(file)
                    if file.link == FileLink.OUTPUT and '.err' in file.name:
                        err_files.append(file)

        # cat_blast task
        task_name = self._generate_task_name('cat_bwa')
        cat_bwa_task = self._generate_task('cat_bwa', task_name,
                                           input_files=out_files,
                                           files_recipe={FileLink.OUTPUT: {'.out': self.num_subsample}})
        workflow.add_node(task_name, task=cat_bwa_task)

        # cat task
        task_name = self._generate_task_name('cat')
        cat_task = self._generate_task('cat', task_name,
                                       input_files=err_files,
                                       files_recipe={FileLink.OUTPUT: {'.err': self.num_subsample}})
        workflow.add_node(task_name, task=cat_task)

        for t in bwa_tasks:
            workflow.add_edge(t, cat_bwa_task.name)
            workflow.add_edge(t, cat_task.name)

        return workflow

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the BWA workflow. Recipes can be
        generated by using the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`.

        :return: A recipe in the form of a dictionary in which keys are task prefixes.
        :rtype: Dict[str, Any]
        """
        return {
            "fastq_reduce": {
                "runtime": {
                    "min": 0.045811,
                    "max": 4.724154,
                    "distribution": {
                        "name": "argus",
                        "params": [
                            1.599540430972508,
                            -0.5837001726233696,
                            1.703310458550342
                        ]
                    }
                },
                "input": {
                    "fastq_reduce": {
                        "distribution": "None",
                        "min": 2,
                        "max": 2
                    },
                    ".fastq": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                1.599540430972508,
                                -0.5837001726233696,
                                1.703310458550342
                            ]
                        },
                        "min": 2438,
                        "max": 247788
                    }
                },
                "output": {
                    ".fastq": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
                            ]
                        },
                        "min": 24,
                        "max": 249
                    }
                }
            },
            "bwa_index": {
                "runtime": {
                    "min": 80.652465,
                    "max": 1158.975678,
                    "distribution": {
                        "name": "arcsine",
                        "params": [
                            -0.2258070520586602,
                            1.2258070520586604
                        ]
                    }
                },
                "input": {
                    "bwa": {
                        "distribution": "None",
                        "min": 1445,
                        "max": 1445
                    },
                    ".fastq": {
                        "distribution": {
                            "name": "arcsine",
                            "params": [
                                -0.2258070520586602,
                                1.2258070520586604
                            ]
                        },
                        "min": 200438,
                        "max": 2047778
                    }
                },
                "output": {
                    ".bwt": {
                        "distribution": {
                            "name": "arcsine",
                            "params": [
                                -0.2258070520586602,
                                1.2258070520586604
                            ]
                        },
                        "min": 100001,
                        "max": 1000001
                    },
                    ".pac": {
                        "distribution": {
                            "name": "arcsine",
                            "params": [
                                -0.2258070520586602,
                                1.2258070520586604
                            ]
                        },
                        "min": 25001,
                        "max": 250001
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 1,
                        "max": 1
                    },
                    ".ann": {
                        "distribution": {
                            "name": "arcsine",
                            "params": [
                                -0.2258070520586602,
                                1.2258070520586604
                            ]
                        },
                        "min": 577,
                        "max": 61667
                    },
                    ".sa": {
                        "distribution": {
                            "name": "arcsine",
                            "params": [
                                -0.2258070520586602,
                                1.2258070520586604
                            ]
                        },
                        "min": 50001,
                        "max": 500001
                    }
                }
            },
            "bwa": {
                "runtime": {
                    "min": 0.18409,
                    "max": 34.912105,
                    "distribution": {
                        "name": "argus",
                        "params": [
                            1.787473923386625e-05,
                            -0.7146525478639822,
                            1.7383310890430437
                        ]
                    }
                },
                "input": {
                    "bwa": {
                        "distribution": "None",
                        "min": 1445,
                        "max": 1445
                    },
                    ".fastq": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 24,
                        "max": 2047778
                    },
                    ".bwt": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 100001,
                        "max": 1000001
                    },
                    ".pac": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 25001,
                        "max": 250001
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 1,
                        "max": 1
                    },
                    ".ann": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
                            ]
                        },
                        "min": 577,
                        "max": 61667
                    },
                    ".sa": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 50001,
                        "max": 500001
                    }
                },
                "output": {
                    ".sam": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 516,
                        "max": 52088
                    },
                    ".err": {
                        "distribution": "None",
                        "min": 1,
                        "max": 1
                    }
                }
            },
            "cat_bwa": {
                "runtime": {
                    "min": 0.540879,
                    "max": 519.460108,
                    "distribution": {
                        "name": "arcsine",
                        "params": [
                            -0.2258070520586602,
                            1.2258070520586604
                        ]
                    }
                },
                "input": {
                    "cat_bwa": {
                        "distribution": "None",
                        "min": 2,
                        "max": 2
                    },
                    ".sam": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 516,
                        "max": 52088
                    }
                },
                "output": {
                    ".sam": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                1.599540430972508,
                                -0.5837001726233696,
                                1.703310458550342
                            ]
                        },
                        "min": 3437,
                        "max": 353943
                    }
                }
            },
            "cat": {
                "runtime": {
                    "min": 0.014307,
                    "max": 0.073623,
                    "distribution": {
                        "name": "argus",
                        "params": [
                            1.599540430972508,
                            -0.5837001726233696,
                            1.703310458550342
                        ]
                    }
                },
                "input": {
                    ".err": {
                        "distribution": "None",
                        "min": 1,
                        "max": 1
                    }
                },
                "output": {
                    ".err": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                1.599540430972508,
                                -0.5837001726233696,
                                1.703310458550342
                            ]
                        },
                        "min": 14,
                        "max": 133
                    }
                }
            }
        }
