#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import math
import random

from typing import Dict, Optional

from .abstract_recipe import WorkflowRecipe
from ...common.file import FileLink
from ...common.workflow import Workflow


class EpigenomicsRecipe(WorkflowRecipe):
    """An Epigenomics workflow recipe class for creating synthetic workflow traces.

    :param num_sequence_files: Number of FASTQ files processed by the workflow.
    :type num_sequence_files: int
    :param num_lines: Number of lines in each FASTQ file.
    :type num_lines: int
    :param bin_size: Number of DNA and protein sequence information to be processed by each computational task.
    :type bin_size: int
    :param data_footprint: The upper bound for the workflow total data footprint (in bytes).
    :type data_footprint: int
    :param num_tasks: The upper bound for the total number of tasks in the workflow.
    :type num_tasks: int
    """

    def __init__(self,
                 num_sequence_files: Optional[int] = 1,
                 num_lines: Optional[int] = 10,
                 bin_size: Optional[int] = 10,
                 data_footprint: Optional[int] = 0,
                 num_tasks: Optional[int] = 9
                 ) -> None:
        """Create an object of the Epigenomics workflow recipe."""
        super().__init__("Epigenomics", data_footprint, num_tasks)

        self.num_sequence_files: int = num_sequence_files
        self.num_lines: int = num_lines
        self.bin_size: int = bin_size

    @classmethod
    def from_num_tasks(cls, num_tasks: int) -> 'EpigenomicsRecipe':
        """
        Instantiate an Epigenomics workflow recipe that will generate synthetic workflows
        up to the total number of tasks provided.

        :param num_tasks: The upper bound for the total number of tasks in the workflow (at least 9).
        :type num_tasks: int

        :return: An Epigenomics workflow recipe object that will generate synthetic workflows up
                 to the total number of tasks provided.
        :rtype: EpigenomicsRecipe
        """
        if num_tasks < 9:
            raise ValueError("The upper bound for the number of tasks should be at least 9.")

        num_sequence_files = random.randint(1, int(math.ceil((num_tasks - 5) / 12)))
        remaining_tasks = num_tasks - (6 * num_sequence_files) - 3
        num_lines = 1

        while remaining_tasks > 0:
            if remaining_tasks >= 4 * num_sequence_files:
                num_lines += 1
                remaining_tasks -= 4 * num_sequence_files
            else:
                break

        return cls(num_sequence_files=num_sequence_files, num_lines=num_lines * 10, bin_size=10,
                   data_footprint=None, num_tasks=num_tasks)

    @classmethod
    def from_sequences(cls,
                       num_sequence_files: int,
                       num_lines: int,
                       bin_size: int,
                       ) -> 'EpigenomicsRecipe':
        """
        Instantiate an Epigenomics workflow recipe that will generate synthetic workflows
        using the defined number of sequence files, lines, and bin size.

        :param num_sequence_files: Number of FASTQ files processed by the workflow.
        :type num_sequence_files: int
        :param num_lines: Number of lines in each FASTQ file.
        :type num_lines: int
        :param bin_size: Number of DNA and protein sequence information to be processed by each computational task.
        :type bin_size: int

        :return: An Epigenomics workflow recipe object that will generate synthetic workflows
                 using the defined number of sequence files, lines, and bin size.
        :rtype: EpigenomicsRecipe
        """
        if num_sequence_files < 1:
            raise ValueError("The number of sequence files should be at least 1.")
        if num_lines < 10:
            raise ValueError("The number of lines in per sequence file should be at least 10.")
        if bin_size < 10:
            raise ValueError("The bin size should be at least 10.")

        return cls(num_sequence_files=num_sequence_files, num_lines=num_lines, bin_size=bin_size, data_footprint=None,
                   num_tasks=None)

    def build_workflow(self, workflow_name: str = None) -> Workflow:
        """Generate a synthetic workflow trace of an Epigenomics workflow.

        :param workflow_name: The workflow name
        :type workflow_name: int

        :return: A synthetic workflow trace object.
        :rtype: Workflow
        """
        workflow = Workflow(name=self.name + "-synthetic-trace" if not workflow_name else workflow_name, makespan=None)
        self.task_id_counter: int = 1

        map_merge_tasks = []

        for _ in range(0, self.num_sequence_files):
            map_tasks = []

            # fastqsplit task
            num_pipelines = int(self.num_lines / self.bin_size)
            task_name = self._generate_task_name("fastqSplit")
            fastqsplit_task = self._generate_task('fastqSplit', task_name,
                                                  files_recipe={FileLink.OUTPUT: {".sfq": num_pipelines}})
            workflow.add_node(task_name, task=fastqsplit_task)

            for seq in range(0, num_pipelines):
                # filterContams task
                input_files = [fastqsplit_task.files[seq + 1]]
                task_name = self._generate_task_name("filterContams")
                filtercontams_task = self._generate_task('filterContams', task_name, input_files)
                workflow.add_node(task_name, task=filtercontams_task)
                workflow.add_edge(fastqsplit_task.name, filtercontams_task.name)

                # sol2sanger task
                input_files = [filtercontams_task.files[1]]
                task_name = self._generate_task_name("sol2sanger")
                sol2sanger_task = self._generate_task('sol2sanger', task_name, input_files)
                workflow.add_node(task_name, task=sol2sanger_task)
                workflow.add_edge(filtercontams_task.name, sol2sanger_task.name)

                # fast2bfq task
                input_files = [sol2sanger_task.files[3]]
                task_name = self._generate_task_name("fast2bfq")
                fast2bfq_task = self._generate_task('fast2bfq', task_name, input_files)
                workflow.add_node(task_name, task=fast2bfq_task)
                workflow.add_edge(sol2sanger_task.name, fast2bfq_task.name)

                # map task
                input_files = [fast2bfq_task.files[3]]
                task_name = self._generate_task_name("map")
                map_task = self._generate_task('map_', task_name, input_files)
                workflow.add_node(task_name, task=map_task)
                workflow.add_edge(fast2bfq_task.name, map_task.name)
                map_tasks.append(map_task)

            # map merge task (per sequence file)
            input_files = []
            for j in map_tasks:
                input_files.append(j.files[4])
            task_name = self._generate_task_name("mapMerge")
            map_merge_task = self._generate_task('mapMerge', task_name, input_files)
            workflow.add_node(task_name, task=map_merge_task)
            for j in map_tasks:
                workflow.add_edge(j.name, map_merge_task.name)
            map_merge_tasks.append(map_merge_task)

        # map merge task
        input_files = []
        for j in map_merge_tasks:
            for f in j.files:
                if f.link == FileLink.OUTPUT:
                    input_files.append(f)
        task_name = self._generate_task_name("mapMerge")
        map_merge_task = self._generate_task('mapMerge', task_name, input_files)
        workflow.add_node(task_name, task=map_merge_task)
        for j in map_merge_tasks:
            workflow.add_edge(j.name, map_merge_task.name)
        map_merge_tasks.append(map_merge_task)

        # chr21 task
        input_files = [map_merge_task.files[len(map_merge_tasks) + 1]]
        task_name = self._generate_task_name("chr21")
        chr21_task = self._generate_task('chr21', task_name, input_files)
        workflow.add_node(task_name, task=chr21_task)
        workflow.add_edge(map_merge_task.name, chr21_task.name)

        # pileup task
        input_files = [chr21_task.files[3]]
        task_name = self._generate_task_name("pileup")
        pileup_task = self._generate_task('pileup', task_name, input_files)
        workflow.add_node(task_name, task=pileup_task)
        workflow.add_edge(chr21_task.name, pileup_task.name)

        self.workflows.append(workflow)
        return workflow

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the Epigenomics workflow. Recipes can be
        generated by using the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`.

        :return: A recipe in the form of a dictionary in which keys are task prefixes.
        :rtype: Dict[str, Any]
        """
        return {
            "chr21": {
                "runtime": {
                    "min": 2.115,
                    "max": 120.272,
                    "distribution": {
                        "name": "alpha",
                        "params": [
                            1.2037882387638316e-08,
                            0.0013824676778065785,
                            1.4109916687079997
                        ]
                    }
                },
                "input": {
                    "maq": {
                        "distribution": "None",
                        "min": 171256,
                        "max": 171256
                    },
                    "maqindex": {
                        "distribution": "None",
                        "min": 118456,
                        "max": 118456
                    },
                    ".map": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.404004449642449,
                                7.999999999999998,
                                0.29643222001452685
                            ]
                        },
                        "min": 8974436,
                        "max": 213287821
                    }
                },
                "output": {
                    ".map": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.404004449642449,
                                7.999999999999998,
                                0.29643222001452685
                            ]
                        },
                        "min": 8974436,
                        "max": 213287821
                    }
                }
            },
            "fast2bfq": {
                "runtime": {
                    "min": 0.033,
                    "max": 10.513,
                    "distribution": {
                        "name": "chi2",
                        "params": [
                            1.3386691664974841,
                            -5.708893978145665e-28,
                            1.851059388057558
                        ]
                    }
                },
                "input": {
                    "maq": {
                        "distribution": "None",
                        "min": 171256,
                        "max": 171256
                    },
                    "maqindex": {
                        "distribution": "None",
                        "min": 118456,
                        "max": 118456
                    },
                    ".fq": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.7956524584974147,
                                -1.1631471106725046e-27,
                                1.7866384876482693
                            ]
                        },
                        "min": 290449,
                        "max": 10074424
                    }
                },
                "output": {
                    ".bfq": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.18924491041245323,
                                -8.340735410821612e-28,
                                2.241150785198654
                            ]
                        },
                        "min": 71135,
                        "max": 2779927
                    }
                }
            },
            "fastqSplit": {
                "runtime": {
                    "min": 1.344,
                    "max": 878.473,
                    "distribution": {
                        "name": "wald",
                        "params": [
                            -0.9319555092279244,
                            12.701940366609303
                        ]
                    }
                },
                "input": {
                    ".sfq": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                11.54431127713843,
                                1.7294835205944359
                            ]
                        },
                        "min": 109431824,
                        "max": 572332024
                    }
                },
                "output": {
                    ".sfq": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.6096730655283322,
                                -1.196048242820172e-25,
                                1.8460713490852902
                            ]
                        },
                        "min": 385300,
                        "max": 15387524
                    }
                }
            },
            "filterContams": {
                "runtime": {
                    "min": 0.067,
                    "max": 47.417,
                    "distribution": {
                        "name": "fisk",
                        "params": [
                            0.6541872041823478,
                            -4.336127245482357e-28,
                            1.938539282462673
                        ]
                    }
                },
                "input": {
                    ".sfq": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.6096730655283322,
                                -1.196048242820172e-25,
                                1.8460713490852902
                            ]
                        },
                        "min": 385300,
                        "max": 15387524
                    }
                },
                "output": {
                    ".sfq": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.577065936564831,
                                -1.0559908824984039e-26,
                                1.6405796858668555
                            ]
                        },
                        "min": 379368,
                        "max": 13605894
                    }
                }
            },
            "mapMerge": {
                "runtime": {
                    "min": 1.574,
                    "max": 37.253,
                    "distribution": {
                        "name": "alpha",
                        "params": [
                            2.5566780151876478e-08,
                            -0.00023397333904237204,
                            2.2836168403441324
                        ]
                    }
                },
                "input": {
                    "maq": {
                        "distribution": "None",
                        "min": 171256,
                        "max": 171256
                    },
                    "maqindex": {
                        "distribution": "None",
                        "min": 118456,
                        "max": 118456
                    },
                    ".map": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.36608100003038957,
                                -3.21850679716644e-26,
                                2.6277130589314748
                            ]
                        },
                        "min": 17429,
                        "max": 38332663
                    }
                },
                "output": {
                    ".map": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.46489086906995225,
                                1.9999999999999998,
                                1.555519982565344
                            ]
                        },
                        "min": 8972748,
                        "max": 213287821
                    }
                }
            },
            "map_": {
                "runtime": {
                    "min": 24.836,
                    "max": 122.642,
                    "distribution": {
                        "name": "trapz",
                        "params": [
                            1.0,
                            1.0,
                            -7.245000000000001,
                            82.8
                        ]
                    }
                },
                "input": {
                    "maq": {
                        "distribution": "None",
                        "min": 171256,
                        "max": 171256
                    },
                    "maqindex": {
                        "distribution": "None",
                        "min": 118456,
                        "max": 118456
                    },
                    ".bfq": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.18924491041245323,
                                -8.340735410821612e-28,
                                2.241150785198654
                            ]
                        },
                        "min": 71135,
                        "max": 2779927
                    },
                    ".bfa": {
                        "distribution": "None",
                        "min": 46944392,
                        "max": 46944392
                    }
                },
                "output": {
                    ".map": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                -0.3008407597354944,
                                0.9942059115945362
                            ]
                        },
                        "min": 17429,
                        "max": 1255720
                    }
                }
            },
            "pileup": {
                "runtime": {
                    "min": 23.112,
                    "max": 99.597,
                    "distribution": {
                        "name": "beta",
                        "params": [
                            1.312969391453295,
                            0.7331939054064136,
                            -8.103797068931092,
                            33.1037970689311
                        ]
                    }
                },
                "input": {
                    "maq": {
                        "distribution": "None",
                        "min": 171256,
                        "max": 171256
                    },
                    "maqindex": {
                        "distribution": "None",
                        "min": 118456,
                        "max": 118456
                    },
                    ".map": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.404004449642449,
                                7.999999999999998,
                                0.29643222001452685
                            ]
                        },
                        "min": 8974436,
                        "max": 213287821
                    },
                    ".bfa": {
                        "distribution": "None",
                        "min": 46944392,
                        "max": 46944392
                    }
                },
                "output": {
                    ".pileup": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.404004449642449,
                                7.999999999999998,
                                0.29643222001452685
                            ]
                        },
                        "min": 4595783,
                        "max": 88008977
                    }
                }
            },
            "sol2sanger": {
                "runtime": {
                    "min": 0.024,
                    "max": 28.568,
                    "distribution": {
                        "name": "chi2",
                        "params": [
                            0.5991026813436918,
                            -1.0785002971672537e-26,
                            3.083881675199752
                        ]
                    }
                },
                "input": {
                    "maq": {
                        "distribution": "None",
                        "min": 171256,
                        "max": 171256
                    },
                    "maqindex": {
                        "distribution": "None",
                        "min": 118456,
                        "max": 118456
                    },
                    ".sfq": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.577065936564831,
                                -1.0559908824984039e-26,
                                1.6405796858668555
                            ]
                        },
                        "min": 379368,
                        "max": 13605894
                    }
                },
                "output": {
                    ".fq": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.7956524584974147,
                                -1.1631471106725046e-27,
                                1.7866384876482693
                            ]
                        },
                        "min": 290449,
                        "max": 10074424
                    }
                }
            }
        }
