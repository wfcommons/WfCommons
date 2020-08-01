#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import networkx as nx

from typing import Dict, List, Optional
from .abstract_recipe import WorkflowRecipe
from ...common.file import FileLink
from ...common.job import Job
from ...common.workflow import Workflow


class GenomeRecipe(WorkflowRecipe):
    def __init__(self,
                 num_chromosomes: Optional[int],
                 num_sequences: Optional[int],
                 num_populations: Optional[int],
                 size: Optional[int],
                 num_jobs: Optional[int]
                 ) -> None:
        super().__init__("1000Genome")

        self.num_chromosomes: int = num_chromosomes
        self.num_sequences: int = num_sequences
        self.num_populations: int = num_populations
        self.size: int = size
        self.num_jobs: int = num_jobs

    @classmethod
    def from_num_chromosomes(cls,
                             num_chromosomes: int,
                             num_sequences: int,
                             num_populations: int,
                             ) -> 'GenomeRecipe':
        """
        :param num_chromosomes:
        :param num_sequences:
        :param num_populations:
        """
        if num_chromosomes < 1 or num_chromosomes > 22:
            raise ValueError("The number of chromosomes should be within the range [1,22].")
        if num_sequences < 1000:
            raise ValueError("The number of sequences should be at least 1000.")
        if num_populations < 1 or num_populations > 7:
            raise ValueError("The number of populations should be within the range [1,7].")

        return cls(num_chromosomes=num_chromosomes, num_sequences=num_sequences, num_populations=num_populations,
                   size=None, num_jobs=None)

    def build_workflow(self) -> Workflow:
        """

        """
        workflow = Workflow(name=self.name + "-synthetic-trace", makespan=None)
        job_id_counter: int = 1

        for _ in range(0, self.num_chromosomes):
            # individuals jobs
            individuals_jobs: List[Job] = []
            for _ in range(0, int(self.num_sequences / 1000)):
                job_name = "individuals_{:08d}".format(job_id_counter)
                individuals_job = self._generate_job('individuals', job_name, None)
                individuals_jobs.append(individuals_job)
                workflow.add_node(job_name, job=individuals_job)
                job_id_counter += 1

            # individuals merge job
            job_name = "individuals_merge_{:08d}".format(job_id_counter)
            input_files = []
            for j in individuals_jobs:
                input_files.extend(self._get_files_by_job_and_link(j.name, FileLink.OUTPUT))
            individuals_merge_job = self._generate_job('individuals_merge', job_name, input_files)
            workflow.add_node(job_name, job=individuals_merge_job)
            for j in individuals_jobs:
                workflow.add_edge(j.name, individuals_merge_job.name)
            job_id_counter += 1

            # sifting job
            job_name = "sifting_{:08d}".format(job_id_counter)
            sifting_job = self._generate_job('sifting', job_name, None)
            workflow.add_node(job_name, job=sifting_job)
            job_id_counter += 1

            populations = ['ALL', 'AFR', 'AMR', 'EAS', 'EUR', 'GBR', 'SAS']

            # mutation overlap jobs
            input_files = self._get_files_by_job_and_link(individuals_merge_job.name, FileLink.OUTPUT)
            input_files.extend(self._get_files_by_job_and_link(sifting_job.name, FileLink.OUTPUT))
            for p in range(0, self.num_populations):
                job_name = "mutation_overlap_{:08d}".format(job_id_counter)
                input_files.append(self._generate_file(populations[p],
                                                       self._workflow_recipe()['mutation_overlap']['input'],
                                                       FileLink.INPUT))
                mutation_overlap_job = self._generate_job('mutation_overlap', job_name, input_files)
                workflow.add_node(job_name, job=mutation_overlap_job)
                workflow.add_edge(sifting_job.name, job_name)
                workflow.add_edge(individuals_merge_job.name, job_name)
                job_id_counter += 1

            # frequency jobs
            input_files = self._get_files_by_job_and_link(individuals_merge_job.name, FileLink.OUTPUT)
            input_files.extend(self._get_files_by_job_and_link(sifting_job.name, FileLink.OUTPUT))
            for p in range(0, self.num_populations):
                job_name = "frequency_{:08d}".format(job_id_counter)
                input_files.append(self._generate_file(populations[p],
                                                       self._workflow_recipe()['frequency']['input'],
                                                       FileLink.INPUT))
                frequency_job = self._generate_job('frequency', job_name, input_files)
                workflow.add_node(job_name, job=frequency_job)
                workflow.add_edge(sifting_job.name, job_name)
                workflow.add_edge(individuals_merge_job.name, job_name)
                job_id_counter += 1

        # TODO: remove printing below
        for j in workflow.nodes:
            print("JOB: {}".format(workflow.nodes[j]['job'].name))
        output = "{0}.dot".format(self.name)
        nx.nx_agraph.write_dot(workflow, output)

        self.workflows.append(workflow)
        return workflow

    def _workflow_recipe(self) -> Dict:
        """

        """
        return {
            "individuals": {
                "runtime": {
                    "min": 48.846,
                    "max": 192.232,
                    "distribution": {
                        "name": "skewnorm",
                        "params": [
                            11115267.652937062,
                            -2.9628504044929433e-05,
                            56.03957070238482
                        ]
                    }
                },
                "input": {
                    ".txt": {
                        "distribution": "None",
                        "min": 20078,
                        "max": 20078
                    },
                    ".vcf": {
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                3.377581164717659,
                                -1.7883202273916273e-16,
                                8.56622206451036e-16
                            ]
                        },
                        "min": 1014213207,
                        "max": 2540066220
                    }
                },
                "output": {
                    ".gz": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                -1.0888749504987145,
                                6.952913441525787
                            ]
                        },
                        "min": 27275,
                        "max": 29311
                    }
                }
            },
            "individuals_merge": {
                "runtime": {
                    "min": 34.471,
                    "max": 157.346,
                    "distribution": {
                        "name": "trapz",
                        "params": [
                            1.0,
                            1.0,
                            -8.190000000000001,
                            105.6
                        ]
                    }
                },
                "input": {
                    ".gz": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                -1.0888749504987145,
                                6.952913441525787
                            ]
                        },
                        "min": 27275,
                        "max": 29311
                    }
                },
                "output": {
                    ".gz": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5579424462170595,
                                -7.62317369524941e-27,
                                2.3113046845460117
                            ]
                        },
                        "min": 24180,
                        "max": 27336
                    }
                }
            },
            "sifting": {
                "runtime": {
                    "min": 0.293,
                    "max": 22.686,
                    "distribution": {
                        "name": "chi2",
                        "params": [
                            1.3512954025717425,
                            0.9999999999999999,
                            1.9181211494208956
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                -1.6067534129634995,
                                9.125715788595995
                            ]
                        },
                        "min": 267806263,
                        "max": 1708156254
                    }
                },
                "output": {
                    ".txt": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                -1.2088789973902925,
                                4.40384271471094
                            ]
                        },
                        "min": 231958,
                        "max": 2126612
                    }
                }
            },
            "mutation_overlap": {
                "runtime": {
                    "min": 1.61,
                    "max": 92.154,
                    "distribution": {
                        "name": "alpha",
                        "params": [
                            1.7105772737618344e-07,
                            -1.5157641888455222,
                            2.675317958811103
                        ]
                    }
                },
                "input": {
                    ".gz": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.6569480734384281,
                                -5.5283127064899054e-27,
                                1.9520478968194923
                            ]
                        },
                        "min": 24180,
                        "max": 27336
                    },
                    ".txt": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.42634352667168895,
                                -3.39060287009396e-26,
                                3.095939096458289
                            ]
                        },
                        "min": 20078,
                        "max": 2126612
                    },
                    "AFR": {
                        "distribution": "None",
                        "min": 8088,
                        "max": 8088
                    },
                    "GBR": {
                        "distribution": "None",
                        "min": 856,
                        "max": 856
                    },
                    "ALL": {
                        "distribution": "None",
                        "min": 28000,
                        "max": 28000
                    },
                    "SAS": {
                        "distribution": "None",
                        "min": 5248,
                        "max": 5248
                    },
                    "EAS": {
                        "distribution": "None",
                        "min": 4896,
                        "max": 4896
                    },
                    "AMR": {
                        "distribution": "None",
                        "min": 4248,
                        "max": 4248
                    },
                    "EUR": {
                        "distribution": "None",
                        "min": 5312,
                        "max": 5312
                    }
                },
                "output": {
                    ".gz": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.1415750994331953,
                                -1.5084913431314967e-26,
                                2.267414573651914
                            ]
                        },
                        "min": 136896,
                        "max": 172605
                    }
                }
            },
            "frequency": {
                "runtime": {
                    "min": 76.023,
                    "max": 238.26,
                    "distribution": {
                        "name": "levy",
                        "params": [
                            -0.8146765116317367,
                            3.4571998225851384
                        ]
                    }
                },
                "input": {
                    ".gz": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.6569480734384281,
                                -5.5283127064899054e-27,
                                1.9520478968194923
                            ]
                        },
                        "min": 24180,
                        "max": 27336
                    },
                    ".txt": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.42634352667168895,
                                -3.39060287009396e-26,
                                3.095939096458289
                            ]
                        },
                        "min": 20078,
                        "max": 2126612
                    },
                    "AFR": {
                        "distribution": "None",
                        "min": 8088,
                        "max": 8088
                    },
                    "GBR": {
                        "distribution": "None",
                        "min": 856,
                        "max": 856
                    },
                    "ALL": {
                        "distribution": "None",
                        "min": 28000,
                        "max": 28000
                    },
                    "SAS": {
                        "distribution": "None",
                        "min": 5248,
                        "max": 5248
                    },
                    "EAS": {
                        "distribution": "None",
                        "min": 4896,
                        "max": 4896
                    },
                    "AMR": {
                        "distribution": "None",
                        "min": 4248,
                        "max": 4248
                    },
                    "EUR": {
                        "distribution": "None",
                        "min": 5312,
                        "max": 5312
                    }
                },
                "output": {
                    ".gz": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                -0.13019549190430535,
                                0.4034357322776298
                            ]
                        },
                        "min": 217409,
                        "max": 298487
                    }
                }
            }
        }
