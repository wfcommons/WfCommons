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


class GenomeRecipe(WorkflowRecipe):
    """A 1000Genome workflow recipe class for creating synthetic workflow traces.

    :param num_chromosomes: The number of chromosomes evaluated in the workflow execution.
    :type num_chromosomes: int
    :param num_sequences: The number of sequences per chromosome file.
    :type num_sequences: int
    :param num_populations: The number of populations being evaluated.
    :type num_populations: int
    :param data_footprint: The upper bound for the workflow total data footprint (in bytes).
    :type data_footprint: int
    :param num_jobs: The upper bound for the total number of jobs in the workflow.
    :type num_jobs: int
    """

    def __init__(self,
                 num_chromosomes: Optional[int] = 1,
                 num_sequences: Optional[int] = 1,
                 num_populations: Optional[int] = 1,
                 data_footprint: Optional[int] = 0,
                 num_jobs: Optional[int] = 5
                 ) -> None:
        """Create an object of the 1000Genome workflow recipe."""
        super().__init__("1000Genome", data_footprint, num_jobs)

        self.num_chromosomes: int = num_chromosomes
        self.num_sequences: int = num_sequences
        self.num_populations: int = num_populations
        self.populations = ['ALL', 'AFR', 'AMR', 'EAS', 'EUR', 'GBR', 'SAS']

    @classmethod
    def from_num_jobs(cls, num_jobs: int) -> 'GenomeRecipe':
        """
        Instantiate a 1000Genome workflow recipe that will generate synthetic workflows up to
        the total number of jobs provided.

        :param num_jobs: The upper bound for the total number of jobs in the workflow (at least 5).
        :type num_jobs: int

        :return: A 1000Genome workflow recipe object that will generate synthetic workflows up
                 to the total number of jobs provided.
        :rtype: GenomeRecipe
        """
        if num_jobs < 5:
            raise ValueError("The upper bound for the number of jobs should be at least 5.")

        num_chromosomes = 1
        num_sequences = 1
        num_populations = 1
        remaining_jobs = num_jobs - 5

        while remaining_jobs > 0:
            added_job = False
            if num_sequences <= num_populations or 2 > remaining_jobs >= num_chromosomes:
                num_sequences += 1
                remaining_jobs -= 1
                added_job = True
            if 7 > num_populations < num_sequences + 1 and remaining_jobs >= 2:
                num_populations += 1
                remaining_jobs -= 2
                added_job = True
            if not added_job:
                jobs_per_branch = 2 + num_populations * 2 + num_sequences
                if num_populations == 7 and num_chromosomes < 22 and remaining_jobs >= jobs_per_branch:
                    num_chromosomes += 1
                    remaining_jobs -= jobs_per_branch
                elif jobs_per_branch > remaining_jobs >= num_chromosomes or \
                        remaining_jobs >= num_chromosomes:
                    num_sequences += 1
                    remaining_jobs -= num_chromosomes
                    jobs_per_branch += 1
                else:
                    break

        return cls(num_chromosomes=num_chromosomes, num_sequences=num_sequences * 1000, num_populations=num_populations,
                   data_footprint=None, num_jobs=num_jobs)

    @classmethod
    def from_num_chromosomes(cls,
                             num_chromosomes: int,
                             num_sequences: int,
                             num_populations: int,
                             ) -> 'GenomeRecipe':
        """
        Instantiate a 1000Genome workflow recipe that will generate synthetic workflows using
        the defined number of chromosomes, sequences, and populations.

        :param num_chromosomes: The number of chromosomes evaluated in the workflow execution.
        :type num_chromosomes: int
        :param num_sequences: The number of sequences per chromosome file.
        :type num_sequences: int
        :param num_populations: The number of populations being evaluated.
        :type num_populations: int

        :return: A 1000Genome workflow recipe object that will generate synthetic workflows
                 using the defined number of chromosomes, sequences, and populations.
        :rtype: GenomeRecipe
        """
        if num_chromosomes < 1 or num_chromosomes > 22:
            raise ValueError("The number of chromosomes should be within the range [1,22].")
        if num_sequences < 1000:
            raise ValueError("The number of sequences should be at least 1000.")
        if num_populations < 1 or num_populations > 7:
            raise ValueError("The number of populations should be within the range [1,7].")

        return cls(num_chromosomes=num_chromosomes, num_sequences=num_sequences, num_populations=num_populations,
                   data_footprint=None, num_jobs=None)

    def build_workflow(self, workflow_name: str = None) -> Workflow:
        """Generate a synthetic workflow trace of a 1000Genome workflow.

        :param workflow_name: The workflow name
        :type workflow_name: int

        :return: A synthetic workflow trace object.
        :rtype: Workflow
        """
        workflow = Workflow(name=self.name + "-synthetic-trace" if not workflow_name else workflow_name, makespan=None)
        self.job_id_counter: int = 1

        for _ in range(0, self.num_chromosomes):
            # individuals jobs
            individuals_jobs: List[Job] = []
            for _ in range(0, int(self.num_sequences / 1000)):
                job_name = self._generate_job_name("individuals")
                individuals_job = self._generate_job('individuals', job_name)
                individuals_jobs.append(individuals_job)
                workflow.add_node(job_name, job=individuals_job)

            # individuals merge job
            job_name = self._generate_job_name("individuals_merge")
            input_files = []
            for j in individuals_jobs:
                input_files.extend(self._get_files_by_job_and_link(j.name, FileLink.OUTPUT))
            individuals_merge_job = self._generate_job('individuals_merge', job_name, input_files)
            workflow.add_node(job_name, job=individuals_merge_job)
            for j in individuals_jobs:
                workflow.add_edge(j.name, individuals_merge_job.name)

            # sifting job
            job_name = self._generate_job_name("sifting")
            sifting_job = self._generate_job('sifting', job_name)
            workflow.add_node(job_name, job=sifting_job)

            # mutation overlap jobs
            input_files = self._get_files_by_job_and_link(individuals_merge_job.name, FileLink.OUTPUT)
            input_files.extend(self._get_files_by_job_and_link(sifting_job.name, FileLink.OUTPUT))
            for p in range(0, self.num_populations):
                job_name = self._generate_job_name("mutation_overlap")
                mutation_overlap_job = self._generate_job('mutation_overlap', job_name, input_files,
                                                          files_recipe=self._get_populations_files_recipe(p))
                workflow.add_node(job_name, job=mutation_overlap_job)
                workflow.add_edge(sifting_job.name, job_name)
                workflow.add_edge(individuals_merge_job.name, job_name)

            # frequency jobs
            input_files = self._get_files_by_job_and_link(individuals_merge_job.name, FileLink.OUTPUT)
            input_files.extend(self._get_files_by_job_and_link(sifting_job.name, FileLink.OUTPUT))
            for p in range(0, self.num_populations):
                job_name = self._generate_job_name("frequency")
                frequency_job = self._generate_job('frequency', job_name, input_files,
                                                   files_recipe=self._get_populations_files_recipe(p))
                workflow.add_node(job_name, job=frequency_job)
                workflow.add_edge(sifting_job.name, job_name)
                workflow.add_edge(individuals_merge_job.name, job_name)

        self.workflows.append(workflow)
        return workflow

    def _get_populations_files_recipe(self, index: int) -> Dict[FileLink, Dict[str, int]]:
        """Get the recipe for generating a population file.

        :param index: Index of the population in the list.
        :type index: int

        :return: Recipe for generating a population file.
        :rtype: Dict[FileLink, Dict[str, int]]
        """
        recipe = {FileLink.INPUT: {}}
        idx_count = 0
        for pop in self.populations:
            if idx_count == index:
                recipe[FileLink.INPUT][pop] = 1
            else:
                recipe[FileLink.INPUT][pop] = 0
            idx_count += 1
        return recipe

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the 1000Genome workflow. Recipes can be
        generated by using the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`.

        :return: A recipe in the form of a dictionary in which keys are job prefixes.
        :rtype: Dict[str, Any]
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
