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

from typing import Dict, List, Optional

from .abstract_recipe import WorkflowRecipe
from ...common.file import FileLink
from ...common.job import Job
from ...common.workflow import Workflow


class SoyKBRecipe(WorkflowRecipe):
    """A SoyKB workflow recipe class for creating synthetic workflow traces.

    :param num_fastq_files: The number of FASTQ files to be analyzed.
    :type num_fastq_files: int
    :param num_chromosomes: The number of chromosomes.
    :type num_chromosomes: int
    :param data_footprint: The upper bound for the workflow total data footprint (in bytes).
    :type data_footprint: int
    :param num_jobs: The upper bound for the total number of jobs in the workflow.
    :type num_jobs: int
    """

    def __init__(self,
                 num_fastq_files: Optional[int] = 2,
                 num_chromosomes: Optional[int] = 1,
                 data_footprint: Optional[int] = 0,
                 num_jobs: Optional[int] = 14
                 ) -> None:
        """Create an object of the SoyKB workflow recipe."""
        super().__init__("SoyKB", data_footprint, num_jobs)

        self.num_fastq_files = num_fastq_files
        self.num_chromosomes = num_chromosomes

    @classmethod
    def from_num_jobs(cls, num_jobs: int) -> 'SoyKBRecipe':
        """
        Instantiate a SoyKB workflow recipe that will generate synthetic workflows up to
        the total number of jobs provided.

        :param num_jobs: The upper bound for the total number of jobs in the workflow (at least 14).
        :type num_jobs: int

        :return: A SoyKB workflow recipe object that will generate synthetic workflows up
                 to the total number of jobs provided.
        :rtype: SoyKBRecipe
        """
        if num_jobs < 14:
            raise ValueError("The upper bound for the number of jobs should be at least 14.")

        num_chromosomes = 1 if num_jobs == 14 else random.randint(1, min(math.ceil((num_jobs - 14) / 2), 22))
        remaining_jobs = num_jobs - (2 * num_chromosomes) - 12
        num_fastq_files = 1

        while remaining_jobs > 0:
            if remaining_jobs >= num_chromosomes + 6:
                num_fastq_files += 1
                remaining_jobs -= num_chromosomes + 6
            else:
                break

        return cls(num_fastq_files=num_fastq_files * 2, num_chromosomes=num_chromosomes, data_footprint=None,
                   num_jobs=num_jobs)

    @classmethod
    def from_sequences(cls, num_fastq_files: int, num_chromosomes: int) -> 'SoyKBRecipe':
        """
        Instantiate a SoyKB workflow recipe that will generate synthetic workflows using
        the defined number of FASTQ files and chromosomes.

        :param num_fastq_files: The number of FASTQ files to be analyzed (at least 2).
        :type num_fastq_files: int
        :param num_chromosomes: The number of chromosomes (range [1,22].
        :type num_chromosomes: int

        :return: A SoyKB workflow recipe object that will generate synthetic workflows
                 using the defined number of FASTQ files and chromosomes.
        :rtype: SoyKBRecipe
        """
        if num_fastq_files < 2 or num_fastq_files % 2 != 0:
            raise ValueError("The number of FASTQ files should be at least 2 and should be an even number.")
        if num_chromosomes < 1 or num_chromosomes > 22:
            raise ValueError("The number of chromosomes should be within range [1, 22].")

        return cls(num_fastq_files=num_fastq_files, num_chromosomes=num_chromosomes, data_footprint=None, num_jobs=None)

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """Generate a synthetic workflow trace of a SoyKB workflow.

        :param workflow_name: The workflow name
        :type workflow_name: int

        :return: A synthetic workflow trace object.
        :rtype: Workflow
        """
        workflow = Workflow(name=self.name + "-synthetic-trace" if not workflow_name else workflow_name, makespan=None)
        self.job_id_counter: int = 1
        haplotype_caller_jobs: Dict[int, List[Job]] = {}
        genotype_gvcfs_jobs: List[Job] = []

        num_pipelines: int = math.ceil(self.num_fastq_files / 2)

        for _ in range(0, num_pipelines):
            # align to reference job
            job_name = self._generate_job_name('alignment_to_reference')
            atr_job = self._generate_job('alignment_to_reference', job_name,
                                         files_recipe={FileLink.INPUT: {".fastq": 2}})
            workflow.add_node(job_name, job=atr_job)

            # sort sam job
            job_name = self._generate_job_name('sort_sam')
            sort_sam_job = self._generate_job('sort_sam', job_name, input_files=[atr_job.files[12]])
            workflow.add_node(job_name, job=sort_sam_job)
            workflow.add_edge(atr_job.name, sort_sam_job.name)

            # dedup job
            job_name = self._generate_job_name('dedup')
            dedup_job = self._generate_job('dedup', job_name,
                                           input_files=[sort_sam_job.files[2], sort_sam_job.files[3]])
            workflow.add_node(job_name, job=dedup_job)
            workflow.add_edge(sort_sam_job.name, dedup_job.name)

            # add replace job
            job_name = self._generate_job_name('add_replace')
            add_replace_job = self._generate_job('add_replace', job_name,
                                                 input_files=[dedup_job.files[3], dedup_job.files[4]])
            workflow.add_node(job_name, job=add_replace_job)
            workflow.add_edge(dedup_job.name, add_replace_job.name)

            # realign target creator replace job
            job_name = self._generate_job_name('realign_target_creator')
            rtc_job = self._generate_job('realign_target_creator', job_name,
                                         input_files=[add_replace_job.files[3], add_replace_job.files[4]])
            workflow.add_node(job_name, job=rtc_job)
            workflow.add_edge(add_replace_job.name, rtc_job.name)

            # indel realign job
            job_name = self._generate_job_name('indel_realign')
            indel_realign_job = self._generate_job('indel_realign', job_name, input_files=[rtc_job.files[12]])
            workflow.add_node(job_name, job=indel_realign_job)
            workflow.add_edge(rtc_job.name, indel_realign_job.name)

            for ch in range(0, self.num_chromosomes):
                # haplotype caller job
                job_name = self._generate_job_name('haplotype_caller')
                haplotype_caller_job = self._generate_job('haplotype_caller', job_name,
                                                          input_files=[indel_realign_job.files[13],
                                                                       indel_realign_job.files[14]])
                workflow.add_node(job_name, job=haplotype_caller_job)
                workflow.add_edge(indel_realign_job.name, haplotype_caller_job.name)
                if ch not in haplotype_caller_jobs:
                    haplotype_caller_jobs[ch] = []
                haplotype_caller_jobs[ch].append(haplotype_caller_job)

        for ch in haplotype_caller_jobs:
            # genotype gvcfs job
            job_name = self._generate_job_name('genotype_gvcfs')
            input_files = []
            for job in haplotype_caller_jobs[ch]:
                workflow.add_edge(job.name, job_name)
                input_files.append(job.files[12])
                input_files.append(job.files[13])
            genotype_gvcfs_job = self._generate_job('genotype_gvcfs', job_name, input_files=input_files)
            workflow.add_node(job_name, job=genotype_gvcfs_job)
            genotype_gvcfs_jobs.append(genotype_gvcfs_job)

        # merge gcvf job
        job_name = self._generate_job_name('merge_gcvf')
        input_files = []
        for ch in haplotype_caller_jobs:
            for job in haplotype_caller_jobs[ch]:
                workflow.add_edge(job.name, job_name)
                input_files.append(job.files[12])
                input_files.append(job.files[13])
        merge_gcvf_job = self._generate_job('merge_gcvf', job_name, input_files=input_files)
        workflow.add_node(job_name, job=merge_gcvf_job)

        # combine variants job
        job_name = self._generate_job_name('combine_variants')
        input_files = []
        for job in genotype_gvcfs_jobs:
            workflow.add_edge(job.name, job_name)
            input_files.extend(self._get_files_by_job_and_link(job.name, FileLink.OUTPUT))
        combine_variants_job = self._generate_job('combine_variants', job_name, input_files=input_files)
        workflow.add_node(job_name, job=combine_variants_job)

        # select variants indel job
        input_files = self._get_files_by_job_and_link(combine_variants_job.name, FileLink.OUTPUT)
        job_name = self._generate_job_name('select_variants_indel')
        svi_job = self._generate_job('select_variants_indel', job_name, input_files=input_files)
        workflow.add_node(job_name, job=svi_job)
        workflow.add_edge(combine_variants_job.name, svi_job.name)

        # select variants snp job
        job_name = self._generate_job_name('select_variants_snp')
        svs_job = self._generate_job('select_variants_snp', job_name, input_files=input_files)
        workflow.add_node(job_name, job=svs_job)
        workflow.add_edge(combine_variants_job.name, svs_job.name)

        # filtering indel job
        job_name = self._generate_job_name('filtering_indel')
        fi_job = self._generate_job('filtering_indel', job_name,
                                    input_files=self._get_files_by_job_and_link(svi_job.name, FileLink.OUTPUT))
        workflow.add_node(job_name, job=fi_job)
        workflow.add_edge(svi_job.name, fi_job.name)

        # filtering indel job
        job_name = self._generate_job_name('filtering_snp')
        fs_job = self._generate_job('filtering_snp', job_name,
                                    input_files=self._get_files_by_job_and_link(svs_job.name, FileLink.OUTPUT))
        workflow.add_node(job_name, job=fs_job)
        workflow.add_edge(svs_job.name, fs_job.name)

        self.workflows.append(workflow)
        return workflow

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the SoyKB workflow. Recipes can be
        generated by using the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`.

        :return: A recipe in the form of a dictionary in which keys are job prefixes.
        :rtype: Dict[str, Any]
        """
        return {
            "alignment_to_reference": {
                "runtime": {
                    "min": 3.305,
                    "max": 133.207,
                    "distribution": {
                        "name": "triang",
                        "params": [
                            0.9999998742580727,
                            -6.69821864227278,
                            61.698354505102785
                        ]
                    }
                },
                "input": {
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".fastq": {
                        "distribution": "None",
                        "min": 23845,
                        "max": 23845
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    },
                    ".txt": {
                        "distribution": {
                            "name": "triang",
                            "params": [
                                0.9999996526041718,
                                -8.777652386868017,
                                83.77775048163616
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".sam": {
                        "distribution": "None",
                        "min": 112561,
                        "max": 112561
                    }
                }
            },
            "sort_sam": {
                "runtime": {
                    "min": 2.082,
                    "max": 89.242,
                    "distribution": {
                        "name": "rayleigh",
                        "params": [
                            -2.2525914610825772,
                            17.64291114779886
                        ]
                    }
                },
                "input": {
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".sam": {
                        "distribution": "None",
                        "min": 112561,
                        "max": 112561
                    }
                },
                "output": {
                    ".bai": {
                        "distribution": "None",
                        "min": 127424,
                        "max": 127424
                    },
                    ".bam": {
                        "distribution": "None",
                        "min": 36576,
                        "max": 36576
                    }
                }
            },
            "dedup": {
                "runtime": {
                    "min": 3.169,
                    "max": 152.024,
                    "distribution": {
                        "name": "triang",
                        "params": [
                            0.9999992247026201,
                            -8.076384057650046,
                            77.07646956305474
                        ]
                    }
                },
                "input": {
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".bai": {
                        "distribution": "None",
                        "min": 127424,
                        "max": 127424
                    },
                    ".bam": {
                        "distribution": "None",
                        "min": 36576,
                        "max": 36576
                    }
                },
                "output": {
                    ".bai": {
                        "distribution": "None",
                        "min": 127424,
                        "max": 127424
                    },
                    ".bam": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.2109632398612863,
                                -2.836686397932103e-26,
                                2.129201183155014
                            ]
                        },
                        "min": 37467,
                        "max": 37470
                    }
                }
            },
            "add_replace": {
                "runtime": {
                    "min": 3.382,
                    "max": 139.801,
                    "distribution": {
                        "name": "argus",
                        "params": [
                            1.1554901234434413e-05,
                            -14.968557591139366,
                            64.32949317987402
                        ]
                    }
                },
                "input": {
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".bai": {
                        "distribution": "None",
                        "min": 127424,
                        "max": 127424
                    },
                    ".bam": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.2109632398612863,
                                -2.836686397932103e-26,
                                2.129201183155014
                            ]
                        },
                        "min": 37467,
                        "max": 37470
                    }
                },
                "output": {
                    ".bam": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5130616667328313,
                                -4.421427043482214e-27,
                                1.9064184774937238
                            ]
                        },
                        "min": 37953,
                        "max": 37956
                    },
                    ".bai": {
                        "distribution": "None",
                        "min": 127424,
                        "max": 127424
                    }
                }
            },
            "realign_target_creator": {
                "runtime": {
                    "min": 78.698,
                    "max": 474.546,
                    "distribution": {
                        "name": "fisk",
                        "params": [
                            0.42030582758525614,
                            -3.368844328154359e-28,
                            2.1958097887228156
                        ]
                    }
                },
                "input": {
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".bai": {
                        "distribution": "None",
                        "min": 127424,
                        "max": 127424
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    },
                    ".bam": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5130616667328313,
                                -4.421427043482214e-27,
                                1.9064184774937238
                            ]
                        },
                        "min": 37953,
                        "max": 37956
                    },
                    ".txt": {
                        "distribution": {
                            "name": "triang",
                            "params": [
                                0.9999996526041718,
                                -8.777652386868017,
                                83.77775048163616
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".list": {
                        "distribution": "None",
                        "min": 40,
                        "max": 40
                    }
                }
            },
            "indel_realign": {
                "runtime": {
                    "min": 5.589,
                    "max": 154.005,
                    "distribution": {
                        "name": "gamma",
                        "params": [
                            0.6031493832274422,
                            4.999999999999999,
                            10.062244059815622
                        ]
                    }
                },
                "input": {
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".bai": {
                        "distribution": "None",
                        "min": 127424,
                        "max": 127424
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".list": {
                        "distribution": "None",
                        "min": 40,
                        "max": 40
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    },
                    ".bam": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5130616667328313,
                                -4.421427043482214e-27,
                                1.9064184774937238
                            ]
                        },
                        "min": 37953,
                        "max": 37956
                    },
                    ".txt": {
                        "distribution": {
                            "name": "triang",
                            "params": [
                                0.9999996526041718,
                                -8.777652386868017,
                                83.77775048163616
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".bai": {
                        "distribution": "None",
                        "min": 127424,
                        "max": 127424
                    },
                    ".bam": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5702069492082148,
                                -1.2181023359846581e-29,
                                2.328292348760873
                            ]
                        },
                        "min": 39550,
                        "max": 39554
                    }
                }
            },
            "haplotype_caller": {
                "runtime": {
                    "min": 30.011,
                    "max": 709.504,
                    "distribution": {
                        "name": "trapz",
                        "params": [
                            1.0,
                            1.0,
                            -14.805000000000001,
                            169.2
                        ]
                    }
                },
                "input": {
                    ".bai": {
                        "distribution": "None",
                        "min": 127424,
                        "max": 127424
                    },
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".bam": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.6414526843394122,
                                -6.597458275218181e-26,
                                1.7270338535457088
                            ]
                        },
                        "min": 39550,
                        "max": 39554
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    },
                    ".txt": {
                        "distribution": {
                            "name": "skewnorm",
                            "params": [
                                810640.4518672728,
                                -0.0010687481525930745,
                                157.76819262166055
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.7384603818061779,
                                -5.18486550777156e-25,
                                2.7633104484562185
                            ]
                        },
                        "min": 59037,
                        "max": 75728
                    },
                    ".idx": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5112505250138766,
                                -7.03556641822723e-26,
                                1.7971380170612603
                            ]
                        },
                        "min": 28447,
                        "max": 28449
                    }
                }
            },
            "merge_gcvf": {
                "runtime": {
                    "min": 2562.381,
                    "max": 37861.356,
                    "distribution": {
                        "name": "levy",
                        "params": [
                            9.99999999999741,
                            5.729323848066606e-13
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.7384603818061779,
                                -5.18486550777156e-25,
                                2.7633104484562185
                            ]
                        },
                        "min": 59037,
                        "max": 75728
                    },
                    ".idx": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5112505250138766,
                                -7.03556641822723e-26,
                                1.7971380170612603
                            ]
                        },
                        "min": 28447,
                        "max": 28449
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    },
                    ".txt": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".list": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 1700,
                        "max": 17000
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    }
                },
                "output": {
                    ".idx": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 28886,
                        "max": 29187
                    },
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 104389,
                        "max": 250300
                    }
                }
            },
            "genotype_gvcfs": {
                "runtime": {
                    "min": 95.618,
                    "max": 899.232,
                    "distribution": {
                        "name": "levy",
                        "params": [
                            -1.237686094683259,
                            5.73897660455604
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.7384603818061779,
                                -5.18486550777156e-25,
                                2.7633104484562185
                            ]
                        },
                        "min": 59037,
                        "max": 75728
                    },
                    ".idx": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5112505250138766,
                                -7.03556641822723e-26,
                                1.7971380170612603
                            ]
                        },
                        "min": 28447,
                        "max": 28449
                    },
                    ".txt": {
                        "distribution": {
                            "name": "uniform",
                            "params": [
                                0.0,
                                100.0
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    },
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.4153135006398344,
                                -1.574105802483078e-26,
                                2.0711234595178167
                            ]
                        },
                        "min": 59101,
                        "max": 83149
                    },
                    ".idx": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.2532229840213113,
                                -5.220068867146425e-26,
                                2.1172225587125273
                            ]
                        },
                        "min": 28496,
                        "max": 28584
                    }
                }
            },
            "combine_variants": {
                "runtime": {
                    "min": 4.156,
                    "max": 163.914,
                    "distribution": {
                        "name": "dweibull",
                        "params": [
                            0.7108618995255725,
                            9.999999999999996,
                            0.052673745569323874
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.4153135006398344,
                                -1.574105802483078e-26,
                                2.0711234595178167
                            ]
                        },
                        "min": 59101,
                        "max": 83149
                    },
                    ".idx": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.2532229840213113,
                                -5.220068867146425e-26,
                                2.1172225587125273
                            ]
                        },
                        "min": 28496,
                        "max": 28584
                    },
                    ".txt": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    },
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 70314,
                        "max": 87692
                    },
                    ".idx": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 28576,
                        "max": 28577
                    }
                }
            },
            "select_variants_snp": {
                "runtime": {
                    "min": 50.953,
                    "max": 197.487,
                    "distribution": {
                        "name": "levy",
                        "params": [
                            9.99999999999741,
                            5.729323848066606e-13
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 70314,
                        "max": 87692
                    },
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    },
                    ".txt": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 71866,
                        "max": 86556
                    }
                }
            },
            "filtering_snp": {
                "runtime": {
                    "min": 4.264,
                    "max": 67.509,
                    "distribution": {
                        "name": "dweibull",
                        "params": [
                            0.7108618995255725,
                            9.999999999999996,
                            0.052673745569323874
                        ]
                    }
                },
                "input": {
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    },
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 71866,
                        "max": 86556
                    },
                    ".txt": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 74686,
                        "max": 89376
                    },
                    ".idx": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 28545,
                        "max": 28546
                    }
                }
            },
            "select_variants_indel": {
                "runtime": {
                    "min": 48.59,
                    "max": 207.667,
                    "distribution": {
                        "name": "levy",
                        "params": [
                            9.99999999999741,
                            5.729323848066606e-13
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 70314,
                        "max": 87692
                    },
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    },
                    ".txt": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 67140,
                        "max": 72804
                    }
                }
            },
            "filtering_indel": {
                "runtime": {
                    "min": 4.028,
                    "max": 338.588,
                    "distribution": {
                        "name": "dweibull",
                        "params": [
                            0.7108618995255725,
                            9.999999999999996,
                            0.052673745569323874
                        ]
                    }
                },
                "input": {
                    ".dict": {
                        "distribution": "None",
                        "min": 142561,
                        "max": 142561
                    },
                    ".pac": {
                        "distribution": "None",
                        "min": 244623820,
                        "max": 244623820
                    },
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 67140,
                        "max": 72804
                    },
                    ".fai": {
                        "distribution": "None",
                        "min": 41379,
                        "max": 41379
                    },
                    ".sa": {
                        "distribution": "None",
                        "min": 489247688,
                        "max": 489247688
                    },
                    ".gz": {
                        "distribution": "None",
                        "min": 108988506,
                        "max": 108988506
                    },
                    ".amb": {
                        "distribution": "None",
                        "min": 259156,
                        "max": 259156
                    },
                    ".fa": {
                        "distribution": "None",
                        "min": 990744229,
                        "max": 990744229
                    },
                    ".ann": {
                        "distribution": "None",
                        "min": 47448,
                        "max": 47448
                    },
                    ".bwt": {
                        "distribution": "None",
                        "min": 978495356,
                        "max": 978495356
                    },
                    ".txt": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 69861,
                        "max": 75525
                    },
                    ".idx": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99999999999741,
                                5.729323848066606e-13
                            ]
                        },
                        "min": 28575,
                        "max": 28576
                    }
                }
            }
        }
