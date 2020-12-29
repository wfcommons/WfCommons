#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WorkflowHub Team.
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
from ...common.task import Task
from ...common.workflow import Workflow


class SoyKBRecipe(WorkflowRecipe):
    """A SoyKB workflow recipe class for creating synthetic workflow traces.

    :param num_fastq_files: The number of FASTQ files to be analyzed.
    :type num_fastq_files: int
    :param num_chromosomes: The number of chromosomes.
    :type num_chromosomes: int
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
                 num_fastq_files: Optional[int] = 2,
                 num_chromosomes: Optional[int] = 1,
                 data_footprint: Optional[int] = 0,
                 num_tasks: Optional[int] = 14,
                 runtime_factor: Optional[float] = 1.0,
                 input_file_size_factor: Optional[float] = 1.0,
                 output_file_size_factor: Optional[float] = 1.0
                 ) -> None:
        """Create an object of the SoyKB workflow recipe."""
        super().__init__("SoyKB",
                         data_footprint,
                         num_tasks,
                         runtime_factor,
                         input_file_size_factor,
                         output_file_size_factor)

        self.num_fastq_files = num_fastq_files
        self.num_chromosomes = num_chromosomes

    @classmethod
    def from_num_tasks(cls,
                       num_tasks: int,
                       runtime_factor: Optional[float] = 1.0,
                       input_file_size_factor: Optional[float] = 1.0,
                       output_file_size_factor: Optional[float] = 1.0
                       ) -> 'SoyKBRecipe':
        """
        Instantiate a SoyKB workflow recipe that will generate synthetic workflows up to
        the total number of tasks provided.

        :param num_tasks: The upper bound for the total number of tasks in the workflow (at least 14).
        :type num_tasks: int
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

        :return: A SoyKB workflow recipe object that will generate synthetic workflows up
                 to the total number of tasks provided.
        :rtype: SoyKBRecipe
        """
        if num_tasks < 14:
            raise ValueError("The upper bound for the number of tasks should be at least 14.")

        num_chromosomes = 1 if num_tasks == 14 else random.randint(1, min(math.ceil((num_tasks - 14) / 2), 22))
        remaining_tasks = num_tasks - (2 * num_chromosomes) - 12
        num_fastq_files = 1

        while remaining_tasks > 0:
            if remaining_tasks >= num_chromosomes + 6:
                num_fastq_files += 1
                remaining_tasks -= num_chromosomes + 6
            else:
                break

        return cls(num_fastq_files=num_fastq_files * 2,
                   num_chromosomes=num_chromosomes,
                   data_footprint=None,
                   num_tasks=num_tasks,
                   runtime_factor=runtime_factor,
                   input_file_size_factor=input_file_size_factor,
                   output_file_size_factor=output_file_size_factor)

    @classmethod
    def from_sequences(cls,
                       num_fastq_files: int,
                       num_chromosomes: int,
                       runtime_factor: Optional[float] = 1.0,
                       input_file_size_factor: Optional[float] = 1.0,
                       output_file_size_factor: Optional[float] = 1.0
                       ) -> 'SoyKBRecipe':
        """
        Instantiate a SoyKB workflow recipe that will generate synthetic workflows using
        the defined number of FASTQ files and chromosomes.

        :param num_fastq_files: The number of FASTQ files to be analyzed (at least 2).
        :type num_fastq_files: int
        :param num_chromosomes: The number of chromosomes (range [1,22].
        :type num_chromosomes: int
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

        :return: A SoyKB workflow recipe object that will generate synthetic workflows
                 using the defined number of FASTQ files and chromosomes.
        :rtype: SoyKBRecipe
        """
        if num_fastq_files < 2 or num_fastq_files % 2 != 0:
            raise ValueError("The number of FASTQ files should be at least 2 and should be an even number.")
        if num_chromosomes < 1 or num_chromosomes > 22:
            raise ValueError("The number of chromosomes should be within range [1, 22].")

        return cls(num_fastq_files=num_fastq_files,
                   num_chromosomes=num_chromosomes,
                   data_footprint=None,
                   num_tasks=None,
                   runtime_factor=runtime_factor,
                   input_file_size_factor=input_file_size_factor,
                   output_file_size_factor=output_file_size_factor)

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """Generate a synthetic workflow trace of a SoyKB workflow.

        :param workflow_name: The workflow name
        :type workflow_name: int

        :return: A synthetic workflow trace object.
        :rtype: Workflow
        """
        workflow = Workflow(name=self.name + "-synthetic-trace" if not workflow_name else workflow_name, makespan=None)
        self.task_id_counter: int = 1
        haplotype_caller_tasks: Dict[int, List[Task]] = {}
        genotype_gvcfs_tasks: List[Task] = []

        num_pipelines: int = math.ceil(self.num_fastq_files / 2)

        for _ in range(0, num_pipelines):
            # align to reference task
            task_name = self._generate_task_name('alignment_to_reference')
            atr_task = self._generate_task('alignment_to_reference', task_name,
                                           files_recipe={FileLink.INPUT: {".fastq": 2}})
            workflow.add_node(task_name, task=atr_task)

            # sort sam task
            task_name = self._generate_task_name('sort_sam')
            sort_sam_task = self._generate_task('sort_sam', task_name, input_files=[atr_task.files[12]])
            workflow.add_node(task_name, task=sort_sam_task)
            workflow.add_edge(atr_task.name, sort_sam_task.name)

            # dedup task
            task_name = self._generate_task_name('dedup')
            dedup_task = self._generate_task('dedup', task_name,
                                             input_files=[sort_sam_task.files[2], sort_sam_task.files[3]])
            workflow.add_node(task_name, task=dedup_task)
            workflow.add_edge(sort_sam_task.name, dedup_task.name)

            # add replace task
            task_name = self._generate_task_name('add_replace')
            add_replace_task = self._generate_task('add_replace', task_name,
                                                   input_files=[dedup_task.files[3], dedup_task.files[4]])
            workflow.add_node(task_name, task=add_replace_task)
            workflow.add_edge(dedup_task.name, add_replace_task.name)

            # realign target creator replace task
            task_name = self._generate_task_name('realign_target_creator')
            rtc_task = self._generate_task('realign_target_creator', task_name,
                                           input_files=[add_replace_task.files[3], add_replace_task.files[4]])
            workflow.add_node(task_name, task=rtc_task)
            workflow.add_edge(add_replace_task.name, rtc_task.name)

            # indel realign task
            task_name = self._generate_task_name('indel_realign')
            indel_realign_task = self._generate_task('indel_realign', task_name, input_files=[rtc_task.files[12]])
            workflow.add_node(task_name, task=indel_realign_task)
            workflow.add_edge(rtc_task.name, indel_realign_task.name)

            for ch in range(0, self.num_chromosomes):
                # haplotype caller task
                task_name = self._generate_task_name('haplotype_caller')
                haplotype_caller_task = self._generate_task('haplotype_caller', task_name,
                                                            input_files=[indel_realign_task.files[13],
                                                                         indel_realign_task.files[14]])
                workflow.add_node(task_name, task=haplotype_caller_task)
                workflow.add_edge(indel_realign_task.name, haplotype_caller_task.name)
                if ch not in haplotype_caller_tasks:
                    haplotype_caller_tasks[ch] = []
                haplotype_caller_tasks[ch].append(haplotype_caller_task)

        for ch in haplotype_caller_tasks:
            # genotype gvcfs task
            task_name = self._generate_task_name('genotype_gvcfs')
            input_files = []
            for task in haplotype_caller_tasks[ch]:
                workflow.add_edge(task.name, task_name)
                input_files.append(task.files[12])
                input_files.append(task.files[13])
            genotype_gvcfs_task = self._generate_task('genotype_gvcfs', task_name, input_files=input_files)
            workflow.add_node(task_name, task=genotype_gvcfs_task)
            genotype_gvcfs_tasks.append(genotype_gvcfs_task)

        # merge gcvf task
        task_name = self._generate_task_name('merge_gcvf')
        input_files = []
        for ch in haplotype_caller_tasks:
            for task in haplotype_caller_tasks[ch]:
                workflow.add_edge(task.name, task_name)
                input_files.append(task.files[12])
                input_files.append(task.files[13])
        merge_gcvf_task = self._generate_task('merge_gcvf', task_name, input_files=input_files)
        workflow.add_node(task_name, task=merge_gcvf_task)

        # combine variants task
        task_name = self._generate_task_name('combine_variants')
        input_files = []
        for task in genotype_gvcfs_tasks:
            workflow.add_edge(task.name, task_name)
            input_files.extend(self._get_files_by_task_and_link(task.name, FileLink.OUTPUT))
        combine_variants_task = self._generate_task('combine_variants', task_name, input_files=input_files)
        workflow.add_node(task_name, task=combine_variants_task)

        # select variants indel task
        input_files = self._get_files_by_task_and_link(combine_variants_task.name, FileLink.OUTPUT)
        task_name = self._generate_task_name('select_variants_indel')
        svi_task = self._generate_task('select_variants_indel', task_name, input_files=input_files)
        workflow.add_node(task_name, task=svi_task)
        workflow.add_edge(combine_variants_task.name, svi_task.name)

        # select variants snp task
        task_name = self._generate_task_name('select_variants_snp')
        svs_task = self._generate_task('select_variants_snp', task_name, input_files=input_files)
        workflow.add_node(task_name, task=svs_task)
        workflow.add_edge(combine_variants_task.name, svs_task.name)

        # filtering indel task
        task_name = self._generate_task_name('filtering_indel')
        fi_task = self._generate_task('filtering_indel', task_name,
                                      input_files=self._get_files_by_task_and_link(svi_task.name, FileLink.OUTPUT))
        workflow.add_node(task_name, task=fi_task)
        workflow.add_edge(svi_task.name, fi_task.name)

        # filtering indel task
        task_name = self._generate_task_name('filtering_snp')
        fs_task = self._generate_task('filtering_snp', task_name,
                                      input_files=self._get_files_by_task_and_link(svs_task.name, FileLink.OUTPUT))
        workflow.add_node(task_name, task=fs_task)
        workflow.add_edge(svs_task.name, fs_task.name)

        self.workflows.append(workflow)
        return workflow

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the SoyKB workflow. Recipes can be
        generated by using the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`.

        :return: A recipe in the form of a dictionary in which keys are task prefixes.
        :rtype: Dict[str, Any]
        """
        return {
            "alignment_to_reference": {
                "runtime": {
                    "min": 3.305,
                    "max": 133.207,
                    "distribution": {
                        "name": "rdist",
                        "params": [
                            1.0497506699055814,
                            0.8751735665957343,
                            0.5778762692984372
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
                            "name": "rdist",
                            "params": [
                                1.7477959990252214,
                                1.0,
                                3.9893803094996296e-22
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
                        "name": "arcsine",
                        "params": [
                            -0.09850603366048286,
                            1.098506033660483
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
                        "name": "rdist",
                        "params": [
                            1.056503552352268,
                            0.8880507079714015,
                            0.5907534106741043
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
                            "name": "rdist",
                            "params": [
                                1.7504259269384672,
                                0.9999999999938114,
                                6.18860928332036e-12
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
                        "name": "arcsine",
                        "params": [
                            -0.09850603366048286,
                            1.098506033660483
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
                            "name": "rdist",
                            "params": [
                                1.7504259269384672,
                                0.9999999999938114,
                                6.18860928332036e-12
                            ]
                        },
                        "min": 37467,
                        "max": 37470
                    }
                },
                "output": {
                    ".bam": {
                        "distribution": {
                            "name": "beta",
                            "params": [
                                0.028007929445023696,
                                0.16801064212397998,
                                -0.08122143521923508,
                                1.0812214352192353
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
                        "name": "arcsine",
                        "params": [
                            -0.08648046099196244,
                            1.0864804609919627
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
                            "name": "beta",
                            "params": [
                                0.028007929445023696,
                                0.16801064212397998,
                                -0.08122143521923508,
                                1.0812214352192353
                            ]
                        },
                        "min": 37953,
                        "max": 37956
                    },
                    ".txt": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.7477959990252214,
                                1.0,
                                3.9893803094996296e-22
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
                        "name": "rdist",
                        "params": [
                            1.054143443775926,
                            0.8779716509043156,
                            0.6087408816735466
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
                            "name": "beta",
                            "params": [
                                0.028007929445023696,
                                0.16801064212397998,
                                -0.08122143521923508,
                                1.0812214352192353
                            ]
                        },
                        "min": 37953,
                        "max": 37956
                    },
                    ".txt": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.7477959990252214,
                                1.0,
                                3.9893803094996296e-22
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
                            "name": "rdist",
                            "params": [
                                1.0354552049234393,
                                0.9230820838539932,
                                0.5764154171873266
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
                        "name": "argus",
                        "params": [
                            0.0008690793369530592,
                            -0.6114566529512571,
                            1.6139824159282967
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
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
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
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
                            ]
                        },
                        "min": 59037,
                        "max": 75728
                    },
                    ".idx": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
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
                        "name": "dgamma",
                        "params": [
                            0.19307932052443602,
                            0.9999999999999999,
                            1.3015975328465914
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
                            ]
                        },
                        "min": 59037,
                        "max": 75728
                    },
                    ".idx": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
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
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 28886,
                        "max": 29187
                    },
                    ".vcf": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                        "name": "arcsine",
                        "params": [
                            -0.09850603366048286,
                            1.098506033660483
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
                            ]
                        },
                        "min": 59037,
                        "max": 75728
                    },
                    ".idx": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
                            ]
                        },
                        "min": 28447,
                        "max": 28449
                    },
                    ".txt": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
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
                            "name": "rdist",
                            "params": [
                                1.7504259269384672,
                                0.9999999999938114,
                                6.18860928332036e-12
                            ]
                        },
                        "min": 59101,
                        "max": 83149
                    },
                    ".idx": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.7504259269384672,
                                0.9999999999938114,
                                6.18860928332036e-12
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
                        "name": "dgamma",
                        "params": [
                            0.19307932052443602,
                            0.9999999999999999,
                            1.3015975328465914
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.7504259269384672,
                                0.9999999999938114,
                                6.18860928332036e-12
                            ]
                        },
                        "min": 59101,
                        "max": 83149
                    },
                    ".idx": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.7504259269384672,
                                0.9999999999938114,
                                6.18860928332036e-12
                            ]
                        },
                        "min": 28496,
                        "max": 28584
                    },
                    ".txt": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 70314,
                        "max": 87692
                    },
                    ".idx": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                        "name": "dgamma",
                        "params": [
                            0.19307932052443602,
                            0.9999999999999999,
                            1.3015975328465914
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                        "name": "dgamma",
                        "params": [
                            0.19307932052443602,
                            0.9999999999999999,
                            1.3015975328465914
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
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 71866,
                        "max": 86556
                    },
                    ".txt": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 74686,
                        "max": 89376
                    },
                    ".idx": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                        "name": "dgamma",
                        "params": [
                            0.19307932052443602,
                            0.9999999999999999,
                            1.3015975328465914
                        ]
                    }
                },
                "input": {
                    ".vcf": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                        "name": "dgamma",
                        "params": [
                            0.19307932052443602,
                            0.9999999999999999,
                            1.3015975328465914
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
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
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
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 60,
                        "max": 15976
                    }
                },
                "output": {
                    ".vcf": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 69861,
                        "max": 75525
                    },
                    ".idx": {
                        "distribution": {
                            "name": "dgamma",
                            "params": [
                                0.19307932052443602,
                                0.9999999999999999,
                                1.3015975328465914
                            ]
                        },
                        "min": 28575,
                        "max": 28576
                    }
                }
            }
        }
