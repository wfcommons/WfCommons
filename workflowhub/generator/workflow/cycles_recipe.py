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
from ...common.workflow import Workflow
from ...utils import ncr


class CyclesRecipe(WorkflowRecipe):
    """A Cycles workflow recipe class for creating synthetic workflow traces.

    :param num_points: The number of points of the spatial grid cell.
    :type num_points: int
    :param num_crops: The number of crops being evaluated.
    :type num_crops: int
    :param num_params: The number of parameter values from the simulation matrix.
    :type num_params: int
    :param data_footprint: The upper bound for the workflow total data footprint (in bytes).
    :type data_footprint: int
    :param num_tasks: The upper bound for the total number of tasks in the workflow.
    :type num_tasks: int
    """

    def __init__(self,
                 num_points: Optional[int] = 1,
                 num_crops: Optional[int] = 1,
                 num_params: Optional[int] = 4,
                 data_footprint: Optional[int] = 0,
                 num_tasks: Optional[int] = 7
                 ) -> None:
        """Create an object of the Cycles workflow recipe."""
        super().__init__("Cycles", data_footprint, num_tasks)

        self.num_points: int = num_points
        self.num_crops: int = num_crops
        self.num_params: int = num_params

    @classmethod
    def from_num_tasks(cls, num_tasks: int) -> 'CyclesRecipe':
        """
        Instantiate a Cycles workflow recipe that will generate synthetic workflows up to
        the total number of tasks provided.

        :param num_tasks: The upper bound for the total number of tasks in the workflow (at least 7).
        :type num_tasks: int

        :return: A Cycles workflow recipe object that will generate synthetic workflows up
                 to the total number of tasks provided.
        :rtype: CyclesRecipe
        """
        if num_tasks < 7:
            raise ValueError("The upper bound for the number of tasks should be at least 7.")

        num_points = 1
        num_crops = 1
        num_params = 4
        remaining_tasks = num_tasks - 7

        while remaining_tasks > 0:
            added_task = False
            cost_param = (ncr(num_params + 1, 4) - ncr(num_params, 4)) * 4 * num_crops * num_points
            if remaining_tasks >= cost_param:
                num_params += 1
                remaining_tasks -= cost_param
                added_task = True

            cost_crop = ncr(num_params, 4) * 4 * num_points + 3
            if remaining_tasks >= cost_crop:
                num_crops += 1
                remaining_tasks -= cost_crop
                added_task = True

            cost_point = ncr(num_params, 4) * 4 * num_crops + 3 * num_crops
            if remaining_tasks >= cost_point and remaining_tasks >= num_tasks / 2:
                num_points += 1
                remaining_tasks -= cost_point
                added_task = True

            if not added_task:
                break

        return cls(num_points=num_points, num_crops=num_crops, num_params=num_params, data_footprint=None,
                   num_tasks=num_tasks)

    @classmethod
    def from_points_and_crops(cls,
                              num_points: int,
                              num_crops: int,
                              num_params: int,
                              ) -> 'CyclesRecipe':
        """
        Instantiate a Cycles workflow recipe that will generate synthetic workflows using
        the defined number of points, crops, and params.

        :param num_points: The number of points of the spatial grid cell.
        :type num_points: int
        :param num_crops: The number of crops being evaluated.
        :type num_crops: int
        :param num_params: The number of parameter values from the simulation matrix.
        :type num_params: int

        :return: A Cycles workflow recipe object that will generate synthetic workflows
                 using the defined number of points, crops, and params.
        :rtype: CyclesRecipe
        """
        if num_points < 1:
            raise ValueError("The number of points should be 1 or higher.")
        if num_crops < 1:
            raise ValueError("The number of crops should be 1 or higher.")
        if num_params < 4:
            raise ValueError("The number of params should be 4 or higher.")

        return cls(num_points=num_points, num_crops=num_crops, num_params=num_params, data_footprint=None,
                   num_tasks=None)

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """Generate a synthetic workflow trace of a Cycles workflow.

        :param workflow_name: The workflow name
        :type workflow_name: int

        :return: A synthetic workflow trace object.
        :rtype: Workflow
        """
        workflow = Workflow(name=self.name + "-synthetic-trace" if not workflow_name else workflow_name, makespan=None)
        self.task_id_counter: int = 1
        num_combinations = ncr(self.num_params, 4)
        summary_tasks_per_crop = {}

        for _ in range(0, self.num_points):
            for crop in range(0, self.num_crops):
                if crop not in summary_tasks_per_crop:
                    summary_tasks_per_crop[crop] = []

                cycles_tasks = []
                cycles_fi_output_tasks = []

                for _ in range(0, num_combinations):
                    # baseline cycles task
                    task_name = self._generate_task_name("baseline_cycles")
                    baseline_cycles_task = self._generate_task('baseline_cycles', task_name)
                    workflow.add_node(task_name, task=baseline_cycles_task)
                    input_files = self._get_files_by_task_and_link(baseline_cycles_task.name, FileLink.OUTPUT)

                    # cycles task
                    task_name = self._generate_task_name("cycles")
                    cycles_task = self._generate_task('cycles', task_name, input_files)
                    workflow.add_node(task_name, task=cycles_task)
                    cycles_tasks.append(cycles_task)

                    # fertilizer increase cycles task
                    task_name = self._generate_task_name("fertilizer_increase_cycles")
                    fi_cycles_task = self._generate_task('fertilizer_increase_cycles', task_name, input_files)
                    workflow.add_node(task_name, task=fi_cycles_task)

                    # fertilizer increase output parser cycles task
                    input_files = self._get_files_by_task_and_link(fi_cycles_task.name, FileLink.OUTPUT)
                    task_name = self._generate_task_name("cycles_fertilizer_increase_output_parser")
                    cycles_fi_output_task = self._generate_task('cycles_fertilizer_increase_output_parser', task_name,
                                                                input_files)
                    workflow.add_node(task_name, task=cycles_fi_output_task)
                    cycles_fi_output_tasks.append(cycles_fi_output_task)

                    # add dependencies
                    workflow.add_edge(baseline_cycles_task.name, cycles_task.name)
                    workflow.add_edge(baseline_cycles_task.name, fi_cycles_task.name)
                    workflow.add_edge(fi_cycles_task.name, cycles_fi_output_task.name)

                # cycles output summary
                input_files = []
                for j in cycles_tasks:
                    input_files.extend(self._get_files_by_task_and_link(j.name, FileLink.OUTPUT))
                task_name = self._generate_task_name("cycles_output_summary")
                cycles_output_summary_task = self._generate_task('cycles_output_summary', task_name, input_files)
                workflow.add_node(task_name, task=cycles_output_summary_task)
                for j in cycles_tasks:
                    workflow.add_edge(j.name, cycles_output_summary_task.name)
                summary_tasks_per_crop[crop].append(cycles_output_summary_task)

                # cycles fertilizer increase output summary
                input_files = []
                for j in cycles_fi_output_tasks:
                    input_files.extend(self._get_files_by_task_and_link(j.name, FileLink.OUTPUT))
                task_name = self._generate_task_name("cycles_fertilizer_increase_output_summary")
                cycles_fi_output_summary_task = self._generate_task('cycles_fertilizer_increase_output_summary',
                                                                    task_name, input_files)
                workflow.add_node(task_name, task=cycles_fi_output_summary_task)
                for j in cycles_fi_output_tasks:
                    workflow.add_edge(j.name, cycles_fi_output_summary_task.name)

        # cycles plots
        for crop in summary_tasks_per_crop:
            input_files = []
            for j in summary_tasks_per_crop[crop]:
                input_files.extend(self._get_files_by_task_and_link(j.name, FileLink.OUTPUT))
            task_name = self._generate_task_name("cycles_plots")
            cycles_plots_task = self._generate_task('cycles_plots', task_name, input_files)
            workflow.add_node(task_name, task=cycles_plots_task)
            for j in summary_tasks_per_crop[crop]:
                workflow.add_edge(j.name, cycles_plots_task.name)

        self.workflows.append(workflow)
        return workflow

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the Cycles workflow. Recipes can be
        generated by using the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`.

        :return: A recipe in the form of a dictionary in which keys are task prefixes.
        :rtype: Dict[str, Any]
        """
        return {
            "baseline_cycles": {
                "runtime": {
                    "min": 1.597,
                    "max": 110.598,
                    "distribution": {
                        "name": "alpha",
                        "params": [
                            8.047751262283954e-09,
                            -1.2756014938873586,
                            2.6812435906886023
                        ]
                    }
                },
                "input": {
                    ".soil": {
                        "distribution": "None",
                        "min": 739,
                        "max": 739
                    },
                    ".weather": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.48810346862081644,
                                -2.2156889517666317e-25,
                                2.1061107035081923
                            ]
                        },
                        "min": 480096,
                        "max": 480097
                    },
                    ".operation": {
                        "distribution": {
                            "name": "chi2",
                            "params": [
                                1.5264675653031392,
                                -9.478040328613063e-27,
                                335.0447820745495
                            ]
                        },
                        "min": 673,
                        "max": 1806
                    },
                    ".ctrl": {
                        "distribution": "None",
                        "min": 766,
                        "max": 766
                    },
                    "cycles_exe": {
                        "distribution": "None",
                        "min": 694648,
                        "max": 694648
                    },
                    ".crop": {
                        "distribution": "None",
                        "min": 14434,
                        "max": 14434
                    }
                },
                "output": {
                    ".dat": {
                        "distribution": {
                            "name": "triang",
                            "params": [
                                3.9109273745053685e-13,
                                -60.360981337099844,
                                7636.433793237831
                            ]
                        },
                        "min": 1056,
                        "max": 1808675
                    },
                    ".csv": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.9882276133785874,
                                279.0769359502969,
                                279.07693595029696
                            ]
                        },
                        "min": 275,
                        "max": 281
                    },
                    ".zip": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.4312876961814035,
                                -4.30712563367272e-26,
                                2.0427130739715498
                            ]
                        },
                        "min": 2267492,
                        "max": 4458415
                    }
                }
            },
            "cycles": {
                "runtime": {
                    "min": 1.277,
                    "max": 101.596,
                    "distribution": {
                        "name": "alpha",
                        "params": [
                            8.484516484670589e-09,
                            -1.2578651677609312,
                            2.4675320024710006
                        ]
                    }
                },
                "input": {
                    ".soil": {
                        "distribution": "None",
                        "min": 739,
                        "max": 739
                    },
                    ".weather": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.48810346862081644,
                                -2.2156889517666317e-25,
                                2.1061107035081923
                            ]
                        },
                        "min": 480096,
                        "max": 480097
                    },
                    ".operation": {
                        "distribution": {
                            "name": "chi2",
                            "params": [
                                1.5264675653031392,
                                -9.478040328613063e-27,
                                335.0447820745495
                            ]
                        },
                        "min": 673,
                        "max": 1806
                    },
                    ".dat": {
                        "distribution": "None",
                        "min": 37728,
                        "max": 37728
                    },
                    ".ctrl": {
                        "distribution": "None",
                        "min": 766,
                        "max": 766
                    },
                    "cycles_exe": {
                        "distribution": "None",
                        "min": 694648,
                        "max": 694648
                    },
                    ".crop": {
                        "distribution": "None",
                        "min": 14434,
                        "max": 14434
                    }
                },
                "output": {
                    ".zip": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.48427962013577297,
                                -1.99155154789682e-28,
                                2.0977079850885403
                            ]
                        },
                        "min": 2928100,
                        "max": 4476265
                    },
                    ".dat": {
                        "distribution": {
                            "name": "triang",
                            "params": [
                                3.9109273745053685e-13,
                                -60.360981337099844,
                                7636.433793237831
                            ]
                        },
                        "min": 1056,
                        "max": 1808675
                    },
                    ".csv": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.9882276133785874,
                                279.0769359502969,
                                279.07693595029696
                            ]
                        },
                        "min": 275,
                        "max": 281
                    }
                }
            },
            "fertilizer_increase_cycles": {
                "runtime": {
                    "min": 1.258,
                    "max": 101.229,
                    "distribution": {
                        "name": "alpha",
                        "params": [
                            1.0434096097164159e-08,
                            -1.2790233222623084,
                            2.4991241569098097
                        ]
                    }
                },
                "input": {
                    ".soil": {
                        "distribution": "None",
                        "min": 739,
                        "max": 739
                    },
                    ".weather": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.48810346862081644,
                                -2.2156889517666317e-25,
                                2.1061107035081923
                            ]
                        },
                        "min": 480096,
                        "max": 480097
                    },
                    ".operation": {
                        "distribution": {
                            "name": "chi2",
                            "params": [
                                1.5264675653031392,
                                -9.478040328613063e-27,
                                335.0447820745495
                            ]
                        },
                        "min": 673,
                        "max": 1806
                    },
                    ".dat": {
                        "distribution": "None",
                        "min": 37728,
                        "max": 37728
                    },
                    ".ctrl": {
                        "distribution": "None",
                        "min": 766,
                        "max": 766
                    },
                    "cycles_exe": {
                        "distribution": "None",
                        "min": 694648,
                        "max": 694648
                    },
                    ".crop": {
                        "distribution": "None",
                        "min": 14434,
                        "max": 14434
                    }
                },
                "output": {
                    ".csv": {
                        "distribution": {
                            "name": "cosine",
                            "params": [
                                83.23949362619993,
                                152.0767648343908
                            ]
                        },
                        "min": 274,
                        "max": 282
                    },
                    ".dat": {
                        "distribution": {
                            "name": "triang",
                            "params": [
                                3.9109273745053685e-13,
                                -60.360981337099844,
                                7636.433793237831
                            ]
                        },
                        "min": 1056,
                        "max": 1808675
                    },
                    ".zip": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.33716386282358723,
                                -1.0674750103622264e-27,
                                2.1701869928746937
                            ]
                        },
                        "min": 2928677,
                        "max": 4486752
                    }
                }
            },
            "cycles_fertilizer_increase_output_parser": {
                "runtime": {
                    "min": 0.035,
                    "max": 1.235,
                    "distribution": {
                        "name": "pareto",
                        "params": [
                            3.178595938028368,
                            -0.5759653856097273,
                            0.5759653845181962
                        ]
                    }
                },
                "input": {
                    ".dat": {
                        "distribution": "None",
                        "min": 5980,
                        "max": 5980
                    },
                    ".csv": {
                        "distribution": {
                            "name": "skewnorm",
                            "params": [
                                1833079.2807219662,
                                -0.0006331045037050451,
                                218.08284593142628
                            ]
                        },
                        "min": 274,
                        "max": 282
                    }
                },
                "output": {
                    ".csv": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5770664167941257,
                                -5.756194313327051e-26,
                                1.6405791319945826
                            ]
                        },
                        "min": 8224,
                        "max": 8604
                    }
                }
            },
            "cycles_output_summary": {
                "runtime": {
                    "min": 0.048,
                    "max": 0.39,
                    "distribution": {
                        "name": "beta",
                        "params": [
                            0.19383300430623274,
                            200.14346240436433,
                            -1.1535111701525267e-29,
                            2746.191287764276
                        ]
                    }
                },
                "input": {
                    ".dat": {
                        "distribution": "None",
                        "min": 5980,
                        "max": 5980
                    },
                    ".csv": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.9882276133785874,
                                279.0769359502969,
                                279.07693595029696
                            ]
                        },
                        "min": 275,
                        "max": 281
                    }
                },
                "output": {
                    ".csv": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.6824174018857059,
                                -1.7507576397743858e-28,
                                1.9958963207829474
                            ]
                        },
                        "min": 74996,
                        "max": 257366
                    }
                }
            },
            "cycles_fertilizer_increase_output_summary": {
                "runtime": {
                    "min": 0.045,
                    "max": 0.538,
                    "distribution": {
                        "name": "chi",
                        "params": [
                            0.15755388349696192,
                            -2.6271332657183503e-27,
                            81.67483343483318
                        ]
                    }
                },
                "input": {
                    ".csv": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5770664167941257,
                                -5.756194313327051e-26,
                                1.6405791319945826
                            ]
                        },
                        "min": 8224,
                        "max": 8604
                    }
                },
                "output": {
                    ".csv": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.6824174018857059,
                                -1.7507576397743858e-28,
                                1.9958963207829474
                            ]
                        },
                        "min": 122964,
                        "max": 422954
                    }
                }
            },
            "cycles_plots": {
                "runtime": {
                    "min": 123.38,
                    "max": 514.566,
                    "distribution": {
                        "name": "fisk",
                        "params": [
                            0.07143329405623092,
                            11.999999999999998,
                            5.486627664253447
                        ]
                    }
                },
                "input": {
                    ".csv": {
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.6824174018857059,
                                -1.7507576397743858e-28,
                                1.9958963207829474
                            ]
                        },
                        "min": 74996,
                        "max": 257366
                    }
                },
                "output": {
                    ".gif": {
                        "distribution": {
                            "name": "levy",
                            "params": [
                                9.99996901653408,
                                9.326668509495445e-05
                            ]
                        },
                        "min": 2126304,
                        "max": 6383929
                    }
                }
            }
        }
