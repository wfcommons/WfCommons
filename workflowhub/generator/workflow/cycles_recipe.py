#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WorkflowHub Team.
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
    :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
    :type runtime_factor: float
    :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
    :type input_file_size_factor: float
    :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
    :type output_file_size_factor: float
    """

    def __init__(self,
                 num_points: Optional[int] = 1,
                 num_crops: Optional[int] = 1,
                 num_params: Optional[int] = 4,
                 data_footprint: Optional[int] = 0,
                 num_tasks: Optional[int] = 7,
                 runtime_factor: Optional[float] = 1.0,
                 input_file_size_factor: Optional[float] = 1.0,
                 output_file_size_factor: Optional[float] = 1.0
                 ) -> None:
        """Create an object of the Cycles workflow recipe."""
        super().__init__("Cycles",
                         data_footprint,
                         num_tasks,
                         runtime_factor,
                         input_file_size_factor,
                         output_file_size_factor)

        self.num_points: int = num_points
        self.num_crops: int = num_crops
        self.num_params: int = num_params

    @classmethod
    def from_num_tasks(cls,
                       num_tasks: int,
                       runtime_factor: Optional[float] = 1.0,
                       input_file_size_factor: Optional[float] = 1.0,
                       output_file_size_factor: Optional[float] = 1.0
                       ) -> 'CyclesRecipe':
        """
        Instantiate a Cycles workflow recipe that will generate synthetic workflows up to
        the total number of tasks provided.

        :param num_tasks: The upper bound for the total number of tasks in the workflow (at least 7).
        :type num_tasks: int
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

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

        return cls(num_points=num_points,
                   num_crops=num_crops,
                   num_params=num_params,
                   data_footprint=None,
                   num_tasks=num_tasks,
                   runtime_factor=runtime_factor,
                   input_file_size_factor=input_file_size_factor,
                   output_file_size_factor=output_file_size_factor)

    @classmethod
    def from_points_and_crops(cls,
                              num_points: int,
                              num_crops: int,
                              num_params: int,
                              runtime_factor: Optional[float] = 1.0,
                              input_file_size_factor: Optional[float] = 1.0,
                              output_file_size_factor: Optional[float] = 1.0
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
        :param runtime_factor: The factor of which tasks runtime will be increased/decreased.
        :type runtime_factor: float
        :param input_file_size_factor: The factor of which tasks input files size will be increased/decreased.
        :type input_file_size_factor: float
        :param output_file_size_factor: The factor of which tasks output files size will be increased/decreased.
        :type output_file_size_factor: float

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

        return cls(num_points=num_points,
                   num_crops=num_crops,
                   num_params=num_params,
                   data_footprint=None,
                   num_tasks=None,
                   runtime_factor=runtime_factor,
                   input_file_size_factor=input_file_size_factor,
                   output_file_size_factor=output_file_size_factor)

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
                        "name": "argus",
                        "params": [
                            5.734652014997823e-05,
                            -0.6805728508716575,
                            1.6811901314620572
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
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
                            ]
                        },
                        "min": 480096,
                        "max": 480097
                    },
                    ".operation": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
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
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
                            ]
                        },
                        "min": 1056,
                        "max": 1808675
                    },
                    ".csv": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 275,
                        "max": 281
                    },
                    ".zip": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                4.427079012898356e-05,
                                -0.7120439234477358,
                                1.712710243275252
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
                        "name": "argus",
                        "params": [
                            2.4601472008810825e-05,
                            -0.760671189673141,
                            1.761325389492184
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
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
                            ]
                        },
                        "min": 480096,
                        "max": 480097
                    },
                    ".operation": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
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
                            "name": "argus",
                            "params": [
                                3.76005481127054e-05,
                                -0.6638469786758867,
                                1.6644768181414857
                            ]
                        },
                        "min": 2928100,
                        "max": 4476265
                    },
                    ".dat": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
                            ]
                        },
                        "min": 1056,
                        "max": 1808675
                    },
                    ".csv": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
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
                        "name": "argus",
                        "params": [
                            0.00024090845979553212,
                            -0.7571715014318916,
                            1.757837187249918
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
                            "name": "trapz",
                            "params": [
                                0.9999999999999999,
                                1.0,
                                -0.10500000000000001,
                                1.1999999999999997
                            ]
                        },
                        "min": 480096,
                        "max": 480097
                    },
                    ".operation": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
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
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 274,
                        "max": 282
                    },
                    ".dat": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.750426650564894,
                                0.9999999999999998,
                                2.3624681987884834e-16
                            ]
                        },
                        "min": 1056,
                        "max": 1808675
                    },
                    ".zip": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                5.012766959293754e-05,
                                -0.6684439199051191,
                                1.6690337988564448
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
                    ".dat": {
                        "distribution": "None",
                        "min": 5980,
                        "max": 5980
                    },
                    ".csv": {
                        "distribution": {
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 274,
                        "max": 282
                    }
                },
                "output": {
                    ".csv": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                0.0002218467867105605,
                                -0.7480684549751317,
                                1.7487201180830585
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
                        "name": "norm",
                        "params": [
                            0.08688656476267097,
                            0.2572832376513094
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
                            "name": "trapz",
                            "params": [
                                1.0,
                                1.0,
                                -0.10500000000000001,
                                1.2
                            ]
                        },
                        "min": 275,
                        "max": 281
                    }
                },
                "output": {
                    ".csv": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.7504259613885509,
                                1.0,
                                4.589985375327013e-23
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
                        "name": "argus",
                        "params": [
                            0.0001699057370632583,
                            -0.703883277161877,
                            1.7272982844402591
                        ]
                    }
                },
                "input": {
                    ".csv": {
                        "distribution": {
                            "name": "argus",
                            "params": [
                                0.0002218467867105605,
                                -0.7480684549751317,
                                1.7487201180830585
                            ]
                        },
                        "min": 8224,
                        "max": 8604
                    }
                },
                "output": {
                    ".csv": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.7504259613885509,
                                1.0,
                                4.589985375327013e-23
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
                        "name": "beta",
                        "params": [
                            0.11165228211118816,
                            0.32785735039410685,
                            -0.09952333240693918,
                            1.0995233324069393
                        ]
                    }
                },
                "input": {
                    ".csv": {
                        "distribution": {
                            "name": "rdist",
                            "params": [
                                1.7504259613885509,
                                1.0,
                                4.589985375327013e-23
                            ]
                        },
                        "min": 74996,
                        "max": 257366
                    }
                },
                "output": {
                    ".gif": {
                        "distribution": {
                            "name": "norm",
                            "params": [
                                0.35,
                                0.36784847423905404
                            ]
                        },
                        "min": 2126304,
                        "max": 6383929
                    }
                }
            }
        }
