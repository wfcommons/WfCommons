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
from ...utils import NoValue


class MontageDataset(NoValue):
    """An enumeration of Montage datasets."""
    TWOMASS = '2mass'
    DSS = 'dss'


class _MontagetaskRatios:
    """An auxiliary class for generating Montage tasks."""
    tasks_ratios = {
        # (min, increase_rate, stddev)
        MontageDataset.TWOMASS: {
            'mProject': (68, 44, 21),
            'mDiffFit': (414, 112, 52),
            'mBackground': (68, 23, 4)
        },
        MontageDataset.DSS: {
            'mProject': (4, 4, 4),
            'mDiffFit': (120, 134, 118),
            'mBackground': (4, 4, 4)
        }
    }

    def _get_num_tasks(self, task_name: str, degree: float, dataset: MontageDataset) -> int:
        """Get a random number of tasks to be generated for a task prefix and workflow degree.

        :param task_name: The task name prefix.
        :type task_name: str
        :param degree: The size (in degrees) to be used for the width/height of the final mosaic.
        :type degree: float
        :param dataset: The dataset to use for the mosaic (e.g., 2mass, dss).
        :type dataset: MontageDataset

        :return: A random number of tasks to be generated for a task prefix and workflow degree.
        :rtype: int
        """
        task_recipe = self.tasks_ratios[dataset][task_name]
        factor = math.ceil((degree - 0.5) * 10)
        return int(
            task_recipe[0] + random.randint(task_recipe[1] - task_recipe[2], task_recipe[1] + task_recipe[2]) * factor)

    def _get_max_rate_increase(self, task_name: str, dataset: MontageDataset) -> int:
        """Get the maximum rate of increase for a task prefix by increasing the workflow degree.

        :param task_name: The task name prefix.
        :type task_name: str
        :param dataset: The dataset to use for the mosaic (e.g., 2mass, dss).
        :type dataset: MontageDataset

        :return: The maximum rate of increase for a task prefix by increasing the workflow degree.
        :rtype: int
        """
        task_recipe = self.tasks_ratios[dataset][task_name]
        return task_recipe[1] + task_recipe[2]

    def _get_max_num_tasks(self, task_name: str, degree: float, dataset: MontageDataset) -> int:
        """Get the maximum number of tasks that can be generated for a defined task.

        :param task_name: The task name prefix.
        :type task_name: str
        :param degree: The size (in degrees) to be used for the width/height of the final mosaic.
        :type degree: float
        :param dataset: The dataset to use for the mosaic (e.g., 2mass, dss).
        :type dataset: MontageDataset

        :return: The maximum number of tasks that can be generated for a defined task.
        :rtype: int
        """
        task_recipe = self.tasks_ratios[dataset][task_name]
        factor = math.ceil((degree - 0.5) * 10)
        return task_recipe[0] + (task_recipe[1] + task_recipe[2]) * factor


class MontageRecipe(WorkflowRecipe, _MontagetaskRatios):
    """
    A Montage workflow recipe class for creating synthetic workflow traces. In this
    workflow recipe, traces will follow different recipes for different
    :class:`~workflowhub.generator.workflow.montage_recipe.MontageDataset`.

    :param dataset: The dataset to use for the mosaic (e.g., 2mass, dss).
    :type dataset: MontageDataset
    :param num_bands: The number of bands (e.g., red, blue, and green) used by the workflow.
    :type num_bands: int
    :param degree: The size (in degrees) to be used for the width/height of the final mosaic.
    :type degree: float
    :param data_footprint: The upper bound for the workflow total data footprint (in bytes).
    :type data_footprint: int
    :param num_tasks: The upper bound for the total number of tasks in the workflow.
    :type num_tasks: int
    """

    def __init__(self,
                 dataset: Optional[MontageDataset] = MontageDataset.DSS,
                 num_bands: Optional[int] = 1,
                 degree: Optional[float] = 0.5,
                 data_footprint: Optional[int] = 0,
                 num_tasks: Optional[int] = 133
                 ) -> None:
        """Create an object of the Montage workflow recipe."""
        super().__init__("Montage", data_footprint, num_tasks)

        if not isinstance(dataset, MontageDataset):
            raise TypeError("the dataset should be an instance of MontageDataset object.")

        self.dataset: MontageDataset = dataset
        self.num_bands: Optional[int] = num_bands
        self.degree: Optional[float] = float(format(degree, '.1f'))

    @classmethod
    def from_num_tasks(cls, num_tasks: int) -> 'MontageRecipe':
        """
        Instantiate a Montage workflow recipe that will generate synthetic workflows up to
        the total number of tasks provided.

        :param num_tasks: The upper bound for the total number of tasks in the workflow (at least 133).
        :type num_tasks: int

        :return: A Montage workflow recipe object that will generate synthetic workflows up
                 to the total number of tasks provided.
        :rtype: MontageRecipe
        """
        if num_tasks < 133:
            raise ValueError("The upper bound for the number of tasks should be at least 133.")

        dataset = MontageDataset.DSS if num_tasks < 555 else random.choice(list(MontageDataset))
        base_num_tasks = 133 if dataset == MontageDataset.DSS else 555
        num_bands = 1
        degree = 0.5
        remaining_tasks = num_tasks - base_num_tasks

        while remaining_tasks > 0:
            added_task = False
            cost_degree = (cls._get_max_rate_increase(cls, 'mProject', dataset) * 2
                           + cls._get_max_rate_increase(cls, 'mDiffFit', dataset) + 5) * num_bands
            if remaining_tasks >= cost_degree:
                degree += 0.1
                remaining_tasks -= cost_degree
                added_task = True

            cost_band = cls._get_max_num_tasks(cls, 'mProject', degree, dataset) * 2 \
                        + cls._get_max_num_tasks(cls, 'mDiffFit', degree, dataset) + 5
            if num_bands < 3 and cost_band <= remaining_tasks / 2:
                num_bands += 1
                remaining_tasks -= cost_band
                added_task = True

            if not added_task:
                break

        return cls(dataset=dataset, num_bands=num_bands, degree=degree, data_footprint=None, num_tasks=num_tasks)

    @classmethod
    def from_degree(cls, dataset: MontageDataset, num_bands: int, degree: float) -> 'MontageRecipe':
        """
        Instantiate a Montage workflow recipe that will generate synthetic workflows using
        the defined dataset, number of bands, and degree.

        :param dataset: The dataset to use for the mosaic (e.g., 2mass, dss).
        :type dataset: MontageDataset
        :param num_bands: The number of bands (e.g., red, blue, and green) used by the workflow (at least 1).
        :type num_bands: int
        :param degree: The size (in degrees) to be used for the width/height of the final mosaic (at least 0.5).
        :type degree: float

        :return: A Montage workflow recipe object that will generate synthetic workflows
                 using the defined dataset, number of bands, and degree.
        :rtype: MontageRecipe
        """
        if num_bands < 1:
            raise ValueError("The number of mosaic bands should be at least 1.")
        if degree < 0.5:
            raise ValueError("The workflow degree should be at least 0.5.")

        return cls(dataset=dataset, num_bands=num_bands, degree=degree, data_footprint=None, num_tasks=None)

    def build_workflow(self, workflow_name: str = None) -> Workflow:
        """Generate a synthetic workflow trace of a Montage workflow.

        :param workflow_name: The workflow name
        :type workflow_name: int

        :return: A synthetic workflow trace object.
        :rtype: Workflow
        """
        workflow = Workflow(
            name=self.name + "-" + str(self.dataset.value) + "-synthetic-trace" if not workflow_name else workflow_name,
            makespan=None)
        self.task_id_counter: int = 1
        madd_tasks = []

        for _ in range(0, self.num_bands):
            # mProject tasks
            mproject_tasks = []
            for _ in range(0, self._get_num_tasks('mProject', self.degree, self.dataset)):
                task_name = self._generate_task_name('mProject')
                mproject_task = self._generate_task('mProject', task_name,
                                                    files_recipe={FileLink.OUTPUT: {'.fits': 2}})
                workflow.add_node(task_name, task=mproject_task)
                mproject_tasks.append(mproject_task)

            # mDiffFit tasks
            count = 0
            mdiff_tasks = []
            for _ in range(0, self._get_num_tasks('mDiffFit', self.degree, self.dataset)):
                input_files = []
                task_name = self._generate_task_name('mDiffFit')
                if count < len(mproject_tasks) - 1:
                    input_files.extend(self._get_files_by_task_and_link(mproject_tasks[count].name, FileLink.OUTPUT))
                    input_files.extend(
                        self._get_files_by_task_and_link(mproject_tasks[count + 1].name, FileLink.OUTPUT))
                    workflow.add_edge(mproject_tasks[count].name, task_name)
                    workflow.add_edge(mproject_tasks[count + 1].name, task_name)
                    count += 1
                else:
                    previous_index = -1
                    for _ in range(0, 2):
                        while True:
                            index = random.randint(0, len(mproject_tasks) - 1)
                            if index != previous_index:
                                break
                        previous_index = index
                        input_files.extend(
                            self._get_files_by_task_and_link(mproject_tasks[index].name, FileLink.OUTPUT))
                        workflow.add_edge(mproject_tasks[index].name, task_name)
                mdiff_task = self._generate_task('mDiffFit', task_name, input_files=input_files)
                workflow.add_node(task_name, task=mdiff_task)
                mdiff_tasks.append(mdiff_task)

            # mConcatFit task
            input_files = []
            task_name = self._generate_task_name('mConcatFit')
            for mdiff_task in mdiff_tasks:
                input_files.extend(self._get_files_by_task_and_link(mdiff_task.name, FileLink.OUTPUT))
                workflow.add_edge(mdiff_task.name, task_name)
            mconcat_task = self._generate_task('mConcatFit', task_name, input_files=input_files)
            workflow.add_node(task_name, task=mconcat_task)

            # mBgModel task
            task_name = self._generate_task_name('mBgModel')
            mbgmodel_task = self._generate_task('mBgModel', task_name,
                                                input_files=self._get_files_by_task_and_link(mconcat_task.name,
                                                                                             FileLink.OUTPUT))
            workflow.add_node(task_name, task=mbgmodel_task)
            workflow.add_edge(mconcat_task.name, task_name)

            # mBackground tasks
            mbg_tasks = []
            for task in mproject_tasks:
                input_files = self._get_files_by_task_and_link(task.name, FileLink.OUTPUT)
                input_files.extend(self._get_files_by_task_and_link(mbgmodel_task.name, FileLink.OUTPUT))
                task_name = self._generate_task_name('mBackground')
                mbg_task = self._generate_task('mBackground', task_name, input_files=input_files,
                                               files_recipe={FileLink.OUTPUT: {'.fits': 2}})
                workflow.add_node(task_name, task=mbg_task)
                workflow.add_edge(mbgmodel_task.name, task_name)
                mbg_tasks.append(mbg_task)

            # mImgtbl task
            task_name = self._generate_task_name('mImgtbl')
            input_files = []
            for task in mbg_tasks:
                input_files.append(self._get_files_by_task_and_link(task.name, FileLink.OUTPUT)[0])
                workflow.add_edge(task.name, task_name)
            mimg_task = self._generate_task('mImgtbl', task_name, input_files=input_files)
            workflow.add_node(task_name, task=mimg_task)

            # mAdd task
            task_name = self._generate_task_name('mAdd')
            input_files = []
            for task in mbg_tasks:
                input_files.extend(self._get_files_by_task_and_link(task.name, FileLink.OUTPUT))
            madd_task = self._generate_task('mAdd', task_name, input_files=input_files,
                                            files_recipe={FileLink.OUTPUT: {'.fits': 2}})
            workflow.add_node(task_name, task=madd_task)
            workflow.add_edge(mimg_task.name, task_name)
            madd_tasks.append(madd_task)

            # mViewer task
            task_name = self._generate_task_name('mViewer')
            mviewer_task = self._generate_task('mViewer', task_name,
                                               input_files=[self._get_files_by_task_and_link(madd_task.name,
                                                                                             FileLink.OUTPUT)[1]])
            workflow.add_node(task_name, task=mviewer_task)
            workflow.add_edge(madd_task.name, task_name)

        # colored mViewer task
        if len(madd_tasks) > 1:
            task_name = self._generate_task_name('mViewer')
            input_files = []
            for task in madd_tasks:
                input_files.append(self._get_files_by_task_and_link(task.name, FileLink.OUTPUT)[1])
                workflow.add_edge(task.name, task_name)
            mviewer_task = self._generate_task('mViewer', task_name, input_files=input_files)
            workflow.add_node(task_name, task=mviewer_task)

        self.workflows.append(workflow)
        return workflow

    def _workflow_recipe(self) -> Dict:
        """
        Recipe for generating synthetic traces of the Montage workflow. Recipes can be
        generated by using the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`.

        :return: A recipe in the form of a dictionary in which keys are task prefixes.
        :rtype: Dict[str, Any]
        """
        if self.dataset == MontageDataset.TWOMASS:
            return {
                "mProject": {
                    "runtime": {
                        "min": 0,
                        "max": 57.262,
                        "distribution": {
                            "name": "levy",
                            "params": [
                                -0.15150224852116434,
                                0.4809103177780148
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "pareto",
                                "params": [
                                    1.8196769868584715,
                                    -6.007916450615754,
                                    6.00791645038787
                                ]
                            },
                            "min": 0,
                            "max": 1594399
                        },
                        ".hdr": {
                            "distribution": {
                                "name": "skewnorm",
                                "params": [
                                    716480.4263935913,
                                    -0.001656034010816449,
                                    216.06420975997213
                                ]
                            },
                            "min": 277,
                            "max": 279
                        }
                    },
                    "output": {
                        ".fits": {
                            "distribution": {
                                "name": "wald",
                                "params": [
                                    -3.1461335409126936,
                                    10.477207652909453
                                ]
                            },
                            "min": 0,
                            "max": 4178880
                        }
                    }
                },
                "mDiffFit": {
                    "runtime": {
                        "min": 0.034,
                        "max": 46.525,
                        "distribution": {
                            "name": "beta",
                            "params": [
                                0.23296662250780056,
                                463.2955191436929,
                                -1.903355217182212e-25,
                                7569.50090953193
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.25284525166127914,
                                    -2.582591851437144e-24,
                                    2.575606368744306
                                ]
                            },
                            "min": 3329280,
                            "max": 4178880
                        },
                        ".hdr": {
                            "distribution": {
                                "name": "skewnorm",
                                "params": [
                                    1237686.3336703498,
                                    -0.001774478204270763,
                                    403.3634536136634
                                ]
                            },
                            "min": 277,
                            "max": 279
                        }
                    },
                    "output": {
                        ".txt": {
                            "distribution": {
                                "name": "cosine",
                                "params": [
                                    145.09419965704984,
                                    262.43239120561213
                                ]
                            },
                            "min": 51,
                            "max": 292
                        }
                    }
                },
                "mConcatFit": {
                    "runtime": {
                        "min": 2.593,
                        "max": 12.386,
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                147.5069177473742,
                                -15.577840787811667,
                                4067.4805108518813
                            ]
                        }
                    },
                    "input": {
                        ".txt": {
                            "distribution": {
                                "name": "cosine",
                                "params": [
                                    145.09419965704984,
                                    262.43239120561213
                                ]
                            },
                            "min": 51,
                            "max": 292
                        },
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 54340,
                            "max": 272870
                        }
                    },
                    "output": {
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 85697,
                            "max": 430977
                        }
                    }
                },
                "mBgModel": {
                    "runtime": {
                        "min": 20.192,
                        "max": 110.529,
                        "distribution": {
                            "name": "chi2",
                            "params": [
                                1.1039936299618875,
                                11.999999999999996,
                                0.7951785715324731
                            ]
                        }
                    },
                    "input": {
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 29868,
                            "max": 430977
                        }
                    },
                    "output": {
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 4456,
                            "max": 23706
                        }
                    }
                },
                "mBackground": {
                    "runtime": {
                        "min": 0.248,
                        "max": 35.961,
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.5358816399287538,
                                -1.4066448661907127e-27,
                                2.1683579540323255
                            ]
                        }
                    },
                    "input": {
                        ".tbl": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.641451719035309,
                                    -2.9474746116185896e-25,
                                    1.7270349932893398
                                ]
                            },
                            "min": 4456,
                            "max": 68710
                        },
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.5112505246020871,
                                    -7.107963785483312e-26,
                                    1.7971380184922583
                                ]
                            },
                            "min": 3329280,
                            "max": 4178880
                        }
                    },
                    "output": {
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.5112505246020871,
                                    -7.107963785483312e-26,
                                    1.7971380184922583
                                ]
                            },
                            "min": 3329280,
                            "max": 4178880
                        }
                    }
                },
                "mImgtbl": {
                    "runtime": {
                        "min": 0.443,
                        "max": 1.762,
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                147.5069177473742,
                                -15.577840787811667,
                                4067.4805108518813
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "pareto",
                                "params": [
                                    0.035118926741593066,
                                    -4.052883707888093e-10,
                                    1.2661614222770254e-10
                                ]
                            },
                            "min": 3329280,
                            "max": 4178880
                        },
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 13060,
                            "max": 68710
                        }
                    },
                    "output": {
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 35772,
                            "max": 188372
                        }
                    }
                },
                "mAdd": {
                    "runtime": {
                        "min": 0.945,
                        "max": 11.53,
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                147.5069177473742,
                                -15.577840787811667,
                                4067.4805108518813
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.5112505246020871,
                                    -7.107963785483312e-26,
                                    1.7971380184922583
                                ]
                            },
                            "min": 3329280,
                            "max": 4178880
                        },
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 35772,
                            "max": 188372
                        },
                        ".hdr": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 275,
                            "max": 277
                        }
                    },
                    "output": {
                        ".fits": {
                            "distribution": {
                                "name": "alpha",
                                "params": [
                                    147.5069177473742,
                                    -15.577840787811667,
                                    4067.4805108518813
                                ]
                            },
                            "min": 25922880,
                            "max": 414722880
                        }
                    }
                },
                "mViewer": {
                    "runtime": {
                        "min": 1.581,
                        "max": 57.913,
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                156.35362804649196,
                                -9.944864151497804,
                                4056.665627724218
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "alpha",
                                "params": [
                                    147.5069177473742,
                                    -15.577840787811667,
                                    4067.4805108518813
                                ]
                            },
                            "min": 25922880,
                            "max": 414722880
                        }
                    },
                    "output": {
                        ".jpg": {
                            "distribution": {
                                "name": "alpha",
                                "params": [
                                    156.35362804649196,
                                    -9.944864151497804,
                                    4056.665627724218
                                ]
                            },
                            "min": 1596605,
                            "max": 101394221
                        }
                    }
                }
            }

        if self.dataset == MontageDataset.DSS:
            return {
                "mProject": {
                    "runtime": {
                        "min": 160.278,
                        "max": 1263.481,
                        "distribution": {
                            "name": "levy",
                            "params": [
                                -0.6390031069644855,
                                2.982938876769154
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    -0.44487345435352355,
                                    1.4573579709548934
                                ]
                            },
                            "min": 5658400,
                            "max": 9655399
                        },
                        ".hdr": {
                            "distribution": {
                                "name": "chi2",
                                "params": [
                                    1.4241375536238945,
                                    -1.9183330928674237e-28,
                                    53.810520087279485
                                ]
                            },
                            "min": 277,
                            "max": 279
                        }
                    },
                    "output": {
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.7954321183916364,
                                    -2.855192723746931e-27,
                                    1.78687751430719
                                ]
                            },
                            "min": 15373440,
                            "max": 55002240
                        }
                    }
                },
                "mDiffFit": {
                    "runtime": {
                        "min": 0.004,
                        "max": 33.974,
                        "distribution": {
                            "name": "chi",
                            "params": [
                                0.04869951514418973,
                                -1.4034132380060528e-24,
                                593.2476321870599
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.5320355812304265,
                                    -1.4818265388258437e-24,
                                    3.2180437671262156
                                ]
                            },
                            "min": 15373440,
                            "max": 55002240
                        },
                        ".hdr": {
                            "distribution": {
                                "name": "skewnorm",
                                "params": [
                                    6310795.774193881,
                                    -0.0002866866728582616,
                                    316.70527916668425
                                ]
                            },
                            "min": 277,
                            "max": 279
                        }
                    },
                    "output": {
                        ".txt": {
                            "distribution": {
                                "name": "rdist",
                                "params": [
                                    2.127543340939409,
                                    33.03674854892499,
                                    36.97398780340566
                                ]
                            },
                            "min": 50,
                            "max": 302
                        }
                    }
                },
                "mConcatFit": {
                    "runtime": {
                        "min": 0.181,
                        "max": 1.754,
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                147.5069177473742,
                                -15.577840787811667,
                                4067.4805108518813
                            ]
                        }
                    },
                    "input": {
                        ".txt": {
                            "distribution": {
                                "name": "rdist",
                                "params": [
                                    2.127543340939409,
                                    33.03674854892499,
                                    36.97398780340566
                                ]
                            },
                            "min": 50,
                            "max": 302
                        },
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 1160,
                            "max": 242400
                        }
                    },
                    "output": {
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 1457,
                            "max": 32657
                        }
                    }
                },
                "mBgModel": {
                    "runtime": {
                        "min": 0.563,
                        "max": 6.696,
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                147.5069177473742,
                                -15.577840787811667,
                                4067.4805108518813
                            ]
                        }
                    },
                    "input": {
                        ".tbl": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.11751004387257374,
                                    5.999999999999999,
                                    0.010205255332337211
                                ]
                            },
                            "min": 1457,
                            "max": 32657
                        }
                    },
                    "output": {
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 276,
                            "max": 2751
                        }
                    }
                },
                "mBackground": {
                    "runtime": {
                        "min": 0.592,
                        "max": 114.103,
                        "distribution": {
                            "name": "fisk",
                            "params": [
                                0.7217049725251219,
                                -2.553644182656991e-27,
                                3.0795613361399266
                            ]
                        }
                    },
                    "input": {
                        ".tbl": {
                            "distribution": {
                                "name": "chi2",
                                "params": [
                                    1.4241375580566633,
                                    -3.280049988544237e-28,
                                    47.100934050350965
                                ]
                            },
                            "min": 276,
                            "max": 10192
                        },
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.7954321183916364,
                                    -2.855192723746931e-27,
                                    1.78687751430719
                                ]
                            },
                            "min": 15373440,
                            "max": 55002240
                        }
                    },
                    "output": {
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.7954321183916364,
                                    -2.855192723746931e-27,
                                    1.78687751430719
                                ]
                            },
                            "min": 15373440,
                            "max": 55002240
                        }
                    }
                },
                "mImgtbl": {
                    "runtime": {
                        "min": 0.152,
                        "max": 0.389,
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                147.5069177473742,
                                -15.577840787811667,
                                4067.4805108518813
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.14660276617638068,
                                    -1.20171432950012e-27,
                                    2.262211223349796
                                ]
                            },
                            "min": 15373440,
                            "max": 55002240
                        },
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 940,
                            "max": 10192
                        }
                    },
                    "output": {
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 2594,
                            "max": 28466
                        }
                    }
                },
                "mAdd": {
                    "runtime": {
                        "min": 0.523,
                        "max": 23.162,
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                147.5069177473742,
                                -15.577840787811667,
                                4067.4805108518813
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "fisk",
                                "params": [
                                    0.7954321183916364,
                                    -2.855192723746931e-27,
                                    1.78687751430719
                                ]
                            },
                            "min": 15373440,
                            "max": 55002240
                        },
                        ".tbl": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 2594,
                            "max": 28466
                        },
                        ".hdr": {
                            "distribution": {
                                "name": "levy",
                                "params": [
                                    11.99999999999986,
                                    1.5547041057555974e-13
                                ]
                            },
                            "min": 275,
                            "max": 277
                        }
                    },
                    "output": {
                        ".fits": {
                            "distribution": {
                                "name": "alpha",
                                "params": [
                                    147.5069177473742,
                                    -15.577840787811667,
                                    4067.4805108518813
                                ]
                            },
                            "min": 25922880,
                            "max": 414722880
                        }
                    }
                },
                "mViewer": {
                    "runtime": {
                        "min": 1.423,
                        "max": 68.368,
                        "distribution": {
                            "name": "alpha",
                            "params": [
                                156.35362804649196,
                                -9.944864151497804,
                                4056.665627724218
                            ]
                        }
                    },
                    "input": {
                        ".fits": {
                            "distribution": {
                                "name": "alpha",
                                "params": [
                                    147.5069177473742,
                                    -15.577840787811667,
                                    4067.4805108518813
                                ]
                            },
                            "min": 25922880,
                            "max": 414722880
                        }
                    },
                    "output": {
                        ".jpg": {
                            "distribution": {
                                "name": "alpha",
                                "params": [
                                    156.35362804649196,
                                    -9.944864151497804,
                                    4056.665627724218
                                ]
                            },
                            "min": 983389,
                            "max": 92116158
                        }
                    }
                }
            }
