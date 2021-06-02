#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import logging

from logging import Logger
from typing import List, Optional

from .abstract_recipe import WorkflowRecipe
from ..common.workflow import Workflow


class WorkflowGenerator:
    """
    A generator of synthetic workflow instances based on workflow recipes obtained from the
    analysis of real workflow execution instances.

    :param workflow_recipe: The workflow recipe to be used for this generator.
    :type workflow_recipe: WorkflowRecipe
    :param logger: The logger where to log information/warning or errors (optional).
    :type logger: Logger
    """

    def __init__(self, workflow_recipe: WorkflowRecipe, logger: Optional[Logger] = None) -> None:
        """Create an object of the workflow generator."""
        # sanity checks
        if not workflow_recipe or not isinstance(workflow_recipe, WorkflowRecipe):
            raise TypeError("A WorkflowRecipe object should be provided.")

        self.logger: Logger = logging.getLogger(__name__) if logger is None else logger
        self.workflow_recipe: WorkflowRecipe = workflow_recipe
        self.workflows: List[Workflow] = []

    def build_workflow(self, workflow_name: Optional[str] = None) -> Workflow:
        """
        Generate a synthetic workflow instance based on the workflow recipe used to instantiate
        the generator.

        :param workflow_name: The workflow name.
        :type workflow_name: str

        :return: A synthetic workflow instance object.
        :rtype: Workflow
        """
        workflow: Workflow = self.workflow_recipe.build_workflow(workflow_name)
        self.workflows.append(workflow)
        self.logger.info(
            "Generated a synthetic workflow with {} tasks".format(self.workflow_recipe.task_id_counter - 1))
        return workflow

    def build_workflows(self, num_workflows: int) -> List[Workflow]:
        """
        Generate a number of synthetic workflow instances based on the workflow recipe used to
        instantiate the generator.

        :param num_workflows: The number of workflows to be generated.
        :type num_workflows: int

        :return: A list of synthetic workflow instance objects.
        :rtype: List[Workflow]
        """
        if num_workflows < 1:
            raise ValueError("The number of workflows should be at least 1.")

        workflows: List[Workflow] = []
        for _ in range(0, num_workflows):
            workflows.append(self.build_workflow())
        return workflows
