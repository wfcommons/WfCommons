#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from typing import List
from .workflow.abstract_recipe import WorkflowRecipe
from ..common.workflow import Workflow
from ..errors import InvalidWorkflowTypeError


class WorkflowGenerator:
    def __init__(self, workflow_recipe: WorkflowRecipe = None) -> None:
        """
        :param workflow_recipe:
        """
        # sanity checks
        if not workflow_recipe or not isinstance(workflow_recipe, WorkflowRecipe):
            raise InvalidWorkflowTypeError("A WorkflowRecipe object should be provided.")

        self.workflow_recipe: WorkflowRecipe = workflow_recipe
        self.workflows: List[Workflow] = []

    def build_workflow(self, workflow_name: str = None) -> Workflow:
        """
        :param workflow_name: workflow name
        """
        workflow: Workflow = self.workflow_recipe.build_workflow(workflow_name)
        self.workflows.append(workflow)
        print("Generated a workflow with {} jobs".format(self.workflow_recipe.job_id_counter - 1))
        return workflow

    def build_workflows(self, num_workflows: int) -> List[Workflow]:
        """
        :param num_workflows: number of worklfows to be generated
        """
        workflows: List[Workflow] = []
        for _ in range(0, num_workflows):
            workflows.append(self.build_workflow())
        return workflows
