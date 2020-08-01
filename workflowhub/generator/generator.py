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
from .errors import InvalidWorkflowTypeError
from .workflow.abstract_recipe import WorkflowRecipe
from ..common.workflow import Workflow


class WorkflowGenerator:
    def __init__(self, workflow_recipe: WorkflowRecipe = None) -> None:
        # sanity checks
        if not workflow_recipe or not isinstance(workflow_recipe, WorkflowRecipe):
            raise InvalidWorkflowTypeError("A WorkflowRecipe object should be provided.")

        self.workflow_recipe: WorkflowRecipe = workflow_recipe
        self.workflows: List[Workflow] = []

    def build_workflow(self) -> Workflow:
        workflow: Workflow = self.workflow_recipe.build_workflow()
        self.workflows.append(workflow)
        return workflow

    def build_workflows(self, num_workflows: int) -> List[Workflow]:
        workflows: List[Workflow] = []
        for _ in range(0, num_workflows):
            workflows.append(self.build_workflow())
        return workflows
