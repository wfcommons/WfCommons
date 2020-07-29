#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from .exceptions import InvalidWorkflowType
from .workflow.abstract_recipe import WorkflowRecipe


class WorkflowGenerator:
	def __init__(self, workflow_recipe=None):
		# sanity checks
		if not workflow_recipe or not isinstance(workflow_recipe, WorkflowRecipe):
			raise InvalidWorkflowType("A WorkflowRecipe object should be provided.")

		self.workflow_recipe = workflow_recipe
