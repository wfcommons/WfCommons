#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from abc import ABC, abstractmethod
from ...common.workflow import Workflow


class WorkflowRecipe(ABC):
	def __init__(self, name: str = None) -> None:
		self.name = name

	@abstractmethod
	def build_workflow(self) -> Workflow:
		pass
