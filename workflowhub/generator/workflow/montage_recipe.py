#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from typing import Dict, Optional
from .abstract_recipe import WorkflowRecipe
from ...common.workflow import Workflow


class MontageRecipe(WorkflowRecipe):
	def __init__(self, num_bands: int, size: Optional[int], degree: Optional[float], num_jobs: Optional[int]):
		super().__init__("Montage")

		# sanity checks
		if num_bands < 1 or num_bands > 3:
			raise ValueError("The number of mosaic bands should be in the range [1,3].")

		self.num_bands: Optional[int] = num_bands
		self.size: Optional[int] = size
		self.degree: Optional[float] = degree
		self.num_jobs: Optional[int] = num_jobs

	@classmethod
	def from_data(cls, num_bands: int, size: int) -> 'MontageRecipe':
		if size < 1000:
			raise ValueError("The total size of the data should be larger than 1000.")
		return cls(num_bands=num_bands, size=size, degree=None, num_jobs=None)

	@classmethod
	def from_degree(cls, num_bands: int, degree: float) -> 'MontageRecipe':
		return cls(num_bands=num_bands, size=None, degree=degree, num_jobs=None)

	@classmethod
	def from_jobs(cls, num_bands: int, num_jobs: int) -> 'MontageRecipe':
		return cls(num_bands=num_bands, size=None, degree=None, num_jobs=num_jobs)

	def build_workflow(self) -> Workflow:
		if self.size:
			pass
		# print(self.recipe)
		return Workflow('Montage')

	def _workflow_recipe(self) -> Dict:
		pass
