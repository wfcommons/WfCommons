#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import random


class UniformDistribution:
	def __init__(self, low, high):
		self.low = low
		self.high = high

	def __call__(self):
		return random.uniform(self.low, self.high)

	def __repr__(self):
		return "UniformDistribution(%s, %s)" % (self.low, self.high)


class NormalDistribution:
	def __init__(self, mu, sigma):
		self.mu = mu
		self.sigma = sigma

	def __call__(self):
		return random.normalvariate(self.mu, self.sigma)

	def __repr__(self):
		return "NormalDistribution(%s, %s)" % (self.mu, self.sigma)
