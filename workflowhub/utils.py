#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 The WorkflowHub Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import logging
import scipy.stats
import warnings
import numpy as np

from logging import Logger
from typing import Optional, List, Tuple
from .types import JsonDict


def read_json(trace_filename: str, logger: Optional[Logger] = None) -> JsonDict:
	"""
		Reads the json from the file path.
		returns the json object loaded with json data from the file
		:param trace_filename:
		:param logger: the logger uses to output debug information
		:return: json object
	"""
	if logger is None:
		logger = logging.getLogger("workflowhub")

	with open(trace_filename) as data:
		data_json: JsonDict = json.load(data)
		logger.info("parsed JSON trace: " + trace_filename)
		return data_json


def best_fit_distribution(data: List[float], bins: int = 100, logger: Optional[Logger] = None):
	if logger is None:
		logger = logging.getLogger("workflowhub")

	# get histogram of original data
	y, x = np.histogram(data, bins=bins, density=True)

	# best holders
	best_distribution = None
	best_params = (0.0, 1.0)
	best_sse = np.inf

	dist_names: List[str] = ['alpha', 'beta', 'cosine', 'gamma', 'lognorm', 'norm', 'pareto', 'rayleigh', 't',
							 'uniform', 'weibull_min', 'weibull_max']

	for dist_name in dist_names:
		# Ignore warnings from data that can't be fit
		with warnings.catch_warnings():
			warnings.filterwarnings('ignore')

			distribution = getattr(scipy.stats, dist_name)
			params = distribution.fit(y)

			# calculate fitted PDF and error with fit in distribution
			pdf = distribution.pdf(x, *params[:-2], loc=params[-2], scale=params[-1])
			sse = np.sum(np.power(y - pdf[0:bins], 2.0))

			# identify if this distribution is better
			if best_sse > sse > 0:
				best_distribution = distribution
				best_params = params
				best_sse = sse

	return best_distribution.name, best_params
