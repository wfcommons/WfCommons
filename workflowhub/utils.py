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
import math
import scipy.stats
import warnings
import numpy as np

from enum import Enum
from logging import Logger
from typing import Dict, Optional, List
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


def best_fit_distribution(data: List[float], logger: Optional[Logger] = None):
    """
        :param data:
        :param logger: the logger uses to output debug information
        :return:
    """
    if logger is None:
        logger = logging.getLogger("workflowhub")

    # get histogram of original data
    bins = math.ceil(len(data) / 20)
    y, x = np.histogram(data, bins=bins)

    # best holders
    best_distribution = None
    best_params = (0.0, 1.0)
    best_sse = np.inf

    distribution_names: List[str] = ['alpha', 'arcsine', 'argus', 'beta', 'chi', 'chi2', 'cosine', 'dgamma', 'dweibull',
                                     'expon', 'fisk', 'gamma', 'gausshyper', 'levy', 'norm', 'pareto', 'rayleigh',
                                     'rdist', 'skewnorm', 'trapz', 'triang', 'uniform', 'wald']

    for dist_name in distribution_names:
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
                best_distribution = dist_name
                best_params = params
                best_sse = sse

    return best_distribution, best_params


def generate_rvs(distribution: Dict, min_value: float, max_value: float) -> float:
    """
    :param distribution:
    :param min_value:
    :param max_value:
    """
    if not distribution or distribution == "None":
        return min_value

    params = distribution['params']
    kwargs = params[:-2]
    rvs: float = getattr(scipy.stats, distribution['name']).rvs(*kwargs, loc=params[-2], scale=params[-1])
    rvs = max(min_value, rvs)
    rvs = min(max_value, rvs)
    return rvs


class NoValue(Enum):
    def __repr__(self):
        return '<%s.%s>' % (self.__class__.__name__, self.name)
