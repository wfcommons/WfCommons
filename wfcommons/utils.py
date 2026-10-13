#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2020-2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
import warnings
warnings.filterwarnings('ignore')

import dateutil.parser
import json
import logging
import math
import pathlib
import re
import scipy.stats
import numpy as np
import operator as op

from datetime import datetime
from enum import Enum
from functools import reduce
from logging import Logger
from typing import Any, Dict, Optional, List, Tuple, Union

# Characters allowed in WfFormat task and file IDs (see "$defs/taskId" in the WfFormat schema)
WFFORMAT_ID_PATTERN = re.compile(r"^[0-9a-zA-Z_.#:/-]+$")
_WFFORMAT_ID_INVALID_CHARS = re.compile(r"[^0-9a-zA-Z_.#:/-]")


class NoValue(Enum):
    def __repr__(self):
        return '<%s.%s>' % (self.__class__.__name__, self.name)


def sanitize_wfformat_id(identifier: Any, replacement: str = "_") -> str:
    """
    Replace the characters that are not allowed in WfFormat task and file IDs.

    :param identifier: The raw identifier (e.g., a task name or a file name from a log).
    :type identifier: Any
    :param replacement: The string used in place of each invalid character.
    :type replacement: str

    :return: An identifier that matches the WfFormat ID pattern.
    :rtype: str
    """
    sanitized = _WFFORMAT_ID_INVALID_CHARS.sub(replacement, str(identifier))
    return sanitized if sanitized else replacement


def to_rfc3339(timestamp: Union[None, str, int, float, datetime]) -> Optional[str]:
    """
    Convert a timestamp into the RFC 3339 date-time format required by WfFormat. Naive
    timestamps are assumed to be in the local timezone.

    :param timestamp: A datetime, a POSIX timestamp, or a date-time string.
    :type timestamp: Union[None, str, int, float, datetime]

    :return: The timestamp in the RFC 3339 format, or None if no timestamp is provided.
    :rtype: Optional[str]
    """
    if timestamp is None or timestamp == "":
        return None
    if isinstance(timestamp, datetime):
        dt = timestamp
    elif isinstance(timestamp, (int, float)):
        dt = datetime.fromtimestamp(timestamp)
    else:
        dt = dateutil.parser.parse(str(timestamp))
    if dt.tzinfo is None:
        dt = dt.astimezone()
    return dt.isoformat()


def read_json(instance_filename: pathlib.Path) -> Dict[str, Any]:
    """
    Read the JSON from the file path.

    :param instance_filename: The absolute path of the instance file.
    :type instance_filename: str

    :return: The json object loaded with json data from the file
    :rtype: Dict[str, Any]
    """
    with open(instance_filename) as data:
        return json.load(data)


def best_fit_distribution(data: List[float], logger: Optional[Logger] = None) -> Tuple:
    """
    Fit a list of values to a distribution.

    :param data: List of values to be fitted to a distribution.
    :type data: List[float]
    :param logger: The logger uses to output debug information.
    :type logger: Logger

    :return: The name of the distribution and its parameters.
    :rtype: Tuple
    """
    if logger is None:
        logger = logging.getLogger("wfcommons")

    # get histogram of original data
    bins = math.ceil(len(data) / 10)
    

    normalized = (data - np.min(data)) / (np.max(data) - np.min(data)) if not np.min(data) == np.max(data) else np.min(data)
    y, x = np.histogram(normalized, bins=bins)
    if np.max(y) - np.min(y) > 0:
        y = (y - np.min(y)) / (np.max(y) - np.min(y))
    else:
        y = np.zeros(len(y))

    # best holders
    best_distribution = None
    best_params = (0.0, 1.0)
    best_sse = np.inf

    distribution_names: List[str] = [
        'alpha', 'arcsine', 'argus', 'beta', 'chi', 'chi2', 'dgamma', 'dweibull',
        'expon', 'fisk', 'gamma', 'gausshyper', 'levy', 'norm', 'pareto', 'rayleigh',
        'rdist', 'skewnorm', 'trapezoid', 'triang', 'uniform', 'wald'
    ]

    for dist_name in distribution_names:
        # Ignore warnings from data that can't be fit
        with warnings.catch_warnings():
            try:
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
            except Exception as e:
                logger.warning(f"WARNING: distribution \"{dist_name}\" failed ({e})")

    logger.debug(f'Best distribution fit: {best_distribution}')
    return best_distribution, best_params


def generate_rvs(distribution: Dict, min_value: float, max_value: float) -> float:
    """
    Generate a random variable from a distribution.

    :param distribution: Distribution dictionary (name and parameters).
    :type distribution: Dict
    :param min_value: Minimum value accepted as a random variable.
    :type min_value: float
    :param max_value: Maximum value accepted as a random variable.
    :type max_value: float

    :return: Random variable generated from a distribution.
    :rtype: float
    """
    if not distribution or distribution == "None":
        return min_value

    params = distribution['params']
    kwargs = params[:-2]
    rvs: float = max(0.1, getattr(scipy.stats, distribution['name']).rvs(*kwargs, loc=params[-2], scale=params[-1]))
    return rvs * max_value


def ncr(n: int, r: int) -> int:
    """
    Calculate the number of combinations.

    :param n: The number of items.
    :type n: int
    :param r: The number of items being chosen at a time.
    :type r: int

    :return: The number of combinations.
    :rtype: int
    """
    r = min(r, n - r)
    numerator = reduce(op.mul, range(n, n - r, -1), 1)
    denominator = reduce(op.mul, range(1, r + 1), 1)
    return numerator // denominator
