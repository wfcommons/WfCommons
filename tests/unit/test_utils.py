#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import pathlib
import pytest
import wfcommons.utils

from typing import Dict


class TestUtils:

    @pytest.mark.unit
    def test_best_fit_distribution(self) -> None:
        pass

    @pytest.mark.parametrize(
        "distribution,min_value,max_value,variable",
        [
            pytest.param(None, 10, 100, 10),
            pytest.param({"name": "norm", "params": [0.08688656476267097, 0.2572832376513094]}, 10, 100, 10),
        ],
    ) 
    def test_generate_rvs(self, distribution: Dict, min_value: float, max_value: float, variable: float) -> None:
        assert(min_value <= wfcommons.utils.generate_rvs(distribution, min_value, max_value) <= max_value )

    @pytest.mark.parametrize(
        "n,r,combinations",
        [
            pytest.param(2, 2, 1),
            pytest.param(10, 2, 45),
            pytest.param(10, 3, 120),
        ],
    )   
    def test_ncr(self, n: int, r: int, combinations: int) -> None:
        assert(wfcommons.utils.ncr(n, r) == combinations)
