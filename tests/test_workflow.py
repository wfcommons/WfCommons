#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pytest

from wfcommons.common import Workflow


class TestWorkflow:
    @pytest.mark.unit
    def test_workflow_creation(self):
        workflow = Workflow(
            name="Workflow Test",
            makespan=100.0
        )
        assert(workflow.name == "Workflow Test")
