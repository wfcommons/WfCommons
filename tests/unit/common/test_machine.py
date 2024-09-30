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

from wfcommons.common import Machine, MachineSystem


cpu = {
    "coreCount": 48,
    "speedInMHz": 1200,
    "vendor": "Vendor Name"
}


class TestMachine:
   
    @pytest.mark.unit
    def test_machine_creation(self) -> None:
        machine = Machine(name="machine_1", cpu=cpu)       
        assert(machine.as_dict() == {
            "nodeName": "machine_1",
            "cpu": cpu
        })

    @pytest.mark.unit
    def test_detailed_machine_creation(self) -> None:
        machine = Machine(name="machine_1", 
                          cpu=cpu,
                          system=MachineSystem.LINUX,
                          architecture="x86_64",
                          memory=1000,
                          release="release_test",
                          hashcode="some_hashcode")
        assert(machine.as_dict() == {
            "nodeName": "machine_1",
            "system": MachineSystem.LINUX.value,
            "architecture": "x86_64",
            "memoryInBytes": 1000,
            "release": "release_test",
            "cpu": cpu
        })
