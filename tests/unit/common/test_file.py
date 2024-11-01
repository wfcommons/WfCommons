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

from wfcommons.common import File, FileLink


class TestFile:
   
    @pytest.mark.unit
    def test_file_creation(self) -> None:
        file = File(file_id="file_123", size=100, link=FileLink.INPUT)
        assert(file.as_dict() == {
            "id": "file_123",
            "sizeInBytes": 100
        })
