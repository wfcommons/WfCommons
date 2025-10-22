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

from wfcommons.common import File


class TestFile:
   
    @pytest.mark.unit
    def test_file_creation(self) -> None:
        file = File(file_id="file_123", size=100)
        assert(file.as_dict() == {
            "id": "file_123",
            "sizeInBytes": 100
        })

    @pytest.mark.unit
    def test_file_equality(self) -> None:
        file1 = File(file_id="file_1", size=100)
        file2 = File(file_id="file_1", size=200)
        file3 = File(file_id="file_2", size=100)
        file4 = File(file_id="file_1", size=100)
        assert (file1 != file2)
        assert (file1 != file3)
        assert (file1 == file4)

    @pytest.mark.unit
    def test_hash(self) -> None:
        file1 = File(file_id="file_1", size=100)
        file2 = File(file_id="file_1", size=200)
        file3 = File(file_id="file_2", size=100)
        file4 = File(file_id="file_1", size=100)
        assert (hash(file1) != hash(file2))
        assert (hash(file1) != hash(file2))
        assert (hash(file1) == hash(file4))

        dict = {file1: "f1", file2: "f2", file3: "f3", file4: "f4"}
        assert(dict[file1] == "f4")
        assert(dict[file2] == "f2")
        assert(dict[file3] == "f3")
        assert(dict[file4] == "f4")
