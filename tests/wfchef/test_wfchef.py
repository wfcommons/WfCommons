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
import pathlib
import shutil
import requests

from wfcommons.wfchef.chef import create_recipe


class TestWfChef:

    #@pytest.mark.unit
    @pytest.mark.skip(reason="Temporarily disabled due to strange ModuleNotFoundError: No module named 'pkg_resources' error")
    def test_create_recipe(self) -> None:
        """
        Just calling the create_recipe function from chef.py directly (i.e., bypassing main())
        """

        dirpath = pathlib.Path("/tmp/recipe/")
        if dirpath.exists():
            shutil.rmtree(dirpath)
        dirpath.mkdir(parents=True, exist_ok=True)

        # Put a few JSON workflows in /tmp
        urls = ["https://raw.githubusercontent.com/wfcommons/WfInstances/refs/heads/main/makeflow/blast/blast-chameleon-small-001.json",
                "https://raw.githubusercontent.com/wfcommons/WfInstances/refs/heads/main/makeflow/blast/blast-chameleon-small-002.json",
                "https://raw.githubusercontent.com/wfcommons/WfInstances/refs/heads/main/makeflow/blast/blast-chameleon-small-003.json",
                ]
        for url in urls:
            response = requests.get(url)
            local_file_name = url.split("/")[-1]
            with open(dirpath / local_file_name, 'wb') as f:
                f.write(response.content)

        # Call create_recipe()
        args = {
            "path": dirpath,
            "out": dirpath,
            "name": "somename",
            "cutoff": 4000
        }
        create_recipe(args["path"], args["out"], args["name"], cutoff=args["cutoff"], verbose=True)

        # Check that some of the expected files are there
        assert((dirpath / "setup.py").exists())
        assert((dirpath / "recipe_recipes" / "__init__.py").exists())
        assert((dirpath / "recipe_recipes" / "__init__.py").exists())
        assert((dirpath / "recipe_recipes" / "somename" / "__init__.py").exists())
        assert((dirpath / "recipe_recipes" / "somename" / "__init__.py").exists())
        assert((dirpath / "recipe_recipes" / "somename" / "recipe.py").exists())
        assert((dirpath / "recipe_recipes" / "somename" / "microstructures").exists())



        # TODO: Do more extensive tests (including USING the recipe)




