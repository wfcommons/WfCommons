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
import subprocess
import sys

from tests.test_helpers import _create_fresh_local_dir
from wfcommons.wfchef.chef import create_recipe
from wfcommons.wfchef.chef import install_recipe
from wfcommons.wfchef.chef import uninstall_recipe
from wfcommons.wfchef.chef import ls_recipe


class TestWfChef:

    @pytest.mark.unit
    def test_create_recipe(self) -> None:
        """
        Just calling the create_recipe function from chef.py directly (i.e., bypassing main())
        """

        dirpath = _create_fresh_local_dir("/tmp/recipe/")

        # Put a few JSON workflows in /tmp
        urls = [
            "https://raw.githubusercontent.com/wfcommons/WfInstances/refs/heads/main/makeflow/blast/blast-chameleon-small-001.json",
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

        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Creating recipe...\n")
        sys.stderr.write("=" * 60 + "\n")

        create_recipe(
            args["path"],
            args["out"],
            args["name"],
            cutoff=args["cutoff"],
            verbose=True
        )

        # Check that expected files are there
        sys.stderr.write("\nVerifying created files...\n")
        assert (dirpath / "pyproject.toml").exists(), "pyproject.toml not found"
        assert (dirpath / "wfchef_recipe_somename" / "__init__.py").exists(), "package __init__.py not found"
        assert (dirpath / "wfchef_recipe_somename" / "recipe.py").exists(), "recipe.py not found"
        assert (dirpath / "wfchef_recipe_somename" / "microstructures").exists(), "microstructures not found"
        sys.stderr.write("✓ All expected files created\n")

        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Calling ls_recipe before the install:\n")
        sys.stderr.write("=" * 60 + "\n")
        ls_recipe()

        # Install the recipe
        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Installing the recipe...\n")
        sys.stderr.write("=" * 60 + "\n")

        success = install_recipe(dirpath, verbose=True)
        assert success, "Recipe installation failed"
        sys.stderr.write("✓ Recipe installed successfully\n")

        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Calling ls_recipe after the install:\n")
        sys.stderr.write("=" * 60 + "\n")
        ls_recipe()

        # Verify the recipe can be loaded
        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Testing the recipe import...\n")
        sys.stderr.write("=" * 60 + "\n")

        try:
            from wfchef_recipe_somename import SomenameRecipe
            sys.stderr.write("✓ Successfully imported SomenameRecipe\n")
            sys.stderr.write(f"  Recipe class: {SomenameRecipe}\n")
            sys.stderr.write(f"  Module: {SomenameRecipe.__module__}\n")
        except ImportError as e:
            sys.stderr.write(f"✗ Failed to import recipe: {e}\n")
            raise

        # Uninstall the recipe
        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Uninstalling the recipe...\n")
        sys.stderr.write("=" * 60 + "\n")

        success = uninstall_recipe("somename")
        assert success, "Recipe uninstallation failed"
        sys.stderr.write("✓ Recipe uninstalled successfully\n")

        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Calling ls_recipe after the uninstall:\n")
        sys.stderr.write("=" * 60 + "\n")
        ls_recipe()

        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("TEST COMPLETED SUCCESSFULLY\n")
        sys.stderr.write("=" * 60 + "\n")


        # TODO: Do more extensive tests
        #  - Install/Uninstall the recipe
        #  - Use the recipe




