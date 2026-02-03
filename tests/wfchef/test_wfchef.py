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
from wfcommons.wfchef.utils import draw
from wfcommons.wfchef.chef import get_recipe
from wfcommons.wfchef.chef import create_recipe
from wfcommons.wfchef.chef import install_recipe
from wfcommons.wfchef.chef import uninstall_recipe
from wfcommons.wfchef.chef import ls_recipe
from wfcommons import WorkflowGenerator


class TestWfChef:

    @pytest.mark.unit
    # @pytest.mark.skip
    def test_recipe_management_functions(self) -> None:
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

        try:
            install_recipe("/bogus/bogus/whatever", verbose=True)
            raise RuntimeError("Should not be able to install a recipe given a bogus path")
        except FileNotFoundError as e:
            pass

        install_recipe(dirpath, verbose=True)
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

        # Verify the recipe can be used
        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Testing the recipe usage...\n")
        sys.stderr.write("=" * 60 + "\n")
        try:
            generator = WorkflowGenerator(SomenameRecipe.from_num_tasks(250))
            generator.build_workflow()
        except Exception as e:
            sys.stderr.write(f"✗ Failed to use installed recipe: {e}\n")
            raise

        # Test the graph drawing utility
        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Testing the graph drawing utility...\n")
        sys.stderr.write("=" * 60 + "\n")
        try:
            generator = WorkflowGenerator(SomenameRecipe.from_num_tasks(50))
            workflow = generator.build_workflow()
            plt_figure, plt_axis = draw(workflow) # Coverage
        except Exception as e:
            sys.stderr.write(f"✗ Failed to use installed recipe: {e}\n")
            raise

        try:
            recipe = get_recipe("somename_recipe")
            generator = WorkflowGenerator(recipe.from_num_tasks(250))
        except Exception as e:
            sys.stderr.write(f"✗ Failed to use get installed recipe by name: {e}\n")
            raise

        sys.stderr.write("✓ Successfully used recipe via direct import and by using get_recipe(name)\n")

        # Uninstall the recipe
        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Uninstalling the recipe...\n")
        sys.stderr.write("=" * 60 + "\n")

        uninstall_recipe("somename")
        sys.stderr.write("✓ Recipe uninstalled successfully\n")

        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Calling ls_recipe after the uninstall:\n")
        sys.stderr.write("=" * 60 + "\n")
        ls_recipe()

        # Verify the recipe can be used
        sys.stderr.write("\n" + "=" * 60 + "\n")
        sys.stderr.write("Testing the recipe usage after an uninstall...\n")
        sys.stderr.write("=" * 60 + "\n")
        try:
            generator = WorkflowGenerator(SomenameRecipe.from_num_tasks(250))
            generator.build_workflow()
            raise Exception("Should not be able to use a recipe after an uninstall")
        except Exception as e:
            pass

    @pytest.mark.unit
    def test_recipe_management_wfchef_main(self, monkeypatch, capsys) -> None:

        # Importing main for easier testing/coverage (compared to subprocessing the executable)
        from wfcommons.wfchef.chef import main

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

        # Calling main with --help
        with capsys.disabled():
            sys.stderr.write("\n" + "=" * 60 + "\n")
            sys.stderr.write("Calling wfchef script with '--help'...\n")
            sys.stderr.write("=" * 60 + "\n")
        monkeypatch.setattr(sys, 'argv', ["wfchef", "--help"])
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        captured = capsys.readouterr()
        if "usage" not in captured.out.lower():
            raise RuntimeError("wfchef script does not print usage information with -h flag")
        with capsys.disabled():
            sys.stderr.write("✓ Script called successfully\n")

        # Calling main with no argument, which should be the same as --help
        with capsys.disabled():
            sys.stderr.write("\n" + "=" * 60 + "\n")
            sys.stderr.write("Calling wfchef script with no argument...\n")
            sys.stderr.write("=" * 60 + "\n")
        monkeypatch.setattr(sys, 'argv', ["wfchef"])
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        captured = capsys.readouterr()
        if "usage" not in captured.out.lower():
            raise RuntimeError("wfchef script does not print usage information with -h flag")
        with capsys.disabled():
            sys.stderr.write("✓ Script called successfully\n")

        # Calling main with 'ls' command
        with capsys.disabled():
            sys.stderr.write("\n" + "=" * 60 + "\n")
            sys.stderr.write("Calling wfchef script with 'ls' command...\n")
            sys.stderr.write("=" * 60 + "\n")
        monkeypatch.setattr(sys, 'argv', ["wfchef", "ls"])
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        captured = capsys.readouterr()
        if "_recipe" not in captured.out.lower():
            raise RuntimeError("wfchef script does not list recipes with 'ls' command")
        num_recipes = sum(["_recipe" in x for x in captured.out.splitlines()])
        with capsys.disabled():
            sys.stderr.write(f"✓ Script called successfully ({num_recipes} recipes listed)\n")

        # Calling main with 'create' command with the --no-install flag
        with capsys.disabled():
            sys.stderr.write("\n" + "=" * 60 + "\n")
            sys.stderr.write("Calling wfchef script with 'create' command withe the --no-install flag...\n")
            sys.stderr.write("=" * 60 + "\n")
        monkeypatch.setattr(sys, 'argv', ["wfchef", "create", str(dirpath), "-o",  str(dirpath), "-n", "SomeRecipe", "--no-install"])
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        capsys.readouterr() # Clear output

        # Calling main with 'ls' command
        with capsys.disabled():
            sys.stderr.write("\n" + "=" * 60 + "\n")
            sys.stderr.write("Calling wfchef script with 'ls' command...\n")
            sys.stderr.write("=" * 60 + "\n")
        monkeypatch.setattr(sys, 'argv', ["wfchef", "ls"])
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        captured = capsys.readouterr()
        if "_recipe" not in captured.out.lower():
            raise RuntimeError("wfchef script does not list recipes with 'ls' command")
        new_num_recipes = sum(["_recipe" in x for x in captured.out.splitlines()])
        with capsys.disabled():
            sys.stderr.write(f"✓ Script called successfully ({new_num_recipes} recipes listed)\n")

        assert(num_recipes == new_num_recipes)

        # Calling main with 'create' command WITHOUT the --no-install flag
        with capsys.disabled():
            sys.stderr.write("\n" + "=" * 60 + "\n")
            sys.stderr.write("Calling wfchef script with 'create' command without the --no-install flag...\n")
            sys.stderr.write("=" * 60 + "\n")
        monkeypatch.setattr(sys, 'argv',
                            ["wfchef", "create", str(dirpath), "-o", str(dirpath), "-n", "SomeRecipe"])
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        capsys.readouterr()  # Clear output

        # Calling main with 'ls' command
        with capsys.disabled():
            sys.stderr.write("\n" + "=" * 60 + "\n")
            sys.stderr.write("Calling wfchef script with 'ls' command...\n")
            sys.stderr.write("=" * 60 + "\n")
        monkeypatch.setattr(sys, 'argv', ["wfchef", "ls"])
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        captured = capsys.readouterr()
        if "_recipe" not in captured.out.lower():
            raise RuntimeError("wfchef script does not list recipes with 'ls' command")
        new_num_recipes = sum(["_recipe" in x for x in captured.out.splitlines()])
        with capsys.disabled():
            sys.stderr.write(f"✓ Script called successfully ({new_num_recipes} recipes listed)\n")

        assert (num_recipes + 1 == new_num_recipes)

        # Calling main with 'uninstall' command
        with capsys.disabled():
            sys.stderr.write("\n" + "=" * 60 + "\n")
            sys.stderr.write("Calling wfchef script with 'uninstall' command...\n")
            sys.stderr.write("=" * 60 + "\n")
        monkeypatch.setattr(sys, 'argv', ["wfchef", "uninstall", "-n", "SomeRecipe"])
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        capsys.readouterr() # Clear output

        # Calling main with 'ls' command
        with capsys.disabled():
            sys.stderr.write("\n" + "=" * 60 + "\n")
            sys.stderr.write("Calling wfchef script with 'ls' command...\n")
            sys.stderr.write("=" * 60 + "\n")
        monkeypatch.setattr(sys, 'argv', ["wfchef", "ls"])
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        captured = capsys.readouterr()
        if "_recipe" not in captured.out.lower():
            raise RuntimeError("wfchef script does not list recipes with 'ls' command")
        new_new_num_recipes = sum(["_recipe" in x for x in captured.out.splitlines()])
        with capsys.disabled():
            sys.stderr.write(f"✓ Script called successfully ({new_new_num_recipes} recipes listed)\n")

        assert(new_new_num_recipes == num_recipes)
