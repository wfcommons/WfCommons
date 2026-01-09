#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import argparse
import math
import sys
import networkx as nx
import numpy as np
import pandas as pd
import pathlib
import pickle
from importlib.metadata import entry_points
import tomli_w
import subprocess
import traceback
import sys

from typing import Dict, Optional, Union
from stringcase import capitalcase

from .duplicate import duplicate, NoMicrostructuresError
from .find_microstructures import save_microstructures
from .utils import create_graph
from ..wfinstances.instance import Instance
from ..wfinstances.instance_analyzer import InstanceAnalyzer

this_dir = pathlib.Path(__file__).resolve().parent
skeleton_path = this_dir.joinpath("skeletons")


def compare_rmse(synth_graph: nx.DiGraph, real_graph: nx.DiGraph) -> float:
    """
    Calculate the Root Mean Square Error of a synthetic instance created
    based on the correspondent (in number of tasks) real-world sample.

    :param synth_graph: a synthetic instance created by WfCommons.
    :type synth_graph: networkX.DiGraph
    :param real_graph: the correspondent (in number of tasks) real-world workflow instance.
    :type real_graph: networkX.DiGraph

    :return: The RMSE between the synthetic instance and the real instance.
    :rtype: float
    """
    synthetic = {}
    real = {}

    for node in synth_graph.nodes:
        _type = synth_graph.nodes[node]['type_hash']
        synthetic.setdefault(_type, 0)
        synthetic[_type] += 1

    for node in real_graph.nodes:
        _type = real_graph.nodes[node]['type_hash']
        real.setdefault(_type, 0)
        real[_type] += 1

    _types = ({*synthetic.keys(), *real.keys()})

    mse = math.sqrt(sum([
        (real.get(_type, 0) - synthetic.get(_type, 0)) ** 2
        for _type in _types
    ]) / len(_types))
    return mse / real_graph.order()


def find_err(workflow: pathlib.Path,
             err_savepath: Optional[pathlib.Path] = None,
             always_update: Optional[bool] = False,
             runs: Optional[int] = 1) -> pd.DataFrame:
    """
    Creates a dataframe with the Root Mean Square Error of the synthetic instances created 
    based on the correspondent, w.r.t. number of tasks, real-world samples available at 
    `WfCommons WfInstances from Pegasus WMS GitHub <https://github.com/wfcommons/pegasus-instances> 
    and from Makeflow WMS GitHub <https://github.com/wfcommons/makeflow-instances>`_
    repositories.

    :param workflow: name (for samples available in WfCommons) or path to the real workflow instances.
    :type workflow: pathlib.Path
    :param err_savepath: path to save the err (rmse) of all instances available into a csv.
    :type err_savepath: Optional[pathlib.Path]
    :param always_update: flag to set if the err needs to be updated or not (True: if new instances are added, False: otherwise).
    :type always_update: Optional[bool]
    :param runs: number of times to repeat the err calculation process (due to randomization).
    :type runs: Optional[bool]

    :return: dataframe with RMSE of all available instances.
    :rtype: pd.DataFrame
    """
    summary = json.loads(workflow.joinpath("summary.json").read_text())
    sorted_graphs = sorted([name for name, _ in summary["base_graphs"].items()],
                           key=lambda name: summary["base_graphs"][name]["order"])

    if err_savepath:
        err_savepath.parent.mkdir(exist_ok=True, parents=True)

    labels = [graph for graph in sorted_graphs]
    rows = [[None for _ in range(len(sorted_graphs))] for _ in range(len(sorted_graphs))]
    df = None
    for i, path in enumerate(sorted_graphs[1:], start=1):
        path = workflow.joinpath(path)
        wf_real = pickle.loads(path.joinpath("base_graph.pickle").read_bytes())

        for j, base in enumerate(sorted_graphs[:i + 1]):
            try:
                dists = []
                for _ in range(runs):
                    wf_synth = duplicate(
                        path=workflow,
                        base=base,
                        num_nodes=wf_real.order()
                    )
                    dists.append(compare_rmse(wf_synth, wf_real))
                rows[j][i] = np.median(dists)
            except NoMicrostructuresError:
                print(f"No microstructures Error")
                continue

            if err_savepath is not None and always_update:
                df = pd.DataFrame(rows, columns=labels, index=labels)
                df = df.dropna(axis=1, how='all')
                df = df.dropna(axis=0, how='all')
                err_savepath.write_text(df.to_csv())

    df = pd.DataFrame(rows, columns=labels, index=labels)
    df = df.dropna(axis=1, how='all')
    df = df.dropna(axis=0, how='all')
    if err_savepath:
        err_savepath.write_text(df.to_csv())
    return df


def analyzer_summary(path_to_instances: pathlib.Path) -> Dict:
    """
    Creates a dataframe with the Root Mean Square Error of the synthetic instances created 
    based on the correspondent, w.r.t. number of tasks, real-world samples available at 
    `WfCommons WfInstances from Pegasus WMS GitHub <https://github.com/wfcommons/pegasus-instances> 
    and from Makeflow WMS GitHub <https://github.com/wfcommons/makeflow-instances>`_
    repositories.

    :param path_to_instances:
    :type path_to_instances: pathlib.Path

    :return:
    :rtype: Dict
    """
    analyzer = InstanceAnalyzer()
    task_types = set()

    for path in path_to_instances.glob("*.json"):
        instance = Instance(input_instance=path)
        analyzer.append_instance(instance)
        graph = create_graph(path)
        for node in graph.nodes:
            task_types.add(graph.nodes[node]["type"])

    stats_dict = analyzer.build_summary(task_types - {"DST", "SRC"}, include_raw_data=False)

    return stats_dict


def get_recipe(recipe: str) -> Optional[type]:
    """
    Load a recipe by name from installed entry points.

    :param recipe: Name of the recipe to load
    :return: Recipe class or None if not found
    """
    # For Python 3.10+, entry_points() returns a more convenient interface
    # For Python 3.9, you may need to use entry_points().get('workflow_recipes', [])
    try:
        eps = entry_points(group='workflow_recipes')
    except TypeError:
        # Python 3.9 compatibility
        eps = entry_points().get('workflow_recipes', [])

    for entry_point in eps:
        # In importlib.metadata, entry points have 'name' instead of 'attrs'
        if entry_point.name == recipe:
            return entry_point.load()

    return None


def get_recipes() -> pd.DataFrame:
    """
    Get a DataFrame of all available workflow recipes.

    :return: DataFrame with columns: name, module, import command
    """
    rows = []

    try:
        eps = entry_points(group='workflow_recipes')
    except TypeError:
        # Python 3.9 compatibility
        eps = entry_points().get('workflow_recipes', [])

    for entry_point in eps:
        try:
            Recipe = entry_point.load()
            # Extract module name from the entry point value
            module_name = entry_point.value.split(':')[0]
            class_name = Recipe.__name__
            rows.append([
                entry_point.name,  # Use entry point name instead of class name
                module_name,
                f"from {module_name} import {class_name}"
            ])
        except Exception as e:
            print(f"Could not load {entry_point.name}: {e}")
            traceback.print_exc()

    return pd.DataFrame(rows, columns=["name", "module", "import command"])


def ls_recipe():
    """
    Inspired by UNIX `ls` command, it lists the recipes already installed
    into the system and how to import it to use.
    """
    print(get_recipes())


def install_recipe(recipe_path: Union[str, pathlib.Path],
                   verbose: bool = False):
    """
    Installs a recipe from a local directory into the system.

    :param recipe_path: Path to the recipe directory (containing setup.py or pyproject.toml)
    :param verbose: If True, show detailed pip output
    """
    recipe_path = pathlib.Path(recipe_path).resolve()

    if not recipe_path.exists():
        raise FileNotFoundError(f"Recipe path does not exist: {recipe_path}")

    # Check for setup.py or pyproject.toml
    has_setup = recipe_path.joinpath("setup.py").exists()
    has_pyproject = recipe_path.joinpath("pyproject.toml").exists()

    if not (has_setup or has_pyproject):
        raise FileNotFoundError(f"No setup.py or pyproject.toml found in {recipe_path}")

    try:
        cmd = [sys.executable, "-m", "pip", "install"]

        # Add verbose flag before -e if needed
        if verbose:
            cmd.append("-v")

        cmd.append(str(recipe_path))

        print(f"Installing recipe from: {recipe_path}")
        print(f"Command: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"Installation failed: {result.stderr}")
        else:
            print(f"Successfully installed recipe from {recipe_path}")
            if verbose:
                print(result.stdout)

    except Exception as e:
        traceback.print_exc()
        raise RuntimeError(f"Could not install recipe from {recipe_path}: {e}")


def uninstall_recipe(recipe_name: str):
    """
    Uninstalls a recipe installed in the system.

    :param recipe_name: Name of the recipe to uninstall (e.g., 'somename' or 'somename_recipe')
    """
    # Remove '_recipe' suffix if present
    if recipe_name.endswith('_recipe'):
        recipe_name = recipe_name[:-7]

    package_name = f"wfchef-recipe-{recipe_name}"

    try:
        cmd = [sys.executable, "-m", "pip", "uninstall", "-y", package_name]
        # print(f"Uninstalling: {package_name}")
        # print(f"Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"Uninstall failed: {result.stderr}")

    except Exception as e:
        traceback.print_exc()
        raise RuntimeError(f"Could not uninstall recipe for {recipe_name}: {e}")


def create_recipe(path_to_instances: Union[str, pathlib.Path],
                  savedir: pathlib.Path,
                  wf_name: str,
                  cutoff: int = 4000,
                  verbose: bool = False,
                  runs: int = 1,
                  author: str = "Workflow Chef",
                  author_email: str = "workflow@example.com",
                  package_version: str = "0.1.0"):
    """
    Creates a standalone recipe package for a workflow application.

    :param path_to_instances: name (for samples available in WfCommons) or
                              path to the real workflow instances.
    :param savedir: path to save the recipe.
    :param wf_name: name of the workflow application.
    :param cutoff: when set, only consider instances of smaller or equal sizes.
    :param verbose: when set, prints status messages (and helpful how-to instructions!)
    :param runs: number of times to repeat the err calculation process
                 (due to randomization).
    :param author: package author name.
    :param author_email: package author email.
    :param package_version: initial package version.
    """
    try:
        import tomli_w
    except ImportError:
        raise ImportError(
            "tomli_w is required for pyproject.toml generation. "
            "Install it with: pip install tomli-w"
        )

    # Import these from your actual modules
    from stringcase import capitalcase

    # Note: You'll need to define these paths in your actual code
    # skeleton_path = pathlib.Path(__file__).parent.joinpath("skeleton")

    camelname = capitalcase(wf_name)

    # Create a standalone package name
    package_name = f"wfchef-recipe-{wf_name}"
    module_name = f"wfchef_recipe_{wf_name}"

    # Create simple directory structure: {module_name}/
    savedir.mkdir(exist_ok=True, parents=True)

    # Create the package directory
    recipe_dir = savedir.joinpath(module_name)
    recipe_dir.mkdir(exist_ok=True, parents=True)

    if verbose:
        print(f"Finding microstructures")

    microstructures_path = recipe_dir.joinpath("microstructures")
    save_microstructures(path_to_instances, microstructures_path,
                         img_type=None, cutoff=cutoff)

    if verbose:
        print(f"Generating Error Table")

    err_savepath = microstructures_path.joinpath("metric", "err.csv")
    err_savepath.parent.mkdir(exist_ok=True, parents=True)
    df = find_err(microstructures_path, runs=runs)
    err_savepath.write_text(df.to_csv())

    # Recipe
    with skeleton_path.joinpath("recipe.py").open() as fp:
        skeleton_str = fp.read()

    if verbose:
        print(f"Generating Recipe Code")

    skeleton_str = skeleton_str.replace("Skeleton", camelname)
    skeleton_str = skeleton_str.replace("skeleton", wf_name)
    with recipe_dir.joinpath("recipe.py").open("w+") as fp:
        fp.write(skeleton_str)

    # Package __init__.py - exports the recipe class
    recipe_dir.joinpath("__init__.py").write_text(
        f"\"\"\"WfChef recipe for {wf_name} workflow.\"\"\"\n\n"
        f"from .recipe import {camelname}Recipe\n\n"
        f"__version__ = '{package_version}'\n"
        f"__all__ = ['{camelname}Recipe']\n"
    )

    # Create pyproject.toml at the package root
    pyproject_path = savedir.joinpath("pyproject.toml")

    if verbose:
        print(f"Generating pyproject.toml at {pyproject_path}")

    # Create config for this standalone recipe package
    config = {
        "build-system": {
            "requires": ["setuptools>=61.0", "wheel"],
            "build-backend": "setuptools.build_meta"
        },
        "project": {
            "name": package_name,
            "version": package_version,
            "description": f"WfChef recipe for {wf_name} workflow",
            "authors": [{"name": author, "email": author_email}],
            "requires-python": ">=3.8",
            "dependencies": [
                "wfcommons>=1.0.0",
                "pandas",
                "numpy",
            ],
        },
        "tool": {
            "setuptools": {
                "packages": [module_name],
                "package-data": {
                    module_name: ["**/*"]
                }
            }
        }
    }

    # Add entry point for this recipe
    config["project"]["entry-points"] = {
        "workflow_recipes": {
            f"{wf_name}_recipe": f"{module_name}:{camelname}Recipe"
        }
    }

    # Add README if it exists
    readme_path = savedir.joinpath("README.md")
    if readme_path.exists():
        config["project"]["readme"] = "README.md"
    elif verbose:
        # Create a basic README
        readme_content = f"""# {package_name}

WfChef recipe for {wf_name} workflow.

## Installation

```bash
pip install -e .
```

## Usage

```python
from {module_name} import {camelname}Recipe

recipe = {camelname}Recipe()
# Use the recipe...
```

## Entry Point

This package registers the following workflow recipe entry point:
- `{wf_name}_recipe` -> `{module_name}:{camelname}Recipe`

You can load it using:
```python
from wfcommons.wfchef import get_recipe

Recipe = get_recipe("{wf_name}_recipe")
recipe = Recipe()
```
"""
        readme_path.write_text(readme_content)
        config["project"]["readme"] = "README.md"

    # Write pyproject.toml
    with pyproject_path.open("wb") as f:
        tomli_w.dump(config, f)

    if verbose:
        print(f"Created pyproject.toml with entry point: {wf_name}_recipe")

    if verbose:
        print(f"Analyzing Workflow Statistics")

    stats = analyzer_summary(path_to_instances)
    recipe_dir.joinpath("task_type_stats.json").write_text(json.dumps(stats))

    if verbose:
        print(f"\n{'=' * 60}")
        print(f"Recipe created successfully!")
        print(f"{'=' * 60}")
        print(f"Recipe location: {recipe_dir}")
        print(f"Package root: {savedir}")
        print(f"Package name: {package_name}")
        print(f"Module name: {module_name}")
        print(f"Entry point: {wf_name}_recipe -> {module_name}:{camelname}Recipe")
        print(f"\nTo install this recipe, run:")
        print(f"  pip install -e {savedir}")
        print(f"\nOr use the install_recipe function:")
        print(f"  install_recipe('{savedir}', editable=True)")
        print(f"\nAfter installation, import it as:")
        print(f"  from {module_name} import {camelname}Recipe")
        print(f"\nOr load via entry point:")
        print(f"  from wfcommons.wfchef import get_recipe")
        print(f"  Recipe = get_recipe('{wf_name}_recipe')")
        print(f"{'=' * 60}\n")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    ls_parser = subparsers.add_parser("ls")
    ls_parser.set_defaults(action=ls_recipe)

    uninstall_parser = subparsers.add_parser("uninstall")
    uninstall_parser.set_defaults(action=uninstall_recipe)
    # uninstall_parser.add_argument("module_name", help="name of recipe module to uninstall")
    uninstall_parser.add_argument("-n", "--name", help="name of the workflow to uninstall")
    uninstall_parser.add_argument("-o", "--out", help="directory where the recipe is located")

    create_parser = subparsers.add_parser("create")
    create_parser.set_defaults(action=create_recipe)
    create_parser.add_argument(
        "path", type=pathlib.Path, help="path to JSON workflow instances"
    )
    create_parser.add_argument(
        "--no-install",
        action="store_true",
        help="if set, automatically installs the package"
    )
    create_parser.add_argument(
        "-r", "--runs",
        default=1, type=int,
        help="number of runs to compute mean RMSE"
    )
    create_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="if set, print verbose logs"
    )
    create_parser.add_argument(
        "-o", "--out",
        required=True,
        type=pathlib.Path,
        help="place to save generated recipe to"
    )
    create_parser.add_argument(
        "-n", "--name",
        required=True,
        help="name of workflow to generate"
    )
    create_parser.add_argument(
        "-c", "--cutoff",
        type=int, default=4000,
        help="ignore workflow instances with greater than cutoff number of nodes. default is 4000."
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    if not hasattr(args, "action"):
        sys.argv.append("--help")
        parser.parse_args()
        
    if args.action == ls_recipe:
        ls_recipe()
    elif args.action == uninstall_recipe:
        uninstall_recipe(args.name, pathlib.Path(args.out))
    elif args.action == create_recipe:
        create_recipe(args.path, args.out, args.name, cutoff=args.cutoff, verbose=True)

        if args.no_install:
            print("Done! To install the package, run: \n")
            print(f"  pip install {args.out}")
            print("\nor, in editable mode:\n")
            print(f"  pip install -e {args.out}")
        else:
            proc = subprocess.Popen([sys.executable, "-m", "pip", "install", str(args.out.resolve())])
            proc.wait()


if __name__ == "__main__":
    main()
