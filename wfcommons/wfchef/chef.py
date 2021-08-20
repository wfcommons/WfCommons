#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import argparse
from email.mime import base
import json
import traceback
from typing import Dict, Optional, Union
import argparse
import networkx as nx
import numpy as np
import pandas as pd
import pathlib
import pickle
import pkg_resources
import subprocess
import traceback
from typing import Dict, Optional, Union
from stringcase import capitalcase
from .duplicate import duplicate, NoMicrostructuresError
from .find_microstructures import save_microstructures
from .utils import create_graph
from ..wfgen.abstract_recipe import WorkflowRecipe
from ..wfinstances.instance import Instance
from ..wfinstances.instance_analyzer import InstanceAnalyzer
import math

this_dir = pathlib.Path(__file__).resolve().parent
skeleton_path = this_dir.joinpath("skeletons")


def compare_rmse(synth_graph: nx.DiGraph, real_graph: nx.DiGraph):
    """
    Calculates the Root Mean Square Error of a synthetic instance created 
    based on the correspondent (in number of tasks) real-world sample.

    :param synth_graph: a synthetic instance created by WfCommons.
    :type synth_graph: networkX DiGraph 
    :param real_graph: the correspondent (in number of tasks) real-world workflow instance.
    :type real_graph: networkX DiGraph 

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


def find_err(workflow: Union[str, pathlib.Path],
             err_savepath: Optional[Union[str, pathlib.Path]] = None,
             always_update: bool = False,
             runs: int = 1) -> pd.DataFrame:
    """
    Creates a dataframe with the Root Mean Square Error of the synthetic instances created 
    based on the correspondent, w.r.t. number of tasks, real-world samples available at 
    `WfCommons WfInstances from Pegasus WMS GitHub <https://github.com/wfcommons/pegasus-instances> 
    and from Makeflow WMS GitHub <https://github.com/wfcommons/makeflow-instances>`_
    repositories.

    :param workflow: name (for samples available in WfCommons) or path to the real workflow instances.
    :type workflow: str or pathlib.Path 
    :param err_savepath: path to save the err (rmse) of all instances available into a csv.
    :type real_graph: str or pathlib.Path
    :param always_update: flag to set if the err needs to be updated or not (True: if new instances are added, False: otherwise).
    :type real_graph: bool 
    :param runs: number of times to repeat the err calculation process (due to randomization).
    :type runs: bool

    :return: dataframe with RMSE of all available instances.
    :rtype: pd.DataFrame
    """
    summary = json.loads(workflow.joinpath("summary.json").read_text())
    sorted_graphs = sorted([name for name, _ in summary["base_graphs"].items()],
                           key=lambda name: summary["base_graphs"][name]["order"])

    if err_savepath:
        err_savepath = pathlib.Path(err_savepath)
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

    :param workflow: name (for samples available in WfCommons) or path to the real workflow instances.
    :type workflow: str or pathlib.Path 
    :param err_savepath: path to save the err (rmse) of all instances available into a csv.
    :type real_graph: str or pathlib.Path
    :param always_update: flag to set if the err needs to be updated or not (True: if new instances are added, False: otherwise).
    :type real_graph: bool 
    :param runs: number of times to repeat the err calculation process (due to randomization).
    :type runs:bool

    :return: dataframe with RMSE of all available instances.
    :rtype: pd.DataFrame
    """
    analyzer = InstanceAnalyzer()
    task_types = set()

    for path in path_to_instances.glob("*.json"):
        instance = Instance(input_instance=str(path))
        analyzer.append_instance(instance)
        graph = create_graph(path)
        for node in graph.nodes:
            task_types.add(graph.nodes[node]["type"])

    stats_dict = analyzer.build_summary(task_types - {"DST", "SRC"}, include_raw_data=False)

    return stats_dict

def get_recipe(recipe: str) -> "Module":
    for entry_point in pkg_resources.iter_entry_points('workflow_recipes'):
        att  = entry_point.attrs[0]
        if att == recipe:
            return entry_point.load()

def get_recipes() -> pd.DataFrame:
    rows = []
    for entry_point in pkg_resources.iter_entry_points('workflow_recipes'):
        try:
            Recipe = entry_point.load()
            rows.append(
                [Recipe.__name__, entry_point.module_name, f"from {entry_point.module_name} import {Recipe.__name__}"])
        except Exception as e:
            traceback.print_exc()
            print(f"Could not load {entry_point.module_name}")
    return pd.DataFrame(rows, columns=["name", "module", "import command"])

def ls_recipe():
    """
    Inspired by UNIX `ls` command, it lists the recipes already installed into the system and 
    how to import it to use.
    """
    print(get_recipes())


def uninstall_recipe(module_name: str):
    """
    Uninstalls a recipe installed in the system.
    """

    for entry_point in pkg_resources.iter_entry_points('workflow_recipes'):
        if entry_point.module_name == module_name:
            print(f"Uninstalling package: wfchef.recipe.{module_name}")
            proc = subprocess.Popen(["pip", "uninstall", f"wfchef.recipe.{module_name}"])
            proc.wait()
            return

    print(f"Could not find recipe with module name {module_name} installed")


def create_recipe(path_to_instances: Union[str, pathlib.Path],
                  savedir: pathlib.Path,
                  wf_name: str,
                  cutoff: int = 4000,
                  verbose: bool = False,
                  runs: int = 1) -> "WorkflowRecipe":
    """
    Creates a recipe for a workflow application by automatically replacing custom information 
    from the recipe skeleton.

    :param path_to_instances: name (for samples available in WfCommons) or path to the real workflow instances.
    :type path_to_instances: str or pathlib.Path 
    :param savedir: path to save the recipe.
    :type savedir: pathlib.Path
    :param wf_name: name of the workflow application.
    :type wf_name: str 
    :param cutoff: when set, only consider instances of smaller or equal sizes.
    :type cutoff: int
    :param verbose: when set, prints status messages.
    :type cutoff: bool
    :param verbose: number of times to repeat the err calculation process (due to randomization).
    :type runs:bool
    """
    camelname = capitalcase(wf_name)
    savedir.mkdir(exist_ok=True, parents=True)
    dst = pathlib.Path(savedir, f"{savedir.stem}_recipes", wf_name).resolve()
    dst.mkdir(exist_ok=True, parents=True)

    if verbose:
        print(f"Finding microstructures")
    microstructures_path = dst.joinpath("microstructures")
    save_microstructures(path_to_instances, microstructures_path, img_type=None, cutoff=cutoff)

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
    with this_dir.joinpath(dst.joinpath("recipe.py")).open("w+") as fp:
        fp.write(skeleton_str)

    # recipe __init__.py
    dst.joinpath("__init__.py").write_text(f"from .recipe import {camelname}Recipe")

    # setup.py 
    with skeleton_path.joinpath("setup.py").open() as fp:
        skeleton_str = fp.read()

    skeleton_str = skeleton_str.replace("PACKAGE_NAME", savedir.stem)
    with this_dir.joinpath(dst.parent.parent.joinpath("setup.py")).open("w+") as fp:
        fp.write(skeleton_str)

    # __init__.py
    dst.parent.joinpath("__init__.py").touch(exist_ok=True)
    with dst.parent.joinpath("__init__.py").open("a") as fp:
        fp.write(f"from .{wf_name} import {camelname}Recipe\n")

    # MANIFEST
    with dst.parent.parent.joinpath("MANIFEST.in").open("a+") as fp:
        fp.write(f"graft {savedir.stem}_recipes/{wf_name}/microstructures/**\n")
        fp.write(f"graft {savedir.stem}_recipes/{wf_name}/microstructures\n")
        fp.write(f"graft {savedir.stem}_recipes/{wf_name}\n")

    # workflow_recipes
    with this_dir.joinpath(dst.parent.parent.joinpath("workflow_recipes.txt")).open("a+") as fp:
        fp.write(f"{wf_name}_recipe = {savedir.stem}_recipes.{wf_name}:{camelname}Recipe\n")

    if verbose:
        print(f"Analyzing Workflow Statistics")
    stats = analyzer_summary(path_to_instances)
    dst.joinpath("task_type_stats.json").write_text(json.dumps(stats))


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    ls_parser = subparsers.add_parser("ls")
    ls_parser.set_defaults(action=ls_recipe)

    uninstall_parser = subparsers.add_parser("uninstall")
    uninstall_parser.set_defaults(action=uninstall_recipe)
    uninstall_parser.add_argument("module_name", help="name of recipe module to uninstall")

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

    if args.action == ls_recipe:
        ls_recipe()
    elif args.action == uninstall_recipe:
        uninstall_recipe(args.module_name)
    elif args.action == create_recipe:
        create_recipe(args.path, args.out, args.name, cutoff=args.cutoff, verbose=True)

        if args.no_install:
            print("Done! To install the package, run: \n")
            print(f"  pip install {args.out}")
            print("\nor, in editable mode:\n")
            print(f"  pip install -e {args.out}")
        else:
            proc = subprocess.Popen(["pip", "install", str(args.out.resolve())])
            proc.wait()


if __name__ == "__main__":
    main()
