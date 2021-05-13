
import pathlib
import json
from re import sub
from ..trace.trace import Trace
from ..trace.trace_analyzer import TraceAnalyzer
from ..generator.workflow.abstract_recipe import WorkflowRecipe, Workflow
from typing import Dict, Optional, Union
import argparse
import pandas as pd
from stringcase import camelcase, snakecase
import pickle
from .duplicate import duplicate, NoMicrostructuresError
from .utils import create_graph, annotate
from .find_microstructures import find_microstructures, save_microstructures
import pandas as pd
import networkx as nx
import math
import subprocess
import numpy as np
import networkx as nx
import pkg_resources

this_dir = pathlib.Path(__file__).resolve().parent
skeleton_path = this_dir.joinpath("skeletons")

def compare_rmse(synth_graph: nx.DiGraph, real_graph: nx.DiGraph):
    synthetic = {}
    real = {}
    
    for node in synth_graph.nodes:
        _type = synth_graph.nodes[node]['type_hash']
        synthetic.setdefault(_type, 0)
        synthetic[_type] +=1 

    for node in real_graph.nodes:
        _type = real_graph.nodes[node]['type_hash']
        real.setdefault(_type, 0)
        real[_type] +=1 
    
    _types = ({*synthetic.keys(), *real.keys()})

    mse = math.sqrt(sum([
        (real.get(_type, 0) - synthetic.get(_type, 0))**2
        for _type in _types
    ]) / len(_types))
    return mse / real_graph.order()

def find_err(workflow: Union[str, pathlib.Path], 
             err_savepath: Optional[Union[str, pathlib.Path]] = None,
             always_update: bool = False,
             runs: int = 1) -> None:
    summary = json.loads(workflow.joinpath("summary.json").read_text())
    sorted_graphs = sorted([name for name, _ in summary["base_graphs"].items()], key=lambda name: summary["base_graphs"][name]["order"])
    
    if err_savepath:
        err_savepath = pathlib.Path(err_savepath)
        err_savepath.parent.mkdir(exist_ok=True, parents=True)
        
    labels = [graph for graph in sorted_graphs]
    rows = [[None for _ in range(len(sorted_graphs))] for _ in range(len(sorted_graphs))]
    df = None 
    for i, path in enumerate(sorted_graphs[1:], start=1):
        path = workflow.joinpath(path)
        wf_real = pickle.loads(path.joinpath("base_graph.pickle").read_bytes())

        for j, base in enumerate(sorted_graphs[:i+1]):             
            try:
                dists = []
                for _ in range(runs):
                    wf_synth = duplicate(
                        path=workflow,
                        base=base,
                        num_nodes=wf_real.order(),
                        interpolate_limit=summary["base_graphs"][base]["order"]
                    )
                    dists.append(compare_rmse(wf_synth, wf_real))
                rows[j][i] = np.median(dists)
            except NoMicrostructuresError:
                print(f"No Microstructures Error")
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
    analyzer = TraceAnalyzer()
    task_types = set()

    for path in path_to_instances.glob("*.json"):
        instance = Trace(input_trace=str(path))
        analyzer.append_trace(instance)
        graph = create_graph(path)
        for node in graph.nodes:
            task_types.add(graph.nodes[node]["type"])

    stats_dict = analyzer.build_summary(task_types - {"DST", "SRC"}, include_raw_data=False)
    return stats_dict 

def ls_recipe():
    rows = []
    for entry_point in pkg_resources.iter_entry_points('worfklow_recipes'):
        Recipe = entry_point.load()
        rows.append([Recipe.__name__, entry_point.module_name, f"from {entry_point.module_name} import {Recipe.__name__}"])
    df = pd.DataFrame(rows, columns=["name", "module", "import command"])
    print(df.to_string(index=None))

def uninstall_recipe(module_name: str):
    for entry_point in pkg_resources.iter_entry_points('worfklow_recipes'):
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
                      runs: int = 1)-> WorkflowRecipe:
    savedir.mkdir(exist_ok=True, parents=True)
    dst = pathlib.Path(savedir, snakecase(wf_name)).resolve()
    dst.mkdir(exist_ok=True, parents=True)
    
    if verbose:
        print(f"Finding Microstructures")
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
    skeleton_str = skeleton_str.replace("Skeleton", wf_name)
    skeleton_str = skeleton_str.replace("skeleton", snakecase(wf_name))
    with this_dir.joinpath(dst.joinpath("__init__.py")).open("w+") as fp:
        fp.write(skeleton_str)

    # setup.py 
    with skeleton_path.joinpath("setup.py").open() as fp:
        skeleton_str = fp.read() 
        
    skeleton_str = skeleton_str.replace("Skeleton", wf_name)
    skeleton_str = skeleton_str.replace("skeleton", snakecase(wf_name))
    with this_dir.joinpath(dst.parent.joinpath("setup.py")).open("w+") as fp:
        fp.write(skeleton_str)

    # MANIFEST
    with skeleton_path.joinpath("MANIFEST.in").open() as fp:
        skeleton_str = fp.read() 
        
    skeleton_str = skeleton_str.replace("Skeleton", wf_name)
    skeleton_str = skeleton_str.replace("skeleton", snakecase(wf_name))
    with this_dir.joinpath(dst.parent.joinpath("MANIFEST.in")).open("w+") as fp:
        fp.write(skeleton_str)

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

