from symbol import varargslist
from typing import Dict, IO, Optional

from importlib_metadata import entry_points
import argparse
from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path
from wfcommons.wfchef.recipes import CyclesRecipe, EpigenomicsRecipe, MontageRecipe, BlastRecipe, BwaRecipe, SeismologyRecipe, SoykbRecipe, GenomeRecipe, SrasearchRecipe
from ..wfchef.chef import get_recipes, get_recipe
from wfcommons import WorkflowGenerator
import pathlib
from typing import List, Dict
import json
import subprocess
from numpy.random import choice
import numpy as np

this_dir = pathlib.Path(__file__).resolve().parent


def create_perf_workflow(recipe: str, num_tasks: int, cpu: float, memory: float, fileio: float, params: List[str]) -> Dict:
    save_dir = pathlib.Path(this_dir.joinpath("benchmark"))
    save_dir.mkdir(exist_ok=True, parents=True)
 
    _recipe = get_recipe(recipe)
    wfname = recipe.split("R")[0]
    generator = WorkflowGenerator(_recipe.from_num_tasks(num_tasks))
    workflow = generator.build_workflow()
    workflow.write_json(f'{save_dir.joinpath(wfname)}.json')

    cat_params = {
        "fileio": [],
        "memory": [],
        "cpu": [],
    }
    cur_action = None 
    for param in params:
        if param.startswith("--memory"):
            cat_params["memory"].append(param)
        elif param.startswith("--file"):
            cat_params["fileio"].append(param)
        elif param.startswith("--cpu"):
            cat_params["cpu"].append(param)
        else:
            raise ValueError(f"{param} is not a valid option")

    with open(f'{save_dir.joinpath(wfname)}.json') as json_file:
        wf = json.load(json_file)
    
    for job in wf["workflow"]["jobs"]:
        job["benchmark"] = choice(["cpu", "fileio", "memory"], p=[cpu, fileio, memory])
        job["benchmark_params"] = cat_params[job["benchmark"]]
    
    with open(f'{save_dir.joinpath(wfname)}.json', 'w') as fp:
        json.dump(wf, fp, indent=4)
        
    save_dir.joinpath("cpu.sh").write_text("\n".join([
        "#!/bin/bash", 
        "",
        f"sysbench cpu {' '.join(cat_params['cpu'])}"
    ]))  
    save_dir.joinpath("fileio.sh").write_text("\n".join([
        "#!/bin/bash", 
        "",
        f"sysbench fileio prepare {' '.join(cat_params['fileio'])}",
        f"sysbench fileio run {' '.join(cat_params['fileio'])}",
        f"sysbench fileio cleanup {' '.join(cat_params['fileio'])}",
    ])) 
    save_dir.joinpath("memory.sh").write_text("\n".join([
        "#!/bin/bash", 
        "",
        f"sysbench memory {' '.join(cat_params['memory'])}"
    ]))
    # if benchmark in ["cpu", "memory"]:
    #     commands = [f"sysbench {benchmark} {' '.join(params)}"]
    # else:
    #     commands = [
    #         f"sysbench {benchmark} prepare {' '.join(params)}",
    #         f"sysbench {benchmark} run {' '.join(params)}",
    #         f"sysbench {benchmark} cleanup {' '.join(params)}",
    #     ]



    # How do I put executables in each node? And do I have to do this?

    # -----------------DOCUMENTATION---------------------------------
    # Refer to https://manpages.debian.org/testing/sysbench/sysbench.1.en.html for params format and options
    # Usage https://linuxtechlab.com/benchmark-linux-systems-install-sysbench-tool/ 

def bench_ls(benchmark: str):
    
    if benchmark == "cpu":
        print(subprocess.run('sysbench cpu help'))
    elif benchmark == "memory":
        print(subprocess.run('sysbench memory help'))
    else:
        print(subprocess.run('sysbench fileio help'))


class SysbenchArgumentParser(argparse.ArgumentParser):
    def print_help(self, file: Optional[IO[str]] = None) -> None:
        super().print_help(file=file)
        try:
            for action in ["cpu", "fileio", "memory"]:
                proc = subprocess.Popen(["sysbench", action, "help"], stdout=subprocess.PIPE)
                out, err = proc.communicate()
                if out:
                    print("\n".join(out.decode("utf-8").splitlines()[1:-1]))
                if err:
                    print(err)
        except ValueError:
            pass 
            
        

def get_parser() -> SysbenchArgumentParser:
    parser = SysbenchArgumentParser()
    # subparsers = parser.add_subparsers()
    
    parser.add_argument(
        "-r", "--recipe",
        choices = list(get_recipes()["name"].values),
        help="Recipe to base workflow on."
    )
    parser.add_argument(
        "-n", "--num-tasks",
        type=int, default=1000,
        help="Total number of tasks of workflow."
    )
    parser.add_argument(
        "--cpu", required=True,
        type=float,
        help="Percentage (0-1) of CPU-intesive nodes in workflow."
    )
    parser.add_argument(
        "--memory", required=True,
        type=float,
        help="Percentage (0-1) of memory-intesive nodes in workflow."
    )
    parser.add_argument(
        "--fileio", required=True,
        type=float,
        help="Percentage (0-1) of IO-intesive nodes in workflow."
    )

    # cpu_parser = subparsers.add_parser("cpu")
    # cpu_parser.set_defaults(benchmark="cpu")
    # fileio_parser = subparsers.add_parser("fileio")
    # fileio_parser.set_defaults(benchmark="fileio")
    # memory_parser = subparsers.add_parser("memory")
    # memory_parser.set_defaults(benchmark="memory")
    

    # for subparser in [cpu_parser, fileio_parser, memory_parser]:
    #     subparser.add_argument(
    #         "-r", "--recipe",
    #         choices = list(get_recipes()["name"].values),
    #         help="Recipe to base workflow on."
    #     )
    #     subparser.add_argument(
    #         "-n", "--num-tasks",
    #         type=int, default=1000,
    #         help="Total number of tasks of workflow."
    #     )

    return parser


def main():
    parser = get_parser()
    args, params = parser.parse_known_args()
    recipe = args.recipe
    num_tasks = args.num_tasks
    
    if not np.isclose(args.cpu + args.fileio + args.memory, 1.0):
        print("--cpu, --fileio, and --memory arguments must sum to 1.")
        return 
    



    create_perf_workflow(recipe, num_tasks, args.cpu, args.fileio, args.memory, params)

    
