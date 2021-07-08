from symbol import varargslist
from typing import Dict

from importlib_metadata import entry_points
import argparse
from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path
import wfcommons

# from wfcommons.wfgen.abstract_recipe import WorkflowRecipe
from wfcommons.wfchef.recipes import CyclesRecipe, EpigenomicsRecipe, MontageRecipe, BlastRecipe, BwaRecipe, SeismologyRecipe, SoykbRecipe, GenomeRecipe, SrasearchRecipe
from ..wfchef.chef import get_recipes, get_recipe
from wfcommons import WorkflowGenerator
import pathlib
from typing import List, Dict
import json

this_dir = pathlib.Path(__file__).resolve().parent


class NotAllParameters(Exception):
    pass

def create_perf_workflow(recipe: str, num_tasks: int, benchmark: str, params:List[str]) -> Dict:
    
    save_dir = pathlib.Path(this_dir.joinpath("benchmark"))
    save_dir.mkdir(exist_ok=True, parents=True)
 
    _recipe = get_recipe(recipe)
    wfname = recipe.split("R")[0]
    generator = WorkflowGenerator(_recipe.from_num_tasks(num_tasks))
    workflow = generator.build_workflow()
    workflow.write_json(f'{save_dir.joinpath(wfname)}.json')

    with open(f'{save_dir.joinpath(wfname)}.json') as json_file:
        wf = json.load(json_file)
    
    for job in wf["workflow"]["jobs"]:
        job["params"] = params
    
    with open(f'{save_dir.joinpath(wfname)}.json', 'w') as fp:
        json.dump(wf, fp, indent=4)
        

    # How do I put executables in each node? And do I have to do this?

    # maybe create perform  ance parameters class on common directory and use object on workflow
    # Parameters to be added to the JSON: 
    # IO: file-total size, mode, max_time, max_requests 
    # Mem: # of threads




    # -----------------DOCUMENTATION---------------------------------
    # Refer to https://manpages.debian.org/testing/sysbench/sysbench.1.en.html for params format and options

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    
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
        "-b", "--benchmark",
        type=str,
        choices =["CPU", "memory", "IO"],
        help="Total number of tasks of workflow."
    )


    return parser


def main():
    parser = get_parser()
    args,params = parser.parse_known_args()
    recipe = args.recipe
    num_tasks = args.num_tasks
    benchmark = args.benchmark

    create_perf_workflow(recipe, num_tasks, benchmark, params)

    
