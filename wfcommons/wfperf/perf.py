from symbol import varargslist
from typing import Dict

from importlib_metadata import entry_points
import argparse

# from wfcommons.wfgen.abstract_recipe import WorkflowRecipe
from wfcommons.wfchef.recipes import CyclesRecipe, EpigenomicsRecipe, MontageRecipe, BlastRecipe, BwaRecipe, SeismologyRecipe, SoykbRecipe, GenomeRecipe, SrasearchRecipe
from ..wfchef.chef import get_recipes, get_recipe
from wfcommons import WorkflowGenerator

def create_perf_workflow(recipe: str, num_tasks: int, type: str) -> Dict:
    
    _recipe = get_recipe(recipe)
    wfname = recipe.split("R")[0]
    generator = WorkflowGenerator(_recipe.from_num_tasks(num_tasks))
    workflow = generator.build_workflow()
    workflow.write_json(f'{wfname}.json')



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
        "-t", "--type",
        type=str,
        choices =["CPU", "memory", "IO"],
        help="Total number of tasks of workflow."
    )


    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    recipe = args.recipe
    num_tasks = args.num_tasks
    _type = args.type

    create_perf_workflow(recipe, num_tasks, _type)

    
