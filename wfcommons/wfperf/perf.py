from typing import Dict
from ..wfgen.abstract_recipe import WorkflowRecipe
import argparse
from ..wfchef.chef import get_recipes

def create_perf_workflow(Recipe: WorkflowRecipe, num_nodes: int, type: str) -> Dict:
    pass



def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "-r", "--recipe",
        type=WorkflowRecipe,
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

