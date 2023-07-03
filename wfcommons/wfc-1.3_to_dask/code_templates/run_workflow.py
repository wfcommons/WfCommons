"""
The contents of this file have been largely generated

You may want to modify this to tune to your taste.
It consists of two parts:
+ The first part is the workflow tasks definition. Feel free to change the values if you need
+ The second part is the dask code to execute the workflow. I wouldn't touch it too much if I were you (it's likely
better to modify the WFCommons JSON workflow definition in my opinion but you hold the chainsaw eventually).
"""
from helpers import execute_task
import random
from workflow_task import WorkflowTask


def run_workflow(client, simulate: bool, seed: int=42) -> list[WorkflowTask]:
# Generated code goes here
    return TASKS
