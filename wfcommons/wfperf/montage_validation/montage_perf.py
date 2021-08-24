from fractions import Fraction
from wfcommons.wfgen.abstract_recipe import WorkflowRecipe
from wfcommons import WorkflowGenerator
from typing import Dict, Union, List, Type, Tuple
from numpy.random import choice
import pathlib
import json
import subprocess

class WorkflowBenchmark():
    def __init__(self, Recipe: Type[WorkflowRecipe], num_tasks: int) -> None:
        self.Recipe = Recipe
        self.num_tasks = num_tasks

    def create(self,
               save_dir: pathlib.Path,
               tasks: Dict[str, Tuple[int, int]],
               mem_total_size: str = "1000T",
               block_size: str = "4096",
               verbose: bool = False) -> Dict:


        if verbose:
            print("Checking if the sysbench is installed.")
        self._check_sysbench()
        if verbose:
            print("Creating directory.")
        save_dir = pathlib.Path(save_dir).resolve()
        save_dir.mkdir(exist_ok=True, parents=True)

        if verbose:
            print("Generating workflow")
        generator = WorkflowGenerator(self.Recipe.from_num_tasks(self.num_tasks))
        workflow = generator.build_workflow()
        workflow.write_json(f'{save_dir.joinpath(workflow.name)}.json')

        with open(f'{save_dir.joinpath(workflow.name)}.json') as json_file:
            wf = json.load(json_file)

        params = {
            job: [
                f"--memory-block-size={block_size}",
                f"--memory-total-size={mem_total_size}",
                f"--cpu-max-prime={cpu_max_prime}",
                f"--percent_cpu={percent_cpu}",  
                f"--forced-shutdown=0",  
                f"--time={time}"
            ]
            for job, (cpu_max_prime, percent_cpu, time) in tasks.items()
        }         
        
        for job in wf["workflow"]["jobs"]:
            job["files"] = []
            job.setdefault("command", {})
            job["command"]["program"] = f"sys_test.py"
            job_name = job["name"].rsplit("_", 1)[0]
            job["command"]["arguments"] = params[job_name]


        with open(f'{save_dir.joinpath(workflow.name)}.json', 'w') as fp:
            json.dump(wf, fp, indent=4)


    def _check_sysbench(self,):
        proc = subprocess.Popen(["which", "sysbench"], stdout=subprocess.PIPE)
        out, _ = proc.communicate()
        if not out:
            raise FileNotFoundError("Sysbench not found. Please install sysbench: https://github.com/aakopytov/sysbench")