from fractions import Fraction
from wfcommons.wfgen.abstract_recipe import WorkflowRecipe
from wfcommons import WorkflowGenerator
from typing import Dict, Union, List, Type
from numpy.random import choice
from .data_gen import generate_sys_data, cleanup_sys_files
import pathlib
import json
import numpy as np
import subprocess

class WorkflowBenchmark():
    def __init__(self, Recipe: Type[WorkflowRecipe], num_tasks: int) -> None:
        self.Recipe = Recipe
        self.num_tasks = num_tasks

    def create(self,
               save_dir: pathlib.Path,
               percent_cpu: float = 0.5,
               data_footprint: int = 100,
               test_mode: str = "seqwr",
               mem_total_size: str = "1000000G",
               block_size: str = "1K",
               scope: str = "global",   
               max_prime: int = 10000,
               file_block_size: int = 16384,
               rw_ratio: float = 1.5,
               max_time: int = 150,
               threads: int = 2, #cpu threads
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
        

        params = [f"save_dir={save_dir}",
                  f"--file-test-mode={test_mode}",
                  f"--file-total-size={data_footprint}G",
                  f"--file-block-size={file_block_size}",
                  f"--file-rw-ratio={rw_ratio}",
                  "--file-num=1",
                  f"--memory-block-size={block_size}",
                  f"--memory-scope={scope}",
                  f"--memory-total-size={mem_total_size}",
                  f"--cpu-max-prime={max_prime}",
                  f"--percent_cpu={percent_cpu}",
                  f"--time={max_time}"]

        for job in wf["workflow"]["jobs"]:
            # job["benchmark"] = choice(["cpu", "fileio", "memory"], p=[cpu, fileio, mem])
            job["files"] = []
            job.setdefault("command", {})
            job["command"]["program"] = f"wfperf_benchmark.py"
            job["command"]["arguments"] = [job["name"]]
            job["command"]["arguments"].extend(params)
            
        num_sys_files, num_total_files = self.input_files(wf)
        
        if verbose:
            print(f"Number of input files to be created by the system: {num_sys_files}")
            print(f"Total number of files used by the workflow: {num_total_files}")

        with open(f'{save_dir.joinpath(workflow.name)}.json', 'w') as fp:
            json.dump(wf, fp, indent=4)
            
        file_size = data_footprint / num_total_files
        if verbose:
            print(f"Every input/output file is of size: {file_size}")
        self.add_io_to_json(wf, file_size)
        
        if verbose:
            print("Generating system files.")
        generate_sys_data(num_sys_files, data_footprint)
        

        with open(f'{save_dir.joinpath(workflow.name)}.json', 'w') as fp:
            json.dump(wf, fp, indent=4)
        
        if verbose:
            print("Removing system files.")
        cleanup_sys_files()

    #fileio
    def input_files(self, wf: Dict[str, Dict]):
        tasks_need_input = 0
        tasks_dont_need_input = 0
        
        for job in wf["workflow"]["jobs"]:
            parents = [parent for parent in job["parents"]] 
            if not parents:
                tasks_need_input +=1                  
            else:
                tasks_dont_need_input += 1 
                            
        total_num_files = tasks_need_input*2 + tasks_dont_need_input
                
        return tasks_need_input, total_num_files         

    def add_io_to_json(self, wf: Dict[str, Dict], file_size: int) -> None:
        i=0
        all_jobs = {
            job["name"]: job
            for job in wf["workflow"]["jobs"]
        }
        
        for job in wf["workflow"]["jobs"]:
            job.setdefault("files", [])
            job["files"].append(
                {
                    "link": "output",
                    "name": f"{job['name']}_test_file.0",
                    "size": file_size
                }
            )

            parents = [parent for parent in job["parents"]] 
            if not parents:
                job["files"].append(
                    {
                        "link": "input",
                        "name": f"sys_test_file.{i}",
                        "size": file_size
                    } 
                )
                i+=1
            else:
                for parent in parents:
                    job["files"].extend(
                        [
                            {
                                "link": "input",
                                "name": item["name"],
                                "size": item["size"]
                            } 
                            for item in all_jobs[parent]["files"] if item["link"] == "output"
                        ]
                    )
        


    def _check_sysbench(self,):
        proc = subprocess.Popen(["which", "sysbench"], stdout=subprocess.PIPE)
        out, _ = proc.communicate()
        if not out:
            raise FileNotFoundError("Sysbench not found. Please install sysbench: https://github.com/aakopytov/sysbench")