from wfcommons.wfgen.abstract_recipe import WorkflowRecipe
from wfcommons import WorkflowGenerator
import pathlib
from typing import Dict, Type, Union, List
import json
from numpy.random import choice
import numpy as np
import subprocess
from .data_gen import generate_sys_data, cleanup_sys_files


class WorkflowBenchmark():
    def __init__(self, Recipe: Type[WorkflowRecipe], num_tasks: int) -> None:        
        self.Recipe = Recipe
        self.num_tasks = num_tasks
       
    def create(self,
               cpu: float, 
               memory: float, 
               fileio: float,
               save_dir: Union[str, pathlib.Path],
               block_size: str,
               total_data_footprint: int,
               total_size: str,
               scope: str,
               max_prime: int,
               num_files: int,
               file_block_size: int,
            #    file_total_size: str = "2G",
               rw_ratio: float,
               test_mode: str) -> Dict:
    
        self._check_sysbench()
        save_dir = pathlib.Path(save_dir).resolve()
        self._check(cpu, fileio, memory)
        save_dir.mkdir(exist_ok=True, parents=True)

        generator = WorkflowGenerator(self.Recipe.from_num_tasks(self.num_tasks))
        workflow = generator.build_workflow()
        
        workflow.write_json(f'{save_dir.joinpath(workflow.name)}.json')
        with open(f'{save_dir.joinpath(workflow.name)}.json') as json_file:
            wf = json.load(json_file)
        
        for job in wf["workflow"]["jobs"]:
            job["benchmark"] = choice(["cpu", "fileio", "memory"], p=[cpu, fileio, memory])
            job["files"] =  []
            job.setdefault("command", {})
            job["command"]["program"] = f"wfperf-{job['benchmark']}.sh"
            job["command"]["arguments"] = []

        num_sys_files, total_files = self._in_out_files(wf) # num_files = system generated files and total_files = number of total files need in this workflow
        file_total_size = total_data_footprint/total_files # figuring out the size of every file 
        
        with open(f'{save_dir.joinpath(workflow.name)}.json', 'w') as fp:
            json.dump(wf, fp, indent=4)
        
        
        if cpu != 0:
            self._cpu_benchmark(save_dir, max_prime)
        
        if memory != 0:
            self._memory_benchmark(save_dir, block_size, total_size, scope)

        if fileio != 0:
            # first call to generate the input for task_need_input
            generate_sys_data(num_sys_files, file_total_size)
            # self._io_benchmark(save_dir, input_files=, num_files, file_block_size, file_total_size, rw_ratio, test_mode) #I don't know how to do this

        cleanup_sys_files()


    def _in_out_files(self, wf: Dict[str]) -> Dict[str]:
        task_need_input = [] # tasks that need input from 
        task_dont_need_input = 0
        for job in wf["workflow"]["jobs"]: 
            if job["benchmark"] == "fileio":
                parents = [parent for parent in job["parents"] if parent["benchmark"] == "fileio"] 
                if not parents:
                    task_need_input.append(job["name"])
                else:
                    task_dont_need_input+=1
                # job["is_root_fileio"] = len(parents) == 0
                for parent in parents:
                    job["files"] = [
                        {
                            "link": "input",
                            "name": item["name"],
                            "size": item["size"]
                        } 
                        for item in parent["files"] if item["link"] == "output"
                    ]
                        
        total_files = len(task_need_input)*2 + task_dont_need_input
        
        return task_need_input, total_files 

    def _check(self, cpu: float, fileio: float, memory: float):
        if not np.isclose(cpu + fileio + memory, 1.0):
            raise ValueError("cpu, fileio, and memory arguments must sum to 1.")

    def _check_sysbench(self,):
        proc = subprocess.Popen(["which", "sysbench"], stdout=subprocess.PIPE)
        out, _ = proc.communicate()
        if not out:
            raise FileNotFoundError("Sysbench not found. Please install sysbench: https://github.com/akopytov/sysbench")
            

    def _cpu_benchmark(self, save_dir: pathlib.Path, max_prime: int = 10000):
        save_dir.joinpath("wfperf-cpu.sh").write_text("\n".join([
            "#!/bin/bash", 
            "",
            f"sysbench cpu --cpu-max-prime={max_prime}"
        ]))
    
    def _memory_benchmark(self, 
                          save_dir: pathlib.Path,
                          block_size: str = "1K",
                          total_size: str = "100G",
                          scope: str = "global"):
        params = [
            f"--memory-block-size={block_size}",
            f"--memory-total-size={total_size}",
            f"--memory-scope={scope}"
        ]

        save_dir.joinpath("wfperf-memory.sh").write_text("\n".join([
            "#!/bin/bash", 
            "",
            f"sysbench memory {' '.join(params)}"
        ]))
    
    def _io_benchmark(self,
                      save_dir: pathlib.Path,
                      file_total_size: str,
                      input_files: List[str],
                      num_files: int = 1, # output to be generated 
                      file_block_size: int = 16384,
                      rw_ratio: float = 1.5,
                      test_mode: str = "seqwr"):
        
        
        _file_total_size =  len(input_files) * int(file_total_size)
        _num_files = f"--file-num={num_files}"

        params = [
            f"--file-total-size={file_total_size}",
            f"--file-test-mode=={test_mode}",
            f"--file-block-size={file_block_size}",
            f"--file-rw-ratio={rw_ratio}"
        ]

        params_run = [
            f"--file-total-size={_file_total_size}G",
            f"--file-test-mode={test_mode}",
            f"--file-block-size={file_block_size}",
            f"--file-rw-ratio={rw_ratio}"
        ]


        save_dir.joinpath("wfperf-fileio.sh").write_text("\n".join([
            "#!/bin/bash", 
            "",
            f" = "
            # f"sysbench fileio prepare {' '.join(params)}",
            f"sysbench fileio run {' '.join(params_run)}",
            f"sysbench fileio prepare {' '.join(_num_files, params)}",
        
        ])) 





