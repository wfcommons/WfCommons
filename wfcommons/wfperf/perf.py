from wfcommons.wfgen.abstract_recipe import WorkflowRecipe
from wfcommons import WorkflowGenerator
import pathlib
from typing import Dict, Type, Union
import json
from numpy.random import choice
import numpy as np



class WorkflowBenchmark():
    def __init__(self, Recipe: Type[WorkflowRecipe], num_tasks: int) -> None:        
        self.Recipe = Recipe
        self.num_tasks = num_tasks
       
    def create(self,
               cpu: float, 
               memory: float, 
               fileio: float,
               save_dir: Union[str, pathlib.Path],
               block_size: str = "1K",
               total_size: str = "100G",
               scope: str = "global",
               max_prime: int = 10000,
               num_files: int = 128,
               file_block_size: int = 16384,
               file_total_size: str = "2G",
               rw_ratio: float = 1.5,
               test_mode: str = "seqwr") -> Dict:
    
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

        with open(f'{save_dir.joinpath(workflow.name)}.json', 'w') as fp:
            json.dump(wf, fp, indent=4)
        
        if cpu != 0:
            self._cpu_benchmark(save_dir, max_prime)
        
        if memory != 0:
            self._memory_benchmark(save_dir, block_size, total_size, scope)

        if fileio != 0:
            self._io_benchmark(save_dir, num_files, file_block_size, file_total_size, rw_ratio, test_mode)


    def _check(self, cpu: float, fileio: float, memory: float):
        if not np.isclose(cpu + fileio + memory, 1.0):
            raise ValueError("cpu, fileio, and memory arguments must sum to 1.")

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
                      num_files: int = 128,
                      file_block_size: int = 16384,
                      file_total_size: str = "2G",
                      rw_ratio: float = 1.5,
                      test_mode: str = "seqwr"):
        
        params = [
            f"--file-total-size={file_total_size}",
            f"--file-test-mode=={test_mode}",
            f"--file-block-size={file_block_size}",
            f"--file-num={num_files}",
            f"--file-rw-ratio={rw_ratio}"
        ]

        save_dir.joinpath("wfperf-fileio.sh").write_text("\n".join([
            "#!/bin/bash", 
            "",
            f"sysbench fileio prepare {' '.join(params)}",
            f"sysbench fileio run {' '.join(params)}",
            f"sysbench fileio cleanup {' '.join(params)}",
        ])) 





