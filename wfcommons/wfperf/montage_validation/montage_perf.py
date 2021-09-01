from wfcommons.wfchef.wfchef_abstract_recipe import WfChefWorkflowRecipe
from wfcommons.wfgen import WorkflowGenerator
from typing import Dict, Optional, Type, Tuple, Union
import pathlib
import json
import subprocess


this_dir = pathlib.Path(__file__).resolve().parent
class WorkflowBenchmark():
    def __init__(self, Recipe: Type[WfChefWorkflowRecipe], num_tasks: int) -> None:
        self.Recipe = Recipe
        self.num_tasks = num_tasks

    def create(self,
               save_dir: pathlib.Path,
               tasks: Dict[str, Tuple[int, int]],
               mem_total_size: str = "1000T",
               block_size: str = "4096",
               create: bool = " True",
               path: Optional[pathlib.Path] = " ",
               verbose: bool = False) -> Dict:


        if verbose:
            print("Checking if sysbench is installed.")
        self._check_sysbench()
        if verbose:
            print("Creating directory.")
        save_dir = pathlib.Path(save_dir).resolve()
        save_dir.mkdir(exist_ok=True, parents=True)

        if create:
            if verbose:
                print("Generating workflow")
            generator = WorkflowGenerator(self.Recipe.from_num_tasks(self.num_tasks))
            workflow = generator.build_workflow()
            workflow_savepath = save_dir.joinpath(f"{workflow.name}_{self.num_tasks}").with_suffix(".json")
            workflow.write_json(workflow_savepath)
            wf = json.loads(workflow_savepath.read_text())
        else:
            wf = json.loads(path.read_text())        
        
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
            job["command"]["program"] = "montage_sys_test.py"
            job_name = job["name"].rsplit("_", 1)[0]
            job["command"]["arguments"] = params[job_name]
            if "runtime" in job:
                del job["runtime"]

        json_path = save_dir.joinpath(f"{wf['name']}_{self.num_tasks}").with_suffix(".json")
        print("SAVING:", json_path)
        json_path.write_text(json.dumps(wf, indent=4))
            


    def _check_sysbench(self,):
        proc = subprocess.Popen(["which", "sysbench"], stdout=subprocess.PIPE)
        out, _ = proc.communicate()
        if not out:
            raise FileNotFoundError("Sysbench not found. Please install sysbench: https://github.com/aakopytov/sysbench")