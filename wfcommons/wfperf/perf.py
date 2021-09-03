import glob
import os
from wfcommons.wfchef.wfchef_abstract_recipe import WfChefWorkflowRecipe
from wfcommons import WorkflowGenerator
from typing import Dict, Optional, Union, List, Type
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
               percent_cpu: float = 0.6,
               percent_mem: float = 0.4,
               percent_io:float = 0.0,
               data_footprint: int = 100,
               test_mode: str = "seqwr",
               mem_total_size: str = "1000000G",
               block_size: str = "1K",
               scope: str = "global",   
               max_prime: int = 10000,
               file_block_size: int = 16384,
               rw_ratio: float = 1.5,
               max_time: int = 150,
               create: bool = " True",
               path: Optional[pathlib.Path] = " ",
               verbose: bool = False) -> Dict:

        save_dir = pathlib.Path(save_dir).resolve()
        save_dir.mkdir(exist_ok=True, parents=True)

        if create:
            if verbose:
                print("Generating workflow")
            generator = WorkflowGenerator(self.Recipe.from_num_tasks(self.num_tasks))
            workflow = generator.build_workflow()
            name = f"{workflow.name.split('-')[0]}-Benchmark"
            workflow_savepath = save_dir.joinpath(f"{name}-{self.num_tasks}").with_suffix(".json")
            workflow.write_json(workflow_savepath)
            wf = json.loads(workflow_savepath.read_text())
        else:
            wf = json.loads(path.read_text())    

        #Creating the lock files
        if verbose:
            print("Creating lock files.")
        with save_dir.joinpath("cores.txt").open("w+") as fp, save_dir.joinpath("cores.txt.lock").open("w+") as fpl:
            pass
    
        lock = save_dir.joinpath("cores.txt.lock")
        cores = save_dir.joinpath("cores.txt")

        #Setting the parameters for the arguments section of the JSON
        params = [f"save_dir={save_dir}",
                  f"--file-test-mode={test_mode}",
                  f"--file-total-size={data_footprint}G",
                  f"--file-block-size={file_block_size}",
                  f"--file-rw-ratio={rw_ratio}",
                  "--file-num=1",
                  "--forced-shutdown=0",  
                  f"--path-lock={lock}",
                  f"--path-cores={cores}",
                  f"--memory-block-size={block_size}",
                  f"--memory-scope={scope}",
                  f"--memory-total-size={mem_total_size}",
                  f"--cpu-max-prime={max_prime}",
                  f"--percent-cpu={percent_cpu}",
                  f"--percent-mem={percent_mem}",
                  f"--percent-io={percent_io}",
                  f"--time={max_time}"]

        wf["name"] = name
        for job in wf["workflow"]["jobs"]:
            job["files"] = []
            job.setdefault("command", {})
            job["command"]["program"] = f"{this_dir.joinpath('wfperf_benchmark.py')}"
            job["command"]["arguments"] = [job["name"]]
            job["command"]["arguments"].extend(params)
            if "runtime" in job:
                del job["runtime"]
            

        with open(f"{name}-{self.num_tasks}.json", 'w') as fp:
            json.dump(wf, fp, indent=4)

        num_sys_files, num_total_files = input_files(wf)
    
        if verbose:
            print(f"Number of input files to be created by the system: {num_sys_files}")
            print(f"Total number of files used by the workflow: {num_total_files}")  
        file_size = data_footprint / num_total_files
        if verbose:
            print(f"Every input/output file is of size: {file_size}")
        add_io_to_json(wf, file_size)
    
        if verbose:
            print("Generating system files.")
        generate_sys_data(num_sys_files, file_size)
        
        if verbose:
            print("Removing system files.")
        cleanup_sys_files()
    
        
        json_path = save_dir.joinpath(f"{name}-{self.num_tasks}").with_suffix(".json")
        print("SAVING:", json_path)
        json_path.write_text(json.dumps(wf, indent=4))
        
        return json_path

       
    def run(self, json_path: Union[str, pathlib.Path], savedir: Union[str, pathlib.Path], verbose: bool =True):
        if verbose:
            print("Running")
        try:
            wf = json.loads(json_path.read_text())
            with savedir.joinpath(f"run.txt").open("w+") as fp:
                procs: List[subprocess.Popen] = []
                for item in wf["workflow"]["jobs"]:
                    exec = item["command"]["program"]
                    arguments = item["command"]["arguments"]
                    prog = [
                        "time", "python", exec, 
                        item["name"].split("_")[0], 
                        *arguments, 
                    ]
                    procs.append(subprocess.Popen(prog), stdout=fp, stderr=fp)
                
                for proc in procs:
                    proc.wait()

        except:
            raise FileNotFoundError("Not able to find the executable.")


#Generates workflow's input data
def generate_sys_data(num_files: int, file_total_size: int, test_mode: str = "seqwr"):
    
    params = [
            f"--file-num={num_files}",
            f"--file-total-size={file_total_size}G",
            f"--file-test-mode={test_mode}"
        ]
    print(" ".join(["sysbench","fileio", *params,"prepare"]))
    proc = subprocess.Popen(["sysbench","fileio", *params,"prepare"], stdout=subprocess.PIPE)
    proc.wait()
    out, _ = proc.communicate()
    if not out:
        raise FileNotFoundError("Couldn't create files. Check parameters.")
    
    # have to change the name back to test_file on bash
    for t in glob.glob("test_file.*"):
        os.rename(t, f"sys_{t}")
        

#Removes files already used
def cleanup_sys_files():
    for t in glob.glob("*test_file.*"):
        new = t.split("_")[1:]
        new = "_".join(new)
        os.rename(t, new)
    proc = subprocess.Popen(["sysbench", "fileio", "cleanup"], stdout=subprocess.PIPE)
    proc.wait()
    out, _ = proc.communicate()
    if not out:
        raise FileNotFoundError("Couldn't delete files.")

#Adds input and output files to JSON
def add_io_to_json(wf: Dict[str, Dict], file_size: int) -> None:
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

#Calculate total number of files needed
def input_files(wf: Dict[str, Dict]):
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
    