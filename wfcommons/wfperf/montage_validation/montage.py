from typing import List
from wfcommons.wfperf.montage_validation.montage_perf import WorkflowBenchmark
from wfcommons.wfchef.recipes import MontageRecipe
import pathlib
import argparse
import json
import subprocess

this_dir = pathlib.Path(__file__).resolve().parent

def get_parser() ->  argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", help="Path to JSON")
    parser.add_argument("-c", "--create", action="store_true", help="Generate Workflow Benchmark when set.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Prints status information when set to true.")
    parser.add_argument("-s", "--save", help="Path to save directory.")

    return parser

def total_tasks():
    num_tasks = {'mProject': 64, 
                 'mDiffFit': 2016, 
                 'mConcatFit': 1, 
                 'mBgModel': 1, 
                 'mBackground': 64, 
                 'mImgtbl': 1, 
                 'mAdd': 1, 
                 'mViewer': 1}
    total = 0
    # for (k, v) in num_tasks.items():
        
        # if k == 'mViewer':
        #     num_tasks[k] = v*3 + 1 
        # else: 
        #     num_tasks[k] = v*3

    for values in num_tasks.values():
        total += values

    return total

def main():
    parser = get_parser()
    args = parser.parse_args()
    savedir = pathlib.Path(args.save)

    print("RUnning")

    
    path_locked = pathlib.Path("/home/tgcoleman/tests/Montage/cores.txt.lock")
    path_cores = pathlib.Path("/home/tgcoleman/tests/Montage/cores.txt")

    path_locked.write_text("")
    path_cores.write_text("")
  
    
    num_tasks = 65 
    #num_tasks = total_tasks()
    tasks = {'mProject': (12800000, 0.7, 120), 
             'mDiffFit': (24900000 , 0.7, 1), 
             'mConcatFit': (24900000 , 0.7, 5), 
             'mBgModel': (1910000, 0.7, 120), 
             'mBackground': (24900000 , 0.7, 1), 
             'mImgtbl': (24900000 , 0.7, 2),
             'mAdd': (1050000, 0.6, 120),
             'mViewer': (7400000, 0.6, 120)}

    if args.create:
        if args.verbose:
            print("Creating Recipe...")
        bench = WorkflowBenchmark(MontageRecipe, num_tasks)
        bench.create(str(savedir), tasks, verbose=True)

     
    json_path = savedir.joinpath(f"Montage-synthetic-instance_{num_tasks}.json")
    
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
                    # "--save", str(savedir)
                ]
                procs.append(subprocess.Popen(prog))
            
            for proc in procs:
                proc.wait()

    except:
        raise FileNotFoundError("Not able to find the executable.")


    

    
    
if __name__ == "__main__":
    main()
