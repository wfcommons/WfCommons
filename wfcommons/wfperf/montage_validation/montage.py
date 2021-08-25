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
    parser.add_argument("path", help="Path to JSON")
    parser.add_argument("-c", "--create", help="Generate Workflow Benchmark when set.")
    parser.add_argument("-v", "--verbose", default=False, help="Prints status information when set to true.")
    parser.add_argument("-s", "--save", help="Path to save directory.")

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
    for (k, v) in num_tasks.items():
        if k == 'mViewer':
            num_tasks[k] = v*3 + 1 
        else: 
            num_tasks[k] = v*3

    for values in num_tasks.values():
        total += values

    return total

def main():
    parser = get_parser()
    args = parser.parse_args()
    json_path = pathlib.Path(args.path) 
    
    
    num_tasks = total_tasks()
    tasks = {'mProject': (12800000, 7, 120), 
             'mDiffFit': (24900000 , 7, 1), 
             'mConcatFit': (24900000 , 7, 5), 
             'mBgModel': (1910000, 7, 120), 
             'mBackground': (24900000 , 7, 1), 
             'mImgtbl': (24900000 , 7, 2),
             'mAdd': (1050000, 6, 120),
             'mViewer': (7400000, 6, 120)}

    if args.create:
        if args.verbose:
            print("Creating Recipe...")
        bench = WorkflowBenchmark(MontageRecipe, num_tasks)
        bench.create("/home/tgcoleman/tests/Montage", tasks, verbose=True)
    
    try:
        with open(json_path) as json_file:
            if args.verbose:
                print("Loading Recipe...")
            wf = json.load(json_file)
            
            with pathlib.Path(args.save).joinpath(f"run.txt").open("w+") as fp:
                procs: List[subprocess.Popen] = []
                for item in wf["workflow"]["jobs"]:
                    exec = item["command"]["program"]
                    arguments = item["command"]["arguments"]
                    if args.verbose:
                        print(f"Executing task:{item['name']}.")
                    procs.append(subprocess.Popen(["time",exec, *arguments], stdout=fp, stderr=fp))
                
                for proc in procs:
                    proc.wait()

    except:
        print("Not able to find the executable.")


    

    
    
if __name__ == "__main__":
    main()