from typing import List
from wfcommons.wfperf.montage_validation.montage_perf import WorkflowBenchmark
from wfcommons.wfchef.recipes import MontageRecipe
from tentativa_recipes.Tentativa import TentativaRecipe
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
    parser.add_argument("-t", "--num-tasks", help="Number os tasks when create is true.")

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
  
    for values in num_tasks.values():
        total += values

    return total

def main():
    parser = get_parser()
    args = parser.parse_args()
    savedir = pathlib.Path(args.save)
    if args.path:
        path = pathlib.Path(args.path)
    num_tasks = int(args.num_tasks)

    print("Running")
    
    # num_tasks = 65 
    # num_tasks = total_tasks()
    tasks = {'mProject': (12800000, 0.7, 120), 
             'mDiffFit': (24900000 , 0.7, 1), 
             'mConcatFit': (24900000 , 0.7, 5), 
             'mBgModel': (1910000, 0.7, 120), 
             'mBackground': (24900000 , 0.7, 1), 
             'mImgtbl': (24900000 , 0.7, 2),
             'mAdd': (1050000, 0.6, 120),
             'mViewer': (7400000, 0.6, 120)}

    bench = WorkflowBenchmark(TentativaRecipe, num_tasks)

    if args.create:
        
        if args.verbose:
            print("Creating Recipe...")
        json_path = bench.create(str(savedir), tasks, verbose=True)
        
    else:
        json_path = bench.create(str(savedir), tasks, create=False, path=path, verbose=True)

    bench.run(json_path, savedir)


    

    
    
if __name__ == "__main__":
    main()
