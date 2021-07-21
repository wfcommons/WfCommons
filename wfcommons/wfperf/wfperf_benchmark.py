#!/usr/bin/python3

import pathlib
import argparse
import subprocess 
import re


this_dir = pathlib.Path(__file__).resolve().parent

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    # parser.add_argument("benchmark", help="Benchmark")
    parser.add_argument("name", help="Task Name")
    return parser 

def main():
    parser = get_parser()
    args, other = parser.parse_known_args()
 
    # if args.benchmark == "fileio":
    for path in this_dir.glob("*test_file*"):
        _, rest = path.name.split("test_file", 1)
        path.rename(path.parent.joinpath(f"test_file{rest}"))
    
    sysbench_file_input_args = [arg for arg in other if arg.startswith("--file")]   
    to_remove = [item for item in sysbench_file_input_args if "num" in item]
    sysbench_file_input_args.remove(to_remove[0])
 
    proc = subprocess.Popen(["sysbench", "fileio", *sysbench_file_input_args, "run"])
    proc.wait()
    proc = subprocess.Popen(["sysbench", "fileio", *sysbench_file_input_args, "cleanup"])
    proc.wait()
   
    # elif args.benchmark == "cpu":
    sysbench_cpu_args = [arg for arg in other if arg.startswith("--cpu")]
    proc = subprocess.Popen(["sysbench", "cpu", "run", *sysbench_cpu_args])
    proc.wait()

    # elif args.benchmark == "memory":
    sysbench_mem_args = [arg for arg in other if arg.startswith("--memory")]
    proc = subprocess.Popen(["sysbench", "memory", "run", *sysbench_mem_args])
    proc.wait()


    sysbench_file_output_args = [arg for arg in other if arg.startswith("--file")]
    proc = subprocess.Popen(["sysbench", "fileio", *sysbench_file_output_args, "prepare"])
    proc.wait()
    
    for path in this_dir.glob("*test_file*"):
        path.rename(path.parent.joinpath(f"{args.name}_{path.name}"))

if __name__ == "__main__":
    main()