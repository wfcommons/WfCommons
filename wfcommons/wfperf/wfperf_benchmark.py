#!/usr/bin/python3

import pathlib
import argparse
import subprocess 

this_dir = pathlib.Path(__file__).resolve().parent

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", help="Benchmark")
    parser.add_argument("name", help="Task Name")
    return parser 

def main():
    parser = get_parser()
    args, other = parser.parse_known_args()

    if args.benchmark == "fileio":
        for path in this_dir.glob("*test_file*"):
            _, rest = path.name.split("test_file", 1)
            path.rename(path.parent.joinpath(f"test_file{rest}"))
        
        sysbench_args = [arg for arg in other if arg.startswith("--file")]
        proc = subprocess.Popen(["sysbench", "fileio", *sysbench_args, "run"])
        proc.wait()
        proc = subprocess.Popen(["sysbench", "fileio", *sysbench_args, "cleanup"])
        proc.wait()
        proc = subprocess.Popen(["sysbench", "fileio", *sysbench_args, "prepare"])
        proc.wait()
        for path in this_dir.glob("*test_file*"):
            path.rename(path.parent.joinpath(f"{args.name}_{path.name}"))

    elif args.benchmark == "cpu":
        sysbench_args = [arg for arg in other if arg.startswith("--cpu")]
        proc = subprocess.Popen(["sysbench", "cpu", *sysbench_args])
        proc.wait()
    elif args.benchmark == "memory":
        sysbench_args = [arg for arg in other if arg.startswith("--memory")]
        proc = subprocess.Popen(["sysbench", "memory", *sysbench_args])
        proc.wait()


if __name__ == "__main__":
    main()