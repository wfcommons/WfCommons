import argparse
from os import mkdir
import subprocess
from typing import List
import pathlib


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="name of the task")
    parser.add_argument("-m", "--max-prime", help="max prime value")
    parser.add_argument("-p", "--percent-cpu", type=float, help="percentage of cpu usage")
    parser.add_argument("-s", "--save", help="Path to save results")

    return parser 


def main():
    parser = get_parser()
    args = parser.parse_args()
    name = args.name
    mp = args.max_prime
    savedir = pathlib.Path(args.save)
    cpu_threads = int(args.percent_cpu*10)
    mem_threads = 10 - cpu_threads
    savedir.mkdir(exist_ok=True, parents=True)
    

    sys_args =["sys_test.py",
           name,
           "--time=120",
           "--save",
           str(savedir),
           "--memory-hugetlb",
           f"--cpu-max-prime={mp}",
           "--memory-total-size=1000T",
           "--memory-block-size=4096",
           f"--percent-cpu={args.percent_cpu}"]
    
    with savedir.joinpath(f"run_10_{name}_{cpu_threads}_{mem_threads}.txt").open("w+") as fp: 
        for i in range(10):
            print(f"Starting {i} process.")
            proc = subprocess.Popen(["time", "python", *sys_args], stdout=fp, stderr=fp)
            proc.wait()
            print(f"Finished {i} process.")

if __name__ == "__main__":
    main()
