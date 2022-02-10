import argparse
from os import mkdir
import subprocess
from typing import List
import pathlib

import pandas as pd


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", help= "csv with thread/max_prime descripton")
    parser.add_argument("-s", "--save", help="Path to save results")

    return parser 


def main():
    parser = get_parser()
    args = parser.parse_args()
    csv = pathlib.Path(args.csv)
    df = pd.read_csv(csv, index_col=0)
    df["cpu"] = df["cpu"].apply(lambda x: x/10)


    savedir = pathlib.Path(args.save)
    savedir.mkdir(exist_ok=True, parents=True)

    
    df = df[~df["real"]]
    
    for i, (task,_, cpu, max_prime,_) in df.iterrows():
        sys_args = [
            "/home/tgcoleman/tests/wfperf_calibration_plus_run_9_22/montage_sys_new.py",
            task,
            f"--max-prime={int(max_prime)}",
            f"--percent-cpu={cpu}"]

        # for i in range(10):
        out = savedir.joinpath("output", task, f"run_{i}_{int(cpu*10)}_{int(10 - cpu*10)}.txt")
        out.parent.mkdir(exist_ok=True, parents=True)
        with out.open("w+") as fp: 
            print(f"Starting {i} process.")
            proc = subprocess.Popen(["time", "python", *sys_args], stdout=fp, stderr=fp)
            proc.wait()
            print(f"Finished {i} process.")
            
             

if __name__ == "__main__":
    main()