import pathlib
import argparse
import subprocess 
import os
from typing import List
from wfcommons.wfperf.montage_validation.lock import lock_core, unlock_core

import threading
import subprocess


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Task Name")
    # parser.add_argument("--save", type=pathlib.Path, required=True, help="directory to save to.")


    return parser

def main():
    parser = get_parser()
    args, other = parser.parse_known_args()

    print(f"Starting {args.name}")

    path_locked = pathlib.Path("/home/tgcoleman/tests/Montage/cores.txt.lock")
    path_cores = pathlib.Path("/home/tgcoleman/tests/Montage/cores.txt")
    
    sysbench_cpu_args = [arg for arg in other if arg.startswith("--cpu") or "time" in arg]
    percent_cpu = [arg for arg in other if arg.startswith("--percent")][0]
    cpu_threads = int(float(percent_cpu.split("=")[1])*10)
    mem_threads = int(10 - cpu_threads)

    print(f"cpu_threads={cpu_threads}, mem_threads={mem_threads}")
    print(sysbench_cpu_args)
    print("Starting CPU benchmark...")
    
    time = [arg for arg in sysbench_cpu_args if arg.startswith("--time")][0]
    time = int(time.split("=")[1])

    #if time is relevant do cpu and memory, else only cpu
    core = lock_core(path_locked, path_cores)
    print(f"{args.name} acquired core {core}")
    if time > 100:
        prog = [
            "sysbench", "cpu",
            *sysbench_cpu_args, f"--threads={cpu_threads}", "run"
        ]
        proc_cpu = subprocess.Popen(prog)
        
        os.sched_setaffinity(proc_cpu.pid, {core})

        print("Starting Memory benchmark...")
        sysbench_mem_args = [arg for arg in other if arg.startswith("--memory") or "time" in arg]
        prog = [
            "sysbench", "memory", "run",
            *sysbench_mem_args, f"--threads={mem_threads}"
        ]
        proc_mem = subprocess.Popen(prog)
        
        os.sched_setaffinity(proc_mem.pid, {core})
    else:
        proc_mem = None
        proc_cpu = subprocess.Popen(
            [
                "sysbench", "cpu",
                *sysbench_cpu_args, f"--threads={cpu_threads}", "run"
            ]
        )
        os.sched_setaffinity(proc_cpu.pid, {core})


    proc_cpu.wait()
    if proc_mem is not None:
        proc_mem.kill()
    unlock_core(path_locked, path_cores, core)
            
    
if __name__ == "__main__":
    main()
