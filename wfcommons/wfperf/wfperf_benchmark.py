#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import argparse
import glob
import pathlib
import os
import subprocess
import time

from filelock import FileLock

this_dir = pathlib.Path(__file__).resolve().parent


def lock_core(path_locked: pathlib.Path, path_cores: pathlib.Path) -> int:
    """Lock cores in use.
    """
    all_cores = set(range(os.cpu_count()))
    while True:
        with FileLock(path_locked) as lock:
            try:
                lock.acquire()
                taken_cores = {int(line) for line in path_cores.read_text().splitlines() if line.strip()}
                available = all_cores - taken_cores
                if available:
                    core = available.pop()
                    taken_cores.add(core)
                    path_cores.write_text("\n".join(map(str, taken_cores)))
                    return core

                print(f"All Cores are taken")
            finally:
                lock.release()
        time.sleep(1)


def unlock_core(path_locked: pathlib.Path, path_cores: pathlib.Path, core: int):
    """Unlock cores after execution is done.
    """
    with FileLock(path_locked) as lock:
        lock.acquire()
        try:
            taken_cores = {
                int(line) for line in path_cores.read_text().splitlines()
                if int(line) != core
            }
            path_cores.write_text("\n".join(map(str, taken_cores)))
        finally:
            lock.release()

def cpu_mem_benchmark(percent_cpu:float = 0.5, percent_mem:float = 0.5, cpu_work:int = 100, core:int = 7):
    """ Runs cpu and memory benchmark.
    """
    cpu_threads = int(percent_cpu*10)
    mem_threads = int(percent_mem*10)

    cpu_procs = []
    cpu_prog = [f"{this_dir.joinpath('cpu-benchmark')}", str(cpu_work)]
    mem_prog = ["stress", "--vm", "1"]

    for i in range(cpu_threads):
        cpu_proc = subprocess.Popen(cpu_prog)
        os.sched_setaffinity(cpu_proc.pid, {core})
        cpu_procs.append(cpu_proc)
    
    for _ in range(mem_threads):
        mem_proc = subprocess.Popen(mem_prog)
        os.sched_setaffinity(mem_proc.pid, {core})

    return cpu_procs
    
def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Task Name")
    parser.add_argument("--percent-cpu", type=float, help="percentage related to the number of cpu threads.")
    parser.add_argument("--path-lock", help="Path to lock file.")
    parser.add_argument("--path-cores", help="Path to cores file.")
    parser.add_argument("--cpu-work", default=100, help="Amount of CPU work.")
    parser.add_argument("--data", action='store_true', default=False, help="Whether to process IO")
    parser.add_argument("--file-size", help="Size of an input/output file.")
    return parser


def main():
    parser = get_parser()
    args, other = parser.parse_known_args()

    path_locked = pathlib.Path(args.path_lock)
    path_cores = pathlib.Path(args.path_cores)
    path_locked.write_text("")
    path_cores.write_text("")

    print(f"Starting {args.name}")

    if args.data:
        print("Starting IO benchmark...")
        for file in other:
            with open(file, "w+") as fp:
                fp.readlines()
                
    print("Starting CPU and Memory benchmark...")
    core = lock_core(path_locked, path_cores)
    print(f"{args.name} acquired core {core}")
    
    cpu_procs = cpu_mem_benchmark(percent_cpu=args.percent_cpu, 
                                  percent_mem=(1-args.percent_cpu),
                                  cpu_work=args.cpu_work,  
                                  core=core)
    for proc in cpu_procs:
        proc.wait()
    subprocess.Popen(["killall", "stress"])
    unlock_core(path_locked, path_cores, core)

    if args.data:
        with open(this_dir.joinpath(f"{args.name}_output.txt"), "wb") as fp:
            fp.write(os.urandom(int(args.file_size)))
        


if __name__ == "__main__":
    main()
