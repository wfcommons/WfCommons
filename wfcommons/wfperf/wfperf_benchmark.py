#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2022 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import argparse
import pathlib
import os
import subprocess
import time
import json 

from filelock import FileLock
from typing import List, Optional, Dict

this_dir = pathlib.Path(__file__).resolve().parent


def lock_core(path_locked: pathlib.Path,
              path_cores: pathlib.Path) -> int:
    """
    Lock cores in use.

    :param path_locked:
    :type path_locked: pathlib.Path
    :param path_cores:
    :type path_cores: pathlib.Path

    :return:
    :rtype: int
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


def unlock_core(path_locked: pathlib.Path,
                path_cores: pathlib.Path,
                core: int) -> None:
    """
    Unlock cores after execution is done.

    :param path_locked:
    :type path_locked: pathlib.Path
    :param path_cores:
    :type path_cores: pathlib.Path
    :param core:
    :type core: int
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


def cpu_mem_benchmark(cpu_threads: Optional[int] = 5,
                      mem_threads: Optional[int] = 5,
                      cpu_work: Optional[int] = 100,
                      core: Optional[int] = 7) -> List:
    """
    Run cpu and memory benchmark.

    :param cpu_threads:
    :type cpu_threads: Optional[int]
    :param mem_threads:
    :type mem_threads: Optional[int]
    :param cpu_work:
    :type cpu_work: Optional[int]
    :param core:
    :type core: Optional[int]

    :return:
    :rtype: List
    """
    total_mem_bytes = 100.0 / os.cpu_count()
    cpu_work_per_thread = int(cpu_work / cpu_threads)

    cpu_procs = []
    cpu_prog = [f"{this_dir.joinpath('cpu-benchmark')}", f"{cpu_work_per_thread}"]
    mem_prog = ["stress-ng", "--vm", f"{mem_threads}", "--vm-bytes", f"{total_mem_bytes}%", "--vm-keep"]

    for i in range(cpu_threads):
        cpu_proc = subprocess.Popen(cpu_prog)
        os.sched_setaffinity(cpu_proc.pid, {core})
        cpu_procs.append(cpu_proc)

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
    # parser.add_argument("--input-data", default=None, help="User input data size from JSON file.")
    # parser.add_argument("--data", action='store_true', default=False, help="Whether to process IO.")
    # parser.add_argument("--outputs-file-size", help="Size of output files that need to be created.")
    parser.add_argument("--out", help="output files name.")
    return parser

def io_read_benchmark_datafootprint(other):
    print("[WfPerf] Starting IO Read Benchmark...")
    for file in other:
        with open(file, "rb") as fp:
            print(f"[WfPerf]   Reading '{file}'")
            fp.readlines()
    
    print("[WfPerf] Completed IO Read Benchmark!\n")

def io_read_benchmark_user_input_data_size(other):
    print("[WfPerf] Starting IO Read Benchmark...")
    for file in other:
        with open(file, "rb") as fp:
            print(f"[WfPerf]   Reading '{file}'")
            fp.readlines()
    print("[WfPerf] Completed IO Read Benchmark!\n")
    
def io_write_benchmark_datafootprint(output_file, file_size):
    print(f"[WfPerf] Writing output file '{output_file}'\n")
    with open(output_file, "wb") as fp:
        fp.write(os.urandom(file_size)) 

def io_write_benchmark_user_input_data_size(output_file, outputs_file_size):
    print(f"[WfPerf] Writing output file '{output_file}'\n")

    for job_name, file_size in outputs_file_size.values:
        with open(f"{output_file}_{job_name}", "wb") as fp:
            fp.write(os.urandom(int(file_size))) 
    

def main():
    """Main program."""
    parser = get_parser()
    args, other = parser.parse_known_args()

    path_locked = pathlib.Path(args.path_lock)
    path_cores = pathlib.Path(args.path_cores)
    core = lock_core(path_locked, path_cores)

    print(f"[WfPerf] Starting {args.name} Benchmark\n")

    # if args.data:
    #     io_read_benchmark_datafootprint(other)
    if args.out:
        io_read_benchmark_user_input_data_size(other)

    print("[WfPerf] Starting CPU and Memory Benchmarks...")
    print(f"[WfPerf]  {args.name} acquired core {core}")

    cpu_procs = cpu_mem_benchmark(cpu_threads=int(10 * args.percent_cpu),
                                  mem_threads=int(10 - 10 * args.percent_cpu),
                                  cpu_work=int(args.cpu_work),
                                  core=core)
    for proc in cpu_procs:
        proc.wait()
    subprocess.Popen(["killall", "stress-ng"])
    print("[WfPerf] Completed CPU and Memory Benchmarks!\n")

    # if args.data:
    #     io_write_benchmark_datafootprint(args.out, args.file_size)
    if args.out:
        outputs_file_size = json.loads(args.outputs_file_size)
        io_write_benchmark_user_input_data_size(args.out, outputs_file_size)
    
    unlock_core(path_locked, path_cores, core)
    print("WfPerf Benchmark completed!")


if __name__ == "__main__":
    main()
