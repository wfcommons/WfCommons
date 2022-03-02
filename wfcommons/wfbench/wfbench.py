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
from io import StringIO
import pandas as pd

from filelock import FileLock
from typing import List, Optional

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
    path_locked.touch(exist_ok=True)
    path_cores.touch(exist_ok=True)

    while True:
        with FileLock(path_locked) as lock:
            try:
                lock.acquire()
                taken_cores = {
                    int(line) for line in path_cores.read_text().splitlines() if line.strip()}
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
                      core: Optional[int] = None) -> List:
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
    cpu_prog = [
        f"{this_dir.joinpath('cpu-benchmark')}", f"{cpu_work_per_thread}"]
    mem_prog = ["stress-ng", "--vm", f"{mem_threads}",
                "--vm-bytes", f"{total_mem_bytes}%", "--vm-keep"]

    for i in range(cpu_threads):
        cpu_proc = subprocess.Popen(cpu_prog)
        if core:
            os.sched_setaffinity(cpu_proc.pid, {core})
        cpu_procs.append(cpu_proc)

    mem_proc = subprocess.Popen(mem_prog)
    if core:
        os.sched_setaffinity(mem_proc.pid, {core})

    return cpu_procs

def get_available_gpus():
    proc = subprocess.Popen(["nvidia-smi", "--query-gpu=utilization.gpu", "--format=csv"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, _ = proc.communicate()
    df = pd.read_csv(StringIO(stdout.decode("utf-8")), sep=" ")
    return df[df["utilization.gpu"] <= 5].index.to_list()

def gpu_benchmark(work, device):
    gpu_prog = [f"CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES={device} {this_dir.joinpath('pi-cuda')} {work}"]
    print(gpu_prog)
    subprocess.Popen(gpu_prog, shell=True)  

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Task Name")
    parser.add_argument("--percent-cpu", default=0.5, type=float,
                        help="percentage related to the number of cpu threads.")
    parser.add_argument("--path-lock", default=None, help="Path to lock file.")
    parser.add_argument("--path-cores", default=None,
                        help="Path to cores file.")
    parser.add_argument("--cpu-work", default=100, help="Amount of CPU work.")
    parser.add_argument("--gpu-work", default=100, help="Amount of GPU work.")
    parser.add_argument("--out", help="output files name.")
    return parser


def io_read_benchmark_user_input_data_size(other):
    print("[WfBench] Starting IO Read Benchmark...")
    for file in other:
        with open(this_dir.joinpath(file), "rb") as fp:
            print(f"[WfBench]   Reading '{file}'")
            fp.readlines()
    print("[WfBench] Completed IO Read Benchmark!\n")


def io_write_benchmark_user_input_data_size(outputs):
    for task_name, file_size in outputs.items():
        print(f"[WfBench] Writing output file '{task_name}'\n")
        with open(this_dir.joinpath(task_name), "wb") as fp:
            fp.write(os.urandom(int(file_size)))


def main():
    """Main program."""
    parser = get_parser()
    args, other = parser.parse_known_args()

    core = None
    if args.path_lock and args.path_cores:
        path_locked = pathlib.Path(args.path_lock)
        path_cores = pathlib.Path(args.path_cores)
        core = lock_core(path_locked, path_cores)

    print(f"[WfBench] Starting {args.name} Benchmark\n")

    if args.out:
        io_read_benchmark_user_input_data_size(other)
    
    print("[WfBench] Starting GPU Benchmark...")
    available_gpus = get_available_gpus() #checking for available GPUs

    if not available_gpus:
        print("No GPU available")
    else:
        device = available_gpus[0]
        print(f"Running on GPU {device}")
        gpu_benchmark(this_dir, args.gpu_work, device)
    
    print("[WfBench] Starting CPU and Memory Benchmarks...")
    if core:
        print(f"[WfBench]  {args.name} acquired core {core}")

    cpu_procs = cpu_mem_benchmark(cpu_threads=int(10 * args.percent_cpu),
                                  mem_threads=int(10 - 10 * args.percent_cpu),
                                  cpu_work=int(args.cpu_work),
                                  core=core)
    for proc in cpu_procs:
        proc.wait()
    mem_kill = subprocess.Popen(["killall", "stress-ng"])
    mem_kill.wait()
    print("[WfBench] Completed CPU and Memory Benchmarks!\n")

    if args.out:
        outputs = json.loads(args.out.replace("'", '"'))
        io_write_benchmark_user_input_data_size(outputs)

    if core:
        unlock_core(path_locked, path_cores, core)

    print("WfBench Benchmark completed!")


if __name__ == "__main__":
    main()
