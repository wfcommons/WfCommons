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
import pathlib
import os
import subprocess
import time

from filelock import FileLock

this_dir = pathlib.Path(__file__).resolve().parent


def check_sysbench():
    """Check if sysbench is installed.
    """
    proc = subprocess.Popen(["which", "sysbench"], stdout=subprocess.PIPE)
    out, _ = proc.communicate()
    if not out:
        raise FileNotFoundError("Sysbench not found. Please install sysbench: https://github.com/aakopytov/sysbench")


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


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Task Name")
    parser.add_argument("--time", type=int, help="maximum compute time in seconds")
    parser.add_argument("--percent-cpu", type=float, help="percent of threads that will be cpu")
    parser.add_argument("--percent-mem", type=float, help="percent of threads that will be mem")
    parser.add_argument("--percent-io", type=float, help="percent of threads that will be io")
    parser.add_argument("--path-lock", help="Path to lock file.")
    parser.add_argument("--path-cores", help="Path to cores file.")
    return parser


def main():
    parser = get_parser()
    args, other = parser.parse_known_args()

    print("Checking if sysbench is installed.")
    check_sysbench()

    save_dir = [item for item in other if "save" in item][0]
    save_dir = pathlib.Path(save_dir.split("=")[1])

    path_locked = pathlib.Path(args.path_lock)
    path_cores = pathlib.Path(args.path_cores)
    path_locked.write_text("")
    path_cores.write_text("")

    print(f"Starting {args.name}")

    assert (0.0 <= args.percent_cpu + args.percent_mem + args.percent_io <= 1.0)

    if args.percent_io:
        for path in save_dir.glob("*test_file*"):
            _, rest = path.name.split("test_file", 1)
            path.rename(path.parent.joinpath(f"test_file{rest}"))

        sysbench_file_input_args = [
            arg for arg in other
            if arg.startswith("--file") and not arg.startswith("--file-num") or "forced" in arg
        ]

        core = lock_core(path_locked, path_cores)
        io_threads = int(args.percent_io * 10)

        print("Starting IO benchmark...")
        proc = subprocess.Popen(
            [
                "sysbench", "fileio", *sysbench_file_input_args, f"--threads={io_threads}", "run"
            ]
        )
        proc.wait()

        proc = subprocess.Popen(
            [
                "sysbench", "fileio", *sysbench_file_input_args, "cleanup"
            ]
        )
        proc.wait()

    else:
        core = lock_core(path_locked, path_cores)

    sysbench_cpu_args = [arg for arg in other if arg.startswith("--cpu") or "time" or arg and "forced" in arg]

    cpu_threads = int(args.percent_cpu * 10)
    mem_threads = int(args.percent_mem * 10)

    print(f"cpu_threads={cpu_threads}, mem_threads={mem_threads}")
    print("Starting CPU benchmark...")

    time = args.time

    print(f"{args.name} acquired core {core}")
    if time > 100:
        prog = [
            "sysbench", "cpu", *sysbench_cpu_args, f"--threads={cpu_threads}", "run"
        ]
        proc_cpu = subprocess.Popen(prog)

        os.sched_setaffinity(proc_cpu.pid, {core})

        print("Starting Memory benchmark...")
        sysbench_mem_args = [arg for arg in other if arg.startswith("--memory") or "time" in arg or "forced" in arg]
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
                "sysbench", "cpu", *sysbench_cpu_args, f"--threads={cpu_threads}", "run"
            ]
        )
        os.sched_setaffinity(proc_cpu.pid, {core})

    proc_cpu.wait()
    if proc_mem is not None:
        proc_mem.kill()

    unlock_core(path_locked, path_cores, core)

    if args.percent_io:
        print("Writing output...")
        sysbench_file_output_args = [arg for arg in other if arg.startswith("--file") or "forced" in arg]
        proc = subprocess.Popen(
            [
                "sysbench", "fileio", *sysbench_file_output_args, f"--threads={io_threads}", "prepare"
            ]
        )
        proc.wait()

        for path in this_dir.glob("*test_file*"):
            path.rename(path.parent.joinpath(f"{save_dir}/{args.name}_{path.name}"))


if __name__ == "__main__":
    main()
