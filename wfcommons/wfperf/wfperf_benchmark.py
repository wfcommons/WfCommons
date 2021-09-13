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
    parser.add_argument("--time", type=int, help="Maximum compute time in seconds")
    parser.add_argument("--percent-cpu", type=float, help="percentage related to the number of cpu threads.")
    parser.add_argument("--path-lock", help="Path to lock file.")
    parser.add_argument("--path-cores", help="Path to cores file.")
    parser.add_argument("--out", help="Output filename.")
    parser.add_argument("--data", action='store_true', help="Whether to process IO")
    parser.set_defaults(data=False)
    return parser


def main():
    parser = get_parser()
    args, other = parser.parse_known_args()

    print("Checking if sysbench is installed.")
    check_sysbench()

    path_locked = pathlib.Path(args.path_lock)
    path_cores = pathlib.Path(args.path_cores)
    path_locked.write_text("")
    path_cores.write_text("")

    print(f"Starting {args.name}")

    if args.data:

        counter = 0
        for file in this_dir.glob("*test_file*"):
            file.rename(file.parent.joinpath(f"test_file.{counter}"))
            counter += 1

        sysbench_file_input_args = [
            arg for arg in other
            if arg.startswith("--file-test-mode") or arg.startswith("--file-block-size=") or arg.startswith(
                "--file-rw-ratio=")
        ]

        core = lock_core(path_locked, path_cores)

        print("Starting IO benchmark...")
        proc = subprocess.Popen(
            [
                "sysbench", "fileio", *sysbench_file_input_args, f"--file-num={counter}", "run"
            ]
        )
        proc.wait()

        proc = subprocess.Popen(
            [
                "sysbench", "fileio", *sysbench_file_input_args, f"--file-num={counter}", "cleanup"
            ]
        )
        proc.wait()

    else:
        core = lock_core(path_locked, path_cores)

    sysbench_cpu_args = [arg for arg in other if arg.startswith("--cpu") or "forced" in arg]

    percent_mem = 1 - args.percent_cpu
    cpu_threads = int(args.percent_cpu * 10)
    mem_threads = int(percent_mem * 10)

    print(f"cpu_threads={cpu_threads}, mem_threads={mem_threads}")
    print("Starting CPU benchmark...")

    print(f"{args.name} acquired core {core}")
    if args.time > 100:
        prog = [
            "sysbench", "cpu", *sysbench_cpu_args, f"--threads={cpu_threads}", f"--time={args.time}", "run"
        ]
        proc_cpu = subprocess.Popen(prog)

        os.sched_setaffinity(proc_cpu.pid, {core})

        print("Starting Memory benchmark...")
        sysbench_mem_args = [arg for arg in other if arg.startswith("--memory") or "forced" in arg]
        prog = [
            "sysbench", "memory", *sysbench_mem_args, f"--threads={mem_threads}", f"--time={args.time}", "run"
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

    if args.data:
        print("Writing output...")
        sysbench_file_output_args = [arg for arg in other if arg.startswith("--file")]
        proc = subprocess.Popen(
            [
                "sysbench", "fileio", *sysbench_file_output_args, "--threads=1", "prepare"
            ]
        )
        proc.wait()

        for path in this_dir.glob("*test_file*"):
            path.rename(path.parent.joinpath(f"{args.out}"))


if __name__ == "__main__":
    main()
