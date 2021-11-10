import argparse
from pprint import pprint
from re import S
import subprocess 
import os
import subprocess
import time
from typing import  List, Tuple
import pandas as pd 
import pathlib 
from concurrent.futures import ThreadPoolExecutor

this_dir = pathlib.Path(__file__).resolve().parent

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Task Name")
    parser.add_argument("--real-time", type=float)
    parser.add_argument("--cpu-work-lower-bound", default=1, type=float)
    parser.add_argument("--mem-work", default=1, type=float)
    parser.add_argument("--cpu-threads", type=int)
    parser.add_argument("--real", action="store_true", help="if doing for real task.")
    parser.add_argument("--stress", action="store_true", help="set if want to evaluate with stress in the other cores.")

    return parser

def cpu_call(work: float, output: pathlib.Path) -> float:
    prog = ["time", "/home/tgcoleman/tests/cpu_benchmark", str(work), "1"]
    output.parent.mkdir(exist_ok=True, parents=True)
    with output.open("a+") as fp:
        start = time.time()
        proc = subprocess.Popen(prog, stdout=fp, stderr=fp)
        os.sched_setaffinity(proc.pid, {15})
        proc.wait()

        return time.time() - start

def mem_call(work: float, output: pathlib.Path) -> float:
    prog = ["time", "/home/tgcoleman/tests/mem_benchmark", str(work)]
    output.parent.mkdir(exist_ok=True, parents=True)
    with output.open("a+") as fp:
        start = time.time()
        proc = subprocess.Popen(prog, stdout=fp, stderr=fp)
        os.sched_setaffinity(proc.pid, {15})
        proc.wait()

        return time.time() - start

def benchmark_call(cpu_work: float,
                   mem_work: float,
                   cpu_threads: int,
                   mem_threads: int) -> Tuple[float, float, float]:
    """Calls benchmarking processes and returns compute times

    Args:
        cpu_work: Work for CPU processes
        mem_work: Work for memory processes
        cpu_threads: Number of threads for CPU
        mem_threads: Number of memory threads

    Returns:
        float: Average execution time per CPU Thread
        float: Average execution time per Memory Thread
        float: Total execution time to execute all threads
    """
    executor = ThreadPoolExecutor(max_workers=cpu_threads+mem_threads)

    logdir = this_dir.joinpath(
        "logs", 
        f"{cpu_threads}_threads", f"{mem_threads}_mem_threads", 
        f"{cpu_work}_cpu_work", f"{mem_work}_mem_work"
    )

    start = time.time() 
    cpu_futures = []
    for i in range(cpu_threads):
        cpu_futures.append(
            executor.submit(
                cpu_call, 
                cpu_work, 
                logdir.joinpath(f"cpu_thread_{i}.txt")
            )
        )

    mem_futures = []
    for i in range(mem_threads):
        mem_futures.append(
            executor.submit(
                mem_call, 
                mem_work, 
                logdir.joinpath(f"mem_thread_{i}.txt")
            )
        )


    cpu_time = sum([future.result() for future in cpu_futures]) / cpu_threads
    mem_time = sum([future.result() for future in mem_futures]) / mem_threads
    duration = time.time() - start 

    return cpu_time, mem_time, duration

def find_work(real_time: int, 
              cpu_threads: int,
              mem_threads: int,
              mem_work_ratio: float,
              cpu_work_lower_bound: float = 1) -> Tuple[float, float, float, float, float]:
    """Finds work that best estimates real time for a given mem to cpu work ratio
    
    Args:
        real_time: Time to find work for
        cpu_threads: Number of CPU threads
        mem_threads: Number of memory threads
        mem_work_ratio: Amount of work per unit of CPU work (mem_work := cpu_work * mem_work_ratio)
        cpu_work_lower_bound: Beginning lower bound for CPU work

    Returns:
        float: Average execution time per CPU Thread
        float: Average execution time per Memory Thread
        float: Total execution time to execute all threads
        float: cpu work
        float: memory work
    """
    lower_bound = cpu_work_lower_bound
    upper_bound = lower_bound * 2
    while True:
        *_, duration = benchmark_call(
            upper_bound,
            upper_bound * mem_work_ratio,
            cpu_threads,
            mem_threads
        )
        if duration >= real_time:
            break
        lower_bound = upper_bound
        upper_bound *= 2
    
    work = (float(upper_bound) + float(lower_bound)) / 2.0
    while max(abs(upper_bound - lower_bound), abs(upper_bound - lower_bound) * mem_work_ratio) > 1:
        *_, duration = benchmark_call(
            work,
            work * mem_work_ratio,
            cpu_threads,
            mem_threads
        )
        if duration < real_time:
            lower_bound = work 
        else:
            upper_bound = work 
        work = (float(upper_bound) + float(lower_bound)) / 2.0

    cpu_time, mem_time, duration = benchmark_call(
        work,
        work * mem_work_ratio,
        cpu_threads,
        mem_threads
    )
    return cpu_time, mem_time, duration, work, work * mem_work_ratio

def find_balanced_work(real_time: int, 
                       cpu_threads: int,
                       mem_threads: int):
    """Finds work that best estimates real time AND balances cpu and memory execution time
    
    Args:
        real_time: Time to find work for
        cpu_threads: Number of CPU threads
        mem_threads: Number of memory threads

    Returns:
        float: Average execution time per CPU Thread
        float: Average execution time per Memory Thread
        float: Total execution time to execute all threads
        float: cpu work
        float: memory work
    """
    lower_bound = 0
    upper_bound = 1 

    while True:
        cpu_time, mem_time, duration, *_ = find_work(real_time, cpu_threads, mem_threads, upper_bound)
        print(f"mem = cpu * {upper_bound:.2f}: {cpu_time:.2f}s CPU & {mem_time:.2f}s MEM -> {duration:.2f}s")
        if mem_time > cpu_time:
            break 
        lower_bound = upper_bound
        upper_bound *= 2


    ratio = (upper_bound + lower_bound) / 2
    while abs(upper_bound - lower_bound) > 0.05:
        cpu_time, mem_time, duration, *_ = find_work(real_time, cpu_threads, mem_threads, ratio)
        print(f"mem = cpu * {ratio:.2f}: {cpu_time:.2f}s CPU & {mem_time:.2f}s MEM -> {duration:.2f}s")
        if cpu_time < mem_time:
            upper_bound = ratio 
        else:
            lower_bound = ratio

        ratio = (upper_bound + lower_bound) / 2

    cpu_time, mem_time, duration, cpu_work, mem_work = find_work(real_time, cpu_threads, mem_threads, ratio)
    return cpu_time, mem_time, duration, cpu_work, mem_work

def main():
    cpu_time, mem_time, duration, cpu_work, mem_work = find_balanced_work(49, 1, 9)

    pprint(dict(
        cpu_time=cpu_time,
        mem_time=mem_time,
        duration=duration,
        cpu_work=cpu_work,
        mem_work=mem_work
    ))




if __name__ == "__main__":
    main()
