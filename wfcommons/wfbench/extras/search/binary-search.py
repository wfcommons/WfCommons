import argparse
import os
import pathlib
import subprocess
import sys
import time

from concurrent.futures import ThreadPoolExecutor
from pprint import pprint
from typing import List, Tuple

this_dir = pathlib.Path(__file__).resolve().parent


program_real = {
                "individuals": ["/home/tgcoleman/1000genome-sequential/bin/individuals.py",
                                "/home/tgcoleman/1000genome-sequential/data/20130502/ALL.chr1.250000.vcf",
                                "1","1","1001", "3000"]
                }

def cpu_call(work: float, output: pathlib.Path) -> float:
    prog = ["./wfperf-benchmark", "--cpu-work", str(work)]
    output.parent.mkdir(exist_ok=True, parents=True)
    with output.open("a+") as fp:
        start = time.time()
        proc = subprocess.Popen(prog, stdout=fp, stderr=fp)
        os.sched_setaffinity(proc.pid, {os.cpu_count() - 1})
        proc.wait()
        return time.time() - start


def mem_call(work: float, output: pathlib.Path) -> float:
    prog = ["./wfperf-benchmark", "--mem-work", str(work)]
    output.parent.mkdir(exist_ok=True, parents=True)
    with output.open("a+") as fp:
        start = time.time()
        proc = subprocess.Popen(prog, stdout=fp, stderr=fp)
        os.sched_setaffinity(proc.pid, {os.cpu_count() - 1})
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
    executor = ThreadPoolExecutor(max_workers=cpu_threads + mem_threads)

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


def stress_test(cpu_work: float,
                mem_work: float,
                cpu_threads: int):
    """Stress remaining cores of the machines in which the benchmark is not running. 
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
    # print("Starting stress...")
    for core in range(os.cpu_count() // 2):
        prog = ["stress", "-m", "1"]
        proc_mem = subprocess.Popen(prog)
        os.sched_setaffinity(proc_mem.pid, {core})

    mem_threads = 10 - cpu_threads
    cpu_time, mem_time, duration = benchmark_call(cpu_work, mem_work, cpu_threads, mem_threads)

    subprocess.Popen(["killall", "stress"])

    return cpu_time, mem_time, duration


def real_stress_test():
    print("--------------------- Starting Real Stress ----------------------------")
    for core in range(os.cpu_count() // 2):
        prog = ["stress", "-m", "1"]
        proc_mem = subprocess.Popen(prog)
        os.sched_setaffinity(proc_mem.pid, {core})
    start = time.time()
    prog = ["/home/cc/Montage/bin/mViewer", "-ct", "1", "-gray",
            "/home/cc/montage-workflow-for-tina/data/1-mosaic.fits",
            "-1s", "max", "gaussian", "-png",
            "/home/cc/montage-workflow-for-tina/data/1-mosaic.png"]
    proc_real = subprocess.Popen(prog)
    os.sched_setaffinity(proc_real.pid, {os.cpu_count() - 1})
    proc_real.wait()
    duration = time.time() - start
    subprocess.Popen(["killall", "stress"])

    print(duration)


def find_balanced_work(real_time: float,
                       cpu_threads: int,
                       mem_threads: int,
                       cpu_work: int = 1000,
                       mem_work: int = 10000,
                       threshold: float = 0.05):
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
    counter = 0
    best_work = (cpu_work, mem_work, 99999)

    while True:
        counter += 1
        # print(f"[{counter}] CPU: {cpu_work} - MEM: {mem_work}")
        cpu_time, mem_time, duration = benchmark_call(cpu_work, mem_work, cpu_threads, mem_threads)
        if real_time == 0:
            return cpu_time, mem_time, duration, cpu_work, mem_work

        max_ratio = max(abs(1 - real_time / mem_time), abs(1 - real_time / cpu_time))
        # print(f"  CPU Time: {cpu_time:.2f}")
        # print(f"  MEM Time: {mem_time:.2f}")
        # print(f"  Max Ratio: {max_ratio:.2f}")

        # times difference within threshold
        if threshold >= max_ratio:
            return cpu_time, mem_time, duration, cpu_work, mem_work

        # if threshold < best_work[2]:
        #     best_work = (cpu_work, mem_work, threshold)
        #
        # else:
        cpu_work = round(cpu_work * (real_time / cpu_time))
        mem_work = round(mem_work * (real_time / mem_time))

        if counter == 100:
            return cpu_time, mem_time, duration, cpu_work, mem_work
        print(".", end=" ")
        # print("")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--real-time", type=float, default=49.0)
    parser.add_argument("--real", action="store_true")
    parser.add_argument("--cpu-threads", type=int, default=1)
    parser.add_argument("--cpu-work", type=int, default=1000)
    parser.add_argument("--mem-work", type=int, default=100000)
    parser.add_argument("--threshold", type=float, default=0.05)
    parser.add_argument("--stress", action="store_true")

    args, other = parser.parse_known_args()

    if args.real:
        real_stress_test()
    else:
        if args.stress:
            cpu_work = args.cpu_work
            mem_work = args.mem_work
            cpu_time, mem_time, duration = stress_test(cpu_work, mem_work, args.cpu_threads)
        else:
            cpu_time, mem_time, duration, cpu_work, mem_work = find_balanced_work(real_time=args.real_time,
                                                                                  cpu_threads=args.cpu_threads,
                                                                                  mem_threads=10 - args.cpu_threads,
                                                                                  cpu_work=args.cpu_work,
                                                                                  mem_work=args.mem_work,
                                                                                  threshold=args.threshold)
        # print("")
        # pprint(dict(
        #     cpu_time=cpu_time,
        #     mem_time=mem_time,
        #     duration=duration,
        #     cpu_work=cpu_work,
        #     mem_work=mem_work
        # ))
        # print(f"{args.cpu_threads}-{10 - args.cpu_threads},{duration},{cpu_work},{mem_work}")
        print(f"--cpu-threads {args.cpu_threads} --cpu-work {cpu_work} --mem-work {mem_work} # {duration}")


if __name__ == "__main__":
    main()
