import argparse
from re import S
import subprocess 
import os
import subprocess
import time
from typing import Generator, List
import pandas as pd 
import pathlib 

this_dir = pathlib.Path(__file__).resolve().parent
DO_TEST = True 

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Task Name")
    parser.add_argument("--work")
    parser.add_argument("--real", action="store_true", help="if doing for real task.")
    parser.add_argument("--stress", action="store_true", help="set if want to evaluate with stress in the other cores.")

    return parser

def cpu_work(work: str, cpu_threads: int):
    procs_cpu = []
    print(f"------------------------- Calibrating ------------------------------")
    procs_cpu = []
    print(f"cpu_threads={cpu_threads}")

    for i in range(cpu_threads):    
        print(f"Starting CPU benchmark {i}...")
        prog = ["/home/tgcoleman/tests/cpu_benchmark", f"{work}", "1"]
        proc_cpu = subprocess.Popen(prog)
        os.sched_setaffinity(proc_cpu.pid, {15})
        procs_cpu.append(proc_cpu)
    

    return procs_cpu 



def benchmark_call(task: str, work: str, cpu_threads: int):
    print(f"------------------------- Starting Benchmark for {task} ------------------------------")
    procs_cpu = []
    mem_threads = int(10 - cpu_threads)
    # cpu_threads = str(cpu_threads)
    # mem_threads = str(mem_threads)
    print(f"cpu_threads={cpu_threads}, mem_threads={mem_threads}")

    for i in range(cpu_threads):    
        print(f"Starting CPU benchmark {i}...")
        prog = ["/home/tgcoleman/tests/cpu_benchmark", work, "1"]
        proc_cpu = subprocess.Popen(prog)
        os.sched_setaffinity(proc_cpu.pid, {15})
        procs_cpu.append(proc_cpu)
        

    print(f"Starting Memory benchmark...")
    prog = ["stress", "--vm", f"{mem_threads}"]     
    proc_mem = subprocess.Popen(prog)
    os.sched_setaffinity(proc_mem.pid, {15})


    return procs_cpu, proc_mem

def real():
    print("--------------------- Starting Real ----------------------------")
    prog = ["/home/tgcoleman/Montage/bin/mViewer", "-ct", "1", "-gray", 
            "/home/tgcoleman/montage-workflow-for-tina/data/1-mosaic.fits", 
            "-1s", "max", "gaussian", "-png", 
            "/home/tgcoleman/montage-workflow-for-tina/data/1-mosaic.png" ]
    proc_real = subprocess.Popen(prog)
    
    return proc_real

def stress_cores() -> List:
    procs = []
    for core in range(os.cpu_count()//2):
        print("Starting Memory stress test...")
        prog = ["stress", "-m", "1"]
        proc_mem = subprocess.Popen(prog)
        os.sched_setaffinity(proc_mem.pid, {core})
        procs.append(proc_mem)
    return procs

def main():
    parser = get_parser()
    args = parser.parse_args()

    if DO_TEST:
        print("Starting Calibration...")
        rows = []
        for cpu_threads in range(1, 10):
            start = time.time()
            procs_cpu = cpu_work(float(args.work), cpu_threads)
            for proc_cpu in procs_cpu:
                proc_cpu.wait()
            duration = time.time() - start
            rows.append([cpu_threads, duration])
        df = pd.DataFrame(rows, columns=["cpu_threads", "duration"])
        print(df)
    else:
        rows = []
        if args.real:
            start = time.time()
            proc_real = real()
            if args.stress:
                other_proc_mems = stress_cores()
            proc_real.wait()
            duration = time.time() - start
            rows.append(duration)
            if args.stress:
                for proc in other_proc_mems:
                    proc.kill()
                proc = subprocess.Popen(["killall", "stress"])
                proc.wait()
            df = pd.DataFrame(rows, columns=["duration"])
            if args.stress:
                df.to_csv(this_dir.joinpath(f"benchmark_{args.name}_real_stress.csv"))
            else:
                df.to_csv(this_dir.joinpath(f"benchmark_{args.name}_real.csv"))
        else:
            for cpu_threads in range(1, 10):
                start = time.time()
                procs_cpu, proc_mem = benchmark_call(args.name, args.work, cpu_threads)
                if args.stress:
                    other_proc_mems = stress_cores()
                for proc_cpu in procs_cpu:
                    proc_cpu.wait()
                duration = time.time() - start
                rows.append([cpu_threads, duration])
                if proc_mem is not None:
                    proc_mem.kill()  
                if args.stress:
                    for proc in other_proc_mems:
                        proc.kill()
                proc = subprocess.Popen(["killall", "stress"])
                proc.wait()
                df = pd.DataFrame(rows, columns=["cpu_threads", "duration"])
                if args.stress:
                    df.to_csv(this_dir.joinpath(f"benchmark_{args.name}_{args.work}_stress.csv"))
                else:
                    df.to_csv(this_dir.joinpath(f"benchmark_{args.name}_{args.work}.csv"))
            


    

if __name__ == "__main__":
    main()
