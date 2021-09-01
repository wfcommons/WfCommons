#!/usr/bin/python3

import pathlib
import argparse
import subprocess 
import os
from wfcommons.wfperf.lock import lock_core, unlock_core



this_dir = pathlib.Path(__file__).resolve().parent

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Task Name")
    parser.add_argument("--time", type=int, help="maximum compute time in seconds")
    parser.add_argument("--percent-cpu", type=float, help="percent of threads which will be cpu heavy")
    parser.add_argument("-l", "--lock", help="Path to lock file.")
    parser.add_argument("-n", "--num-cores", help="Path to cores file.")
    return parser 


def main():
    parser = get_parser()
    args, other = parser.parse_known_args()
    name = args.name
    assert(args.percent_cpu >= 0.0 and args.percent_cpu <= 1.0)

    save_dir = [item for item in other if "save" in item][0]
    save_dir = pathlib.Path(save_dir.split("=")[1]) 

    path_locked = args.lock
    path_cores = args.num_cores
    path_locked.write_text("")
    path_cores.write_text("")

    for path in save_dir.glob("*test_file*"):
        _, rest = path.name.split("test_file", 1)
        path.rename(path.parent.joinpath(f"test_file{rest}"))

    sysbench_file_input_args = [
        arg for arg in other 
        if arg.startswith("--file") and not arg.startswith("--file-num")
    ] + ["--threads=1"] 
    
    core = lock_core(path_locked, path_cores)

    with save_dir.joinpath(f"{name}_fileio_run.txt").open("w+") as fp:
        print("Starting IO benchmark...")
        proc = subprocess.Popen(["sysbench", "fileio", *sysbench_file_input_args, "run"], stdout=fp, stderr=fp)
        proc.wait()
    
    with save_dir.joinpath(f"{name}_fileio_cleanup.txt").open("w+") as fp:
        proc = subprocess.Popen(["sysbench", "fileio", *sysbench_file_input_args, "cleanup"], stdout=fp, stderr=fp)
        proc.wait()

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

    # with save_dir.joinpath(f"{name}_cpu.txt").open("w+") as fp_cpu, save_dir.joinpath(f"{name}_memory.txt").open("w+") as fp_mem, save_dir.joinpath(f"{name}_ps_log.txt").open("w+") as fp_ps:
    #     num_cores = os.cpu_count()
    #     cpu_threads = int(args.percent_cpu*10)
    #     mem_threads = int(10 - cpu_threads)
        
    #     print(cpu_threads, mem_threads)
        
    #     for i in range(num_cores):
    #         print("Starting CPU benchmark...")
    #         sysbench_cpu_args = [arg for arg in other if arg.startswith("--cpu")] + [f"--time={args.time}", f"--threads={cpu_threads}"]
    #         proc_cpu = subprocess.Popen(
    #             [
    #                 "sysbench", "cpu", "--forced-shutdown=0", #time", "timeout", f"{args.time}s",
    #                 *sysbench_cpu_args, "run"
    #             ], 
    #             stdout=fp_cpu, stderr=fp_cpu
    #         )
    #         os.sched_setaffinity(proc_cpu.pid, {i})
        
    #         print("Starting Memory benchmark...")     
    #         sysbench_mem_args = [arg for arg in other if arg.startswith("--memory")] + [f"--time={args.time}", f"--threads={mem_threads}"]
    #         proc_mem = subprocess.Popen(
    #             ["sysbench", "memory", "run",*sysbench_mem_args], #"time",  
    #             stdout=fp_mem, stderr=fp_mem
    #         )
    #         os.sched_setaffinity(proc_mem.pid, {i})

    #     proc = subprocess.Popen(["ps", "-o","pid,psr,comm,thcount"], stdout=fp_ps)
    #     proc.wait()    
    #     proc_cpu.wait()
    #     proc_mem.wait()

    with save_dir.joinpath(f"{name}_fileio_output.txt").open("w+") as fp:
        print("Writing output...")
        sysbench_file_output_args = [arg for arg in other if arg.startswith("--file")] + ["--threads=1"]
        proc = subprocess.Popen(["sysbench", "fileio", *sysbench_file_output_args, "prepare"], stdout=fp, stderr=fp)
        proc.wait()
    
    for path in this_dir.glob("*test_file*"):
        path.rename(path.parent.joinpath(f"{save_dir}/{args.name}_{path.name}"))

if __name__ == "__main__":
    main()


