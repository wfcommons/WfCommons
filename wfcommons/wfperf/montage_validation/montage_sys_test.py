import pathlib
import argparse
import subprocess 
import os
from typing import List
from wfcommons.wfperf.montage_validation.lock import file_lock, file_unlock


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Task Name")
    parser.add_argument("--save", type=pathlib.Path, help="directory to save to.")


    return parser

def set_affinity(lockfile_path: pathlib.Path, 
                 corefile_path: pathlib.Path, 
                 pid: int,
                 core: int) -> None:
        file_lock(lockfile_path, corefile_path, core)
        os.sched_setaffinity(pid, {core})   
        file_unlock(lockfile_path, corefile_path, core)

def main():
    parser = get_parser()
    args, other = parser.parse_known_args()
    name = args.name
    save_dir = pathlib.Path(args.save)

    path_locked = "/home/tgcoleman/tests/Montage/cores.txt.lock"
    path_cores = "/home/tgcoleman/tests/Montage/cores.txt"
    
    with save_dir.joinpath(f"{name}_cpu.txt").open("w+") as fp_cpu, save_dir.joinpath(f"{name}_memory.txt").open("w+") as fp_mem:
        num_cores = os.cpu_count()
        sysbench_cpu_args = [arg for arg in other if "cpu" in arg or "time" in arg]
        percent_cpu = [arg for arg in sysbench_cpu_args if arg.startswith("--percent")][0]
        cpu_threads = int(float(percent_cpu.split("=")[1])*10)
        mem_threads = int(10 - cpu_threads)
        print(f"cpu_threads={cpu_threads}, mem_threads={mem_threads}")

        print(sysbench_cpu_args)

        print("Starting CPU benchmark...")
        
        proc_cpus: List[subprocess.Popen] = []
        proc_mems: List[subprocess.Popen] = []
        
        for i in range(num_cores):
            time =[arg for arg in sysbench_cpu_args if arg.startswith("--time")][0]
            time = int(time.split("=")[1])

            #if time is relevant do cpu and memory, else only cpu
            if time > 100:
                proc_cpus.append(subprocess.Popen(
                    [
                        "sysbench", "cpu",
                        *sysbench_cpu_args, f"--threads={cpu_threads}", "run"
                    ], 
                    stdout=fp_cpu, stderr=fp_cpu, 
                ))
                
                set_affinity(path_locked, path_cores, proc_cpus[-1].pid, i)

                print("Starting Memory benchmark...")
                sysbench_mem_args = [arg for arg in other if arg.startswith("--memory") or "time" in arg]
                proc_mems.append(subprocess.Popen(
                    [
                        "sysbench", "memory","run",
                        *sysbench_mem_args, f"--threads={mem_threads}"
                    ], 
                    stdout=fp_mem, stderr=fp_mem
                ))
                
                set_affinity(path_locked, path_cores, proc_mems[-1].pid, i)
            else:
                proc_cpus.append(subprocess.Popen(
                    [
                        "sysbench", "cpu",
                        *sysbench_cpu_args, f"--threads={cpu_threads}", "run"
                    ], 
                    stdout=fp_cpu, stderr=fp_cpu, 
                ))
                
                set_affinity(path_locked, path_cores, proc_cpus[-1].pid, i)


        for proc_cpu in proc_cpus:
            proc_cpu.wait()
        for proc_mem in proc_mems:
            proc_mem.kill()
 
    
if __name__ == "__main__":
    main()
