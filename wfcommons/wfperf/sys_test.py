import pathlib
import argparse
import subprocess 
import os
import time

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Task Name")
    parser.add_argument("--time", type=int, help="maximum compute time in seconds")
    parser.add_argument("--percent-cpu", type=float, help="percent of threads which will be cpu heavy")
    parser.add_argument("--save", type=pathlib.Path, help="directory to save to.")

    return parser

def main():
    parser = get_parser()
    args, other = parser.parse_known_args()
    name = args.name
    save_dir = pathlib.Path(args.save)

    
    
    # save_dir = [item for item in other if "save" in item][0]
    # save_dir = pathlib.Path(save_dir.split("=")[1]) 

    
    with save_dir.joinpath(f"{name}_cpu.txt").open("w+") as fp_cpu, save_dir.joinpath(f"{name}_memory.txt").open("w+") as fp_mem, save_dir.joinpath(f"{name}_ps.txt").open("w+") as fp_ps:
        num_cores = 1 #os.cpu_count()
        cpu_threads = 1 #int(args.percent_cpu*10)
        mem_threads = 1 #int(10 - cpu_threads)


        print("Starting CPU benchmark...")
        sysbench_cpu_args = [arg for arg in other if arg.startswith("--cpu")] + [f"--threads={cpu_threads}"]
        
        for i in range(num_cores):

            proc_cpu = subprocess.Popen(
                [
                    "sysbench", "cpu",
                    *sysbench_cpu_args, "run"
                ], 
                stdout=fp_cpu, stderr=fp_cpu, 
            )
            
            pid_1 = proc_cpu.pid
            print(pid_1)
            os.sched_setaffinity(pid_1, {i})   
            

            print("Starting Memory benchmark...")
            sysbench_mem_args = [arg for arg in other if arg.startswith("--memory")] + [f"--time={args.time}", f"--threads={mem_threads}"]
            proc_mem = subprocess.Popen(
                [
                    "sysbench", "memory","run",
                    *sysbench_mem_args
                ], 
                stdout=fp_mem, stderr=fp_mem
            )
            
            pid_2 = proc_mem.pid
            print(pid_2)
            os.sched_setaffinity(pid_2, {i})

        proc = subprocess.Popen(["ps", "-o","pid,psr,comm,lstart"], stdout=fp_ps)
        proc.wait()
        proc_cpu.wait()
        print(time.time)
        subprocess.Popen(["killall", "sysbench"])
        
 
    
if __name__ == "__main__":
    main()