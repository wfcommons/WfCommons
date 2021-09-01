import pathlib
from filelock import FileLock
import os
import time 

def lock_core(path_locked: pathlib.Path, path_cores: pathlib.Path) -> int:
    all_cores = set(range(os.cpu_count()))
    available = set()
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

