import pathlib
from typing import Union
from filelock import FileLock

def file_lock(path_locked: Union[str, pathlib.Path], path_cores: Union[str, pathlib.Path], core: int) -> None:
    
    with FileLock(path_locked) as lock:        
        lock.acquire()
        print("Lock acquired.")
        
        try:
            cores = {int(line) for line in path_cores.read_text().splitlines()}
            if core in cores:
                print(f"Sorry core {core} is taken.")
            else:
                cores.add(core) 
                print(f"Core {core} is yours!")

            path_cores.write_text("\n".join(map(str, cores)))

        finally:
            lock.release()
            print("Lock released.")

def file_unlock(path_locked, path_cores, core):
    with FileLock(path_locked) as lock: 
        lock.acquire()
        try:
            cores = {int(line) for line in path_cores.read_text().splitlines()}
            if core in cores:
                cores.remove(core) 
            else:
                print(f"This core is not being used.")

            path_cores.write_text("\n".join(map(str, cores)))

        finally:
            lock.release()

