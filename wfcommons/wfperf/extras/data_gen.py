import subprocess
from typing import List, Union
import glob
import os



def generate_sys_data(num_files: int, file_total_size: int, test_mode: str = "seqwr"):
    
    params = [
            f"--file-num={num_files}",
            f"--file-total-size={file_total_size}G",
            f"--file-test-mode={test_mode}"
        ]
    print(" ".join(["sysbench","fileio", *params,"prepare"]))
    proc = subprocess.Popen(["sysbench","fileio", *params,"prepare"], stdout=subprocess.PIPE)
    proc.wait()
    out, _ = proc.communicate()
    if not out:
        raise FileNotFoundError("Couldn't create files. Check parameters.")
    
    # have to change the name back to test_file on bash
    for t in glob.glob("test_file.*"):
        os.rename(t, f"sys_{t}")
        


def cleanup_sys_files():
    for t in glob.glob("*test_file.*"):
        new = t.split("_")[1:]
        new = "_".join(new)
        os.rename(t, new)
    proc = subprocess.Popen(["sysbench", "fileio", "cleanup"], stdout=subprocess.PIPE)
    proc.wait()
    out, _ = proc.communicate()
    if not out:
        raise FileNotFoundError("Couldn't delete files.")