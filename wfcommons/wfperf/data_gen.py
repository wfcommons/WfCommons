import subprocess


def generate_sys_data(num_files: int, file_total_size: int):
    
    params = [
            f"--file-num={num_files}",
            f"--file-total-size={file_total_size}",
        ]
    proc = subprocess.Popen(["sysbench fileio", param for param in params], stdout=subprocess.PIPE)
    out, _ = proc.communicate()
    if not out:
        raise FileNotFoundError("Couldn't create files. Check parameters.")

def cleanup_sys_files():

    proc = subprocess.Popen("sysbench fileio cleanup", stdout=subprocess.PIPE)
    out, _ = proc.communicate()
    if not out:
        raise FileNotFoundError("Couldn't delete files.")
