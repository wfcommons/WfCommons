from io import StringIO
import pathlib 
import subprocess
import argparse
import pandas as pd

this_dir = pathlib.Path(__file__).resolve().parent

def gpu_benchmark(path, work, device):
    gpu_prog = [f"CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES={device} {path.joinpath('pi-cuda')} {work}"]
    print(gpu_prog)
    gpu_proc = subprocess.Popen(gpu_prog, shell=True)    

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("work", help="Amount of work for the gpu.")
    return parser

def get_available_gpus():
    proc = subprocess.Popen(["nvidia-smi", "--query-gpu=utilization.gpu", "--format=csv"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, _ = proc.communicate()
    df = pd.read_csv(StringIO(stdout.decode("utf-8")), sep=" ")
    return df[df["utilization.gpu"] <= 5].index.to_list()

def main(): 
    parser = get_parser()
    args = parser.parse_args()
    work = args.work

    available_gpus = get_available_gpus()
    print(available_gpus)

    if not available_gpus:
        print("No GPU available")
    else:
        device = available_gpus[0]
        print(f"Running on GPU {device}")
        gpu_benchmark(this_dir, work, device)
    # print(torch.cuda.is_available())    


if __name__=="__main__":
    main()