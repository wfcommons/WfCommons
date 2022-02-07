from wfcommons.wfchef.recipes import BlastRecipe
from wfcommons.wfperf.perf import WorkflowBenchmark
import pathlib
import glob
import os
import argparse

this_dir = pathlib.Path(__file__).parent.resolve()
# Rememeber to run as PYTHONUNBUFFERED=1 python test.py 

def cleanup_sys_files() -> None:
    """Remove files already used"""
    input_files = glob.glob("*input*.txt")
    output_files = glob.glob("*output.txt")
    all_files = input_files + output_files
    for t in all_files:
        os.remove(t)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datafootprint", action="store_true", help="Runs with datafootprint.")
    parser.add_argument("--input-file", action="store_true", help="Runs from input-file.")
    parser.add_argument("--clean", action="store_true", help="Cleans the files.")
   
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()

    bench = WorkflowBenchmark(BlastRecipe, 100)
    if args.datafootprint:
        bench.create_benchmark(
            save_dir=this_dir, 
            percent_cpu=0.5, 
            data=10
        )
        bench.run(pathlib.Path("/workspace/wfcommons/wfcommons/wfperf/Blast-Benchmark-100.json"),this_dir)

    
    if args.input_file:
        bench.generate_input_file(this_dir.joinpath("test.json"))
        bench.create_benchmark_from_input_file(
            this_dir, 
            this_dir.joinpath("test.json"), 
            lock_files_folder=this_dir.joinpath("tmp")
        )
    
        bench.run(pathlib.Path("/workspace/wfcommons/wfcommons/wfperf/Blast-Benchmark-100.json"),this_dir)
    
    if args.clean:
        cleanup_sys_files()

if __name__ == "__main__":
    main()  