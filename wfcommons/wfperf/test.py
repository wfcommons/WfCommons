from wfcommons.wfchef.recipes import BlastRecipe
from wfcommons.wfperf.perf import WorkflowBenchmark
import pathlib
import glob
import os

this_dir = pathlib.Path(__file__).parent.resolve()
# Rememeber to run as PYTHONUNBUFFERED=1 python test.py 

def cleanup_sys_files() -> None:
    """Remove files already used"""
    input_files = glob.glob("*input*.txt")
    output_files = glob.glob("*output.txt")
    all_files = input_files + output_files
    for t in all_files:
        os.remove(t)
def main():
    bench = WorkflowBenchmark(BlastRecipe, 100)
    # bench.create_benchmark(this_dir, 0.5, input_data=1000)
    # bench.generate_input_file(this_dir.joinpath("test.json"))
    bench.create_benchmark_from_input_file(this_dir, this_dir.joinpath("test.json"), lock_files_folder=this_dir.joinpath("tmp"))
    bench.run(pathlib.Path("/workspace/wfcommons/wfcommons/wfperf/Blast-Benchmark-100.json"),this_dir)
    # cleanup_sys_files()

if __name__ == "__main__":
    main()  