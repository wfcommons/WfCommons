from wfcommons.wfchef.recipes import BlastRecipe
from wfcommons.wfperf.perf import WorkflowBenchmark
import pathlib

this_dir = pathlib.Path(__file__).parent.resolve()

def main():
    bench = WorkflowBenchmark(BlastRecipe, 300)
    # bench.create_benchmark(this_dir, 0.5, data_footprint=1000)
    # bench.generate_input_file(this_dir.joinpath("test.txt"))
    bench.create_benchmark_from_input_file(this_dir, this_dir.joinpath("test.txt"))

if __name__ == "__main__":
    main()