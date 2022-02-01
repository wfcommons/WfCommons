from wfcommons.wfchef.recipes import BlastRecipe
from wfcommons.wfperf.perf import WorkflowBenchmark
import pathlib

this_dir = pathlib.Path(__file__).parent.resolve()

def main():
    bench = WorkflowBenchmark(BlastRecipe, 100)
    # bench.create_benchmark(this_dir, 0.5, data_footprint=1000)
    # bench.generate_input_file(this_dir.joinpath("test.json"))
    # bench.create_benchmark_from_input_file(this_dir, this_dir.joinpath("test.json"))
    # bench.create_benchmark(save_dir=pathlib.Path(this_dir),
    #                         percent_cpu=0.5,
    #                         cpu_work=150,
    #                         data_footprint=1000,
    #                         lock_files_folder=pathlib.Path("/tmp"))
    bench.create_benchmark_from_input_file(this_dir, this_dir.joinpath("test.json"), lock_files_folder=this_dir.joinpath("tmp"))
    bench.run(pathlib.Path("/workspace/wfcommons/wfcommons/wfperf/Blast-Benchmark-100.json"),this_dir)

if __name__ == "__main__":
    main()  