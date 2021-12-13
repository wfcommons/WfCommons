from wfcommons.wfperf.perf import WorkflowBenchmark
from wfcommons.wfchef.recipes import BlastRecipe
import pathlib
import logging
from wfcommons.wfperf.translator import PegasusTranslator

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

this_dir = pathlib.Path(__file__).resolve().parent

def main():
    num_tasks = 200
    bench = WorkflowBenchmark(BlastRecipe, num_tasks, logger=logger)
    
    benchmark_inputs = this_dir.joinpath("benchmark_input.json")
    bench.generate_input_file(benchmark_inputs)

    input(f"Fill out the benchmark input parameter file {benchmark_inputs.resolve()} and press ENTER when finished:")

    json_path = bench.create_benchmark_from_input_file(this_dir, benchmark_inputs)
    bench.run(json_path, this_dir)






if __name__ == "__main__":
    main()
