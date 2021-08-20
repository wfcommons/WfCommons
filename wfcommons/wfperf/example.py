from wfcommons.wfperf.perf import WorkflowBenchmark
from wfcommons.wfchef.recipes import BlastRecipe, MontageRecipe
import pathlib

this_dir = pathlib.Path(__file__).resolve().parent

def main():
    bench = WorkflowBenchmark(BlastRecipe, 100)
    bench.create(this_dir.joinpath("Blast"), 4, verbose=True)

if __name__=="__main__":
    main()