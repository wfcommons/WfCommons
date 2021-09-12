from wfcommons.wfperf.perf import WorkflowBenchmark
from wfcommons.wfchef.recipes import BlastRecipe
import pathlib
import argparse
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

this_dir = pathlib.Path(__file__).resolve().parent


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", help="Path to JSON")
    parser.add_argument("-c", "--create", action="store_true", help="Generate Workflow Benchmark when set.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Prints status information when set to true.")
    parser.add_argument("-s", "--save", help="Path to save directory.")
    parser.add_argument("-t", "--num-tasks", help="Number os tasks when create is true.")

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    num_tasks = int(args.num_tasks)
    save_dir = pathlib.Path(args.save)

    print("Running")

    bench = WorkflowBenchmark(BlastRecipe, num_tasks, logger=logger)

    if args.create:
        if args.verbose:
            print("Creating Recipe...")
        json_path = bench.create(save_dir, percent_cpu=0.5, data_footprint=10000)

    else:
        json_path = bench.create(save_dir, create=False, path=pathlib.Path(args.path))

    bench.run(json_path, save_dir)


if __name__ == "__main__":
    main()
